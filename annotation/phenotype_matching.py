import logging
import multiprocessing as mp
import time
from typing import Optional, Iterable

import nltk
from django.conf import settings
from django.db import connections

from annotation.models.models_phenotype_match import TextPhenotypeMatch, PhenotypeDescription, TextPhenotype, \
    TextPhenotypeSentence, filter_ambiguous_acronym_matches
from annotation.phenotype_matcher import PhenotypeMatcher, SkipAllPhenotypeMatchException
from annotation.phenotype_tokenizer import PhenotypeTokenizer
from patients.models import Patient

MAX_COMBO_LENGTH = 14  # Checked HPO words in DB


def _get_word_combos_and_spans(words_and_spans: list, max_combo_length: Optional[int]) -> list:
    combos_and_spans = []
    for i in range(len(words_and_spans)):
        combo_length = len(words_and_spans) - i
        if max_combo_length:
            combo_length = min(max_combo_length, combo_length)
        for j in range(i, i + combo_length):
            combo_and_span = words_and_spans[i:j + 1]
            combos_and_spans.append(combo_and_span)

    return combos_and_spans


def _get_word_combos_and_spans_sorted_by_length(words_and_spans, max_combo_length: Optional[int] = None) -> Iterable:
    word_combos_and_spans = _get_word_combos_and_spans(words_and_spans, max_combo_length)
    return reversed(sorted(word_combos_and_spans, key=lambda item: sum(len(i[0]) for i in item)))


def _get_terms_from_words(text_phenotype, words_and_spans_subset, phenotype_matcher: PhenotypeMatcher):
    """ Returns unsaved TextPhenotypeMatch instances; caller bulk_creates them. """
    offset_start = words_and_spans_subset[0][1][0]
    offset_end = words_and_spans_subset[-1][1][1]

    results = []
    for ontology_term_id in phenotype_matcher.get_matches(words_and_spans_subset):
        results.append(TextPhenotypeMatch(text_phenotype=text_phenotype,
                                          ontology_term_id=ontology_term_id,
                                          offset_start=offset_start,
                                          offset_end=offset_end))

    return results


def _sub_array_index(array, sub_array):
    for i in range(0, len(array) - len(sub_array) + 1):
        if array[i:i + len(sub_array)] == sub_array:
            return i
    return None


def parse_words(text_phenotype, input_words_and_spans, phenotype_matcher) -> list:
    results = []

    word_combos_and_spans = _get_word_combos_and_spans_sorted_by_length(input_words_and_spans, max_combo_length=MAX_COMBO_LENGTH)
    for words_and_spans_subset in word_combos_and_spans:
        words_results = _get_terms_from_words(text_phenotype, words_and_spans_subset, phenotype_matcher)
        if words_results:
            results.extend(words_results)
            if words_and_spans_subset != input_words_and_spans:  # More to match
                i = _sub_array_index(input_words_and_spans, words_and_spans_subset)
                before_words = input_words_and_spans[:i]
                after_words = input_words_and_spans[i + len(words_and_spans_subset):]

                if before_words:
                    before_results = parse_words(text_phenotype, before_words, phenotype_matcher)
                    results.extend(before_results)

                if after_words:
                    after_results = parse_words(text_phenotype, after_words, phenotype_matcher)
                    results.extend(after_results)

            break

    return results


def _split_adj_noun_and_noun(words, tags, spans):
    """ This is used to transform e.g. "long fingers and toes" into "long fingers and long toes"
        I looked and this pattern doesn't appear in HPO terms so is ok to split up """

    new_words_and_spans = []
    if tags == ['JJ', 'NNS', 'CC', 'NNS'] and words[2].lower() == 'and':
        new_order = [0, 1, 2, 0, 3]
        for i in new_order:
            new_words_and_spans.append((words[i], spans[i]))

    return new_words_and_spans


def _transform_words(words, tags, spans):
    # TODO:
    # X disorder => Abnormality of X
    # X Ca => X Cancer?

    TRANSFORMATIONS = [_split_adj_noun_and_noun]

    for t in TRANSFORMATIONS:
        new_words_and_spans = t(words, tags, spans)
        if new_words_and_spans:
            return new_words_and_spans

    return list(zip(words, spans))


def _process_text_phenotype(text_phenotype, phenotype_tokenizer, phenotype_matcher):
    tokenized_text_and_spans = phenotype_tokenizer.word_tokenise_and_spans(text_phenotype.text)

    tokenized_text = []
    tokenized_spans = []
    for ts in tokenized_text_and_spans:
        tokenized_text.append(ts[0])
        tokenized_spans.append(ts[1:])

    tagged = nltk.pos_tag(tokenized_text)

    words = []
    tags = []
    spans = []
    commented_out = False
    for (w, t), span in zip(tagged, tokenized_spans):
        if not commented_out:
            if w == t:
                if w.startswith('--'):
                    commented_out = True
            else:
                words.append(w)
                tags.append(t)
                offset_start = span[0]
                offset_end = span[1]
                spans.append((offset_start, offset_end))

    words_and_spans = _transform_words(words, tags, spans)

    try:
        matches = parse_words(text_phenotype, words_and_spans, phenotype_matcher)
        matches = filter_ambiguous_acronym_matches(matches)
        if matches:
            TextPhenotypeMatch.objects.bulk_create(matches)
    except SkipAllPhenotypeMatchException:
        logging.info("Completely skipping: %s", text_phenotype.text)

    text_phenotype.processed = True
    text_phenotype.save()


def _replace_comments_with_spaces(text):
    """ Replace with spaces instead of removing so that characters offsets remain the same """
    COMMENT = '-'
    COMMENT_REPLACE = ' '
    cleaned_chars = []

    last_char = None
    in_comment = False
    for c in text:
        if c == COMMENT and last_char == COMMENT:
            in_comment = True
        elif c == '\n':
            in_comment = False

        if in_comment:
            cleaned_chars.append(COMMENT_REPLACE)
        else:
            cleaned_chars.append(c)
        last_char = c

    return ''.join(cleaned_chars)


def create_phenotype_description(text, phenotype_matcher=None, defer_processing=False):
    """ PhenotypeMatcher takes a while to initialize, so can create once and pass in for
        faster bulk processing.
        If defer_processing is True the PhenotypeDescription/TextPhenotype rows are created but NLP matching is
        skipped (used by bulk_patient_phenotype_matching to batch the heavy work for parallel execution). """
    if not defer_processing and phenotype_matcher is None:
        phenotype_matcher = PhenotypeMatcher()
    phenotype_tokenizer = PhenotypeTokenizer()
    phenotype_description = PhenotypeDescription.objects.create(original_text=text)
    known_sentences = []
    unknown_sentences = []

    # Strip carriage returns - some browsers send \r\n in forms, Ajax usually sends \n
    # These need to be the same so save/ajax preview look the same
    text = text.replace('\r', '')
    text = _replace_comments_with_spaces(text)  # Need to keep length the same for offsets
    for sentence, sentence_offset in phenotype_tokenizer.sentences_and_offsets(text):
        text_phenotype, tp_created = TextPhenotype.objects.get_or_create(text=sentence)
        tps = TextPhenotypeSentence.objects.create(phenotype_description=phenotype_description,
                                                   text_phenotype=text_phenotype,
                                                   sentence_offset=sentence_offset)
        if tp_created:
            has_alpha_numeric = any(x.isalnum() for x in sentence)
            if not has_alpha_numeric:  # Impossible to match anything so skip
                text_phenotype.processed = True
                text_phenotype.save()
                known_sentences.append(tps)
            else:
                unknown_sentences.append(tps)
        else:
            known_sentences.append(tps)

    if not defer_processing and unknown_sentences:
        for sentence in unknown_sentences:
            try:
                _process_text_phenotype(sentence.text_phenotype, phenotype_tokenizer, phenotype_matcher)
            except Exception:
                logging.error("Problem processing phenotype for '%s'", sentence.text_phenotype.text)
                raise
            known_sentences.append(sentence)

    return phenotype_description


_WORKER_PHENOTYPE_TOKENIZER: Optional[PhenotypeTokenizer] = None
_WORKER_PHENOTYPE_MATCHER: Optional[PhenotypeMatcher] = None


def _phenotype_worker_init():
    """Pool initializer: drop inherited DB connections, warm NLTK models, build a per-worker PhenotypeMatcher."""
    connections.close_all()
    # Warm NLTK so the averaged-perceptron tagger + tokenizer aren't loaded
    # lazily on the first sentence each worker processes.
    nltk.pos_tag(nltk.word_tokenize("warm up"))
    global _WORKER_PHENOTYPE_TOKENIZER, _WORKER_PHENOTYPE_MATCHER
    _WORKER_PHENOTYPE_TOKENIZER = PhenotypeTokenizer()
    _WORKER_PHENOTYPE_MATCHER = PhenotypeMatcher()


def _process_text_phenotype_by_pk(text_phenotype_pk):
    text_phenotype = TextPhenotype.objects.get(pk=text_phenotype_pk)
    _process_text_phenotype(text_phenotype, _WORKER_PHENOTYPE_TOKENIZER, _WORKER_PHENOTYPE_MATCHER)


def bulk_patient_phenotype_matching(patients=None, cores=1):
    if patients is None:
        patients = Patient.objects.filter(phenotype__isnull=False).exclude(phenotype='')
        exclude_string = getattr(settings, "PATIENT_PHENOTYPE_EXCLUDE_STRING", None)
        if exclude_string:
            patients = patients.exclude(phenotype__contains=exclude_string)

    patients = list(patients)
    num_patients = len(patients)
    if not num_patients:
        logging.info("No patients")
        return

    # Phase 1: walk patients and register PhenotypeDescription/TextPhenotype/TextPhenotypeSentence rows.
    # No NLP runs here, so this stays single-threaded and avoids races on TextPhenotype.get_or_create.
    start = time.time()
    num_parsed_phenotypes = 0
    for i, p in enumerate(patients):
        parsed_phenotypes = p.process_phenotype_if_changed(defer_processing=True)
        num_parsed_phenotypes += parsed_phenotypes
        if not i % 50:
            perc_complete = 100.0 * i / num_patients
            logging.info("Registering patient sentences: %d/%d (%0.2f%%)", i, num_patients, perc_complete)
    logging.info("bulk_patient_phenotype_matching: register sentences took %.2f secs", time.time() - start)

    # Phase 2: NLP-process each unique unprocessed sentence exactly once, optionally in parallel.
    unprocessed_pks = list(TextPhenotype.objects.filter(processed=False).values_list("pk", flat=True))
    num_sentences = len(unprocessed_pks)
    if not num_sentences:
        logging.info("No unprocessed sentences to match")
        return

    logging.info("Matching %d unique sentences across %d core(s)", num_sentences, cores)
    start = time.time()
    if cores > 1:
        connections.close_all()  # ensure no live parent connection is inherited across fork
        ctx = mp.get_context("fork")
        with ctx.Pool(cores, initializer=_phenotype_worker_init) as pool:
            for done, _ in enumerate(
                pool.imap_unordered(_process_text_phenotype_by_pk, unprocessed_pks, chunksize=10),
                start=1,
            ):
                if not done % 50:
                    perc_complete = 100.0 * done / num_sentences
                    logging.info("Sentences processed: %d/%d (%0.2f%%)", done, num_sentences, perc_complete)
    else:
        phenotype_tokenizer = PhenotypeTokenizer()
        phenotype_matcher = PhenotypeMatcher()
        for j, pk in enumerate(unprocessed_pks):
            text_phenotype = TextPhenotype.objects.get(pk=pk)
            _process_text_phenotype(text_phenotype, phenotype_tokenizer, phenotype_matcher)
            if not j % 50:
                perc_complete = 100.0 * j / num_sentences
                logging.info("Sentences processed: %d/%d (%0.2f%%)", j, num_sentences, perc_complete)

    ts = time.time() - start
    logging.info("bulk_patient_phenotype_matching: NLP took %.2f secs", ts)
    time_per_sentence = ts / num_sentences
    sps = 1.0 / time_per_sentence if time_per_sentence else 0
    logging.info("%d sentences matched - %.2f per second/%.4f seconds per sentence",
                 num_sentences, sps, time_per_sentence)
