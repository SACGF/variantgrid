import logging
import time
from typing import List, Optional, Iterable

import nltk
from django.db.models import Count

from annotation.models.models_phenotype_match import PhenotypeMatchTypes, \
    TextPhenotypeMatch, PhenotypeDescription, TextPhenotype, TextPhenotypeSentence
from annotation.phenotype_matcher import PhenotypeMatcher, SkipAllPhenotypeMatchException
from library.utils import get_and_log_time_since, invert_dict_of_lists
from patients.models import Patient

MAX_COMBO_LENGTH = 14  # Checked HPO words in DB


def get_word_combos_and_spans(words_and_spans: List, max_combo_length: Optional[int]) -> List:
    combos_and_spans = []
    for i in range(len(words_and_spans)):
        combo_length = len(words_and_spans) - i
        if max_combo_length:
            combo_length = min(max_combo_length, combo_length)
        for j in range(i, i + combo_length):
            combo_and_span = words_and_spans[i:j + 1]
            combos_and_spans.append(combo_and_span)

    return combos_and_spans


def get_word_combos_and_spans_sorted_by_length(words_and_spans, max_combo_length: Optional[int] = None) -> Iterable:
    word_combos_and_spans = get_word_combos_and_spans(words_and_spans, max_combo_length)
    return reversed(sorted(word_combos_and_spans, key=lambda item: sum([len(i[0]) for i in item])))


def get_terms_from_words(text_phenotype, words_and_spans_subset, phenotype_matcher: PhenotypeMatcher):
    hpo_list, omim_list, gene_symbols = phenotype_matcher.get_matches(words_and_spans_subset)

    offset_start = words_and_spans_subset[0][1][0]
    offset_end = words_and_spans_subset[-1][1][1]

    results = []
    for ontology_term_id in hpo_list:
        tpm = TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype,
                                                match_type=PhenotypeMatchTypes.HPO,
                                                ontology_term_id=ontology_term_id,
                                                offset_start=offset_start,
                                                offset_end=offset_end)
        results.append(tpm)

    for ontology_term_id in omim_list:
        tpm = TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype,
                                                match_type=PhenotypeMatchTypes.OMIM,
                                                ontology_term_id=ontology_term_id,
                                                offset_start=offset_start,
                                                offset_end=offset_end)
        results.append(tpm)

    for gene_symbol_id in gene_symbols:
        tpm = TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype,
                                                match_type=PhenotypeMatchTypes.GENE,
                                                gene_symbol_id=gene_symbol_id,
                                                offset_start=offset_start,
                                                offset_end=offset_end)
        results.append(tpm)

    return results


def sub_array_index(array, sub_array):
    for i in range(0, len(array) - len(sub_array) + 1):
        if array[i:i + len(sub_array)] == sub_array:
            return i
    return None


def parse_words(text_phenotype, input_words_and_spans, phenotype_matcher) -> List:
    results = []

    word_combos_and_spans = get_word_combos_and_spans_sorted_by_length(input_words_and_spans, max_combo_length=MAX_COMBO_LENGTH)
    for words_and_spans_subset in word_combos_and_spans:
        words_results = get_terms_from_words(text_phenotype, words_and_spans_subset, phenotype_matcher)
        if words_results:
            results.extend(words_results)
            if words_and_spans_subset != input_words_and_spans:  # More to match
                i = sub_array_index(input_words_and_spans, words_and_spans_subset)
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


def split_adj_noun_and_noun(words, tags, spans):
    """ This is used to transform eg "long fingers and toes" into "long fingers and long toes"
        I looked and this pattern doesn't appear in HPO terms so is ok to split up """

    new_words_and_spans = []
    if tags == ['JJ', 'NNS', 'CC', 'NNS'] and words[2].lower() == 'and':
        new_order = [0, 1, 2, 0, 3]
        for i in new_order:
            new_words_and_spans.append((words[i], spans[i]))

    return new_words_and_spans


def transform_words(words, tags, spans):
    # TODO:
    # X disorder => Abnormality of X
    # X Ca => X Cancer?

    TRANSFORMATIONS = [split_adj_noun_and_noun]

    for t in TRANSFORMATIONS:
        new_words_and_spans = t(words, tags, spans)
        if new_words_and_spans:
            return new_words_and_spans

    return list(zip(words, spans))


def tokens_and_spans(txt, tokens):
    offset = 0

    tokenized_text_and_spans = []
    for token in tokens:
        offset = txt.find(token, offset)
        tokenized_text_and_spans.append((token, offset, offset + len(token)))
        offset += len(token)

    return tokenized_text_and_spans


def handle_unsplit(unsplit_digits_and_spans, split_words_and_spans):
    if unsplit_digits_and_spans:
        unsplit_word = '/'.join([d[0] for d in unsplit_digits_and_spans])
        my_offset_start = unsplit_digits_and_spans[0][1]
        my_offset_end = unsplit_digits_and_spans[-1][2]
        split_words_and_spans.append((unsplit_word, my_offset_start, my_offset_end))
        unsplit_digits_and_spans[:] = []


def initial_replace_characters(text):
    REPLACE_CHARS = {';': ' '}

    cleaned_text = text
    for from_char, to_char in REPLACE_CHARS.items():
        cleaned_text = cleaned_text.replace(from_char, to_char)

    return cleaned_text


def word_tokenise_and_spans(txt):
    """ Split by nltk.sent_tokenize and / (only if alphanumeric not for dates) """

    cleaned_text = initial_replace_characters(txt)
    tokens = nltk.word_tokenize(cleaned_text)
    split_words_and_spans = []

    # {TO: FROM} - These have to be the same size! Then switch them out later
    PRE_SPLIT_REPLACE = {"dv/p": "dv_p"}

    # replace word : List of alternate words to replace
    REPLACE_ALTERNATES = {
        "bartter": ["bartters", "bartter's"],
        "bull's eye": ["bullseye"],
        "developmental": ["dev", "devt", "dv_p", "devel"],
        "dysmorphic": ["dysm"],
        "gaucher": ["gauchers", "gaucher's"],
        "hemophagocytic": ["haemophagcytic", "haemophogocytic", "haemophagocytic"],
        "leukemia": ["leuk"],
        "seizures": ["sz"],
        "sanfilippo": ["sanfilipo"],
        "sclerosis": ["schlerosis", "scleosis"],
        "respiratory": ["rep", "resp"],
        "syndrome": ["synd", "syndromes"],
        "undescended": ["undesc"],
    }
    replace_words = invert_dict_of_lists(REPLACE_ALTERNATES)

    for word, offset_start, _ in tokens_and_spans(cleaned_text, tokens):
        new_offset_start = offset_start
        unsplit_digits_and_spans = []  # Don't split up dates

        word = PRE_SPLIT_REPLACE.get(word, word)
        for w in word.split('/'):
            size = len(w)  # calculate on old size re-replace

            replace = replace_words.get(w.lower())
            if not replace and word.endswith("."):
                replace = replace_words.get(w[:-1].lower())

            if replace:
                w = replace

            new_offset_end = new_offset_start + size
            if w:
                data = (w, new_offset_start, new_offset_end)
                if w.isdigit():
                    unsplit_digits_and_spans.append(data)
                else:
                    handle_unsplit(unsplit_digits_and_spans, split_words_and_spans)
                    split_words_and_spans.append(data)
            new_offset_start = new_offset_end
        handle_unsplit(unsplit_digits_and_spans, split_words_and_spans)

    return split_words_and_spans


def sentences_and_offsets(txt):
    """ Split by nltk.sent_tokenize and new line """
    tokens = nltk.sent_tokenize(txt)

    split_sentences_and_offsets = []
    for sentence, sentence_offset, _ in tokens_and_spans(txt, tokens):
        new_offset = sentence_offset
        for line in sentence.split('\n'):
            if line:
                split_sentences_and_offsets.append((line, new_offset))
            new_offset += len(line) + 1

    return split_sentences_and_offsets


def process_text_phenotype(text_phenotype, phenotype_matcher):
    tokenized_text_and_spans = word_tokenise_and_spans(text_phenotype.text)

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

    words_and_spans = transform_words(words, tags, spans)

    try:
        parse_words(text_phenotype, words_and_spans, phenotype_matcher)
    except SkipAllPhenotypeMatchException:
        logging.info("Completely skipping: %s", text_phenotype.text)

    text_phenotype.processed = True
    text_phenotype.save()


def replace_comments_with_spaces(text):
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


def create_phenotype_description(text, phenotype_matcher=None):
    if phenotype_matcher is None:
        phenotype_matcher = PhenotypeMatcher()

    phenotype_description = PhenotypeDescription.objects.create(original_text=text)
    known_sentences = []
    unknown_sentences = []

    # Strip carriage returns - some browsers send \r\n in forms, Ajax usually sends \n
    # These need to be the same so save/ajax preview look the same
    text = text.replace('\r', '')
    text = replace_comments_with_spaces(text)  # Need to keep length the same for offsets
    for sentence, sentence_offset in sentences_and_offsets(text):
        text_phenotype, tp_created = TextPhenotype.objects.get_or_create(text=sentence)
        tps = TextPhenotypeSentence.objects.create(phenotype_description=phenotype_description,
                                                   text_phenotype=text_phenotype,
                                                   sentence_offset=sentence_offset)
        if tp_created:
            has_alpha_numeric = any([x.isalnum() for x in sentence])
            if not has_alpha_numeric:  # Impossible to match anything so skip
                text_phenotype.processed = True
                text_phenotype.save()
                known_sentences.append(tps)
            else:
                unknown_sentences.append(tps)
        else:
            known_sentences.append(tps)

    if unknown_sentences:
        for sentence in unknown_sentences:
            try:
                process_text_phenotype(sentence.text_phenotype, phenotype_matcher)
            except Exception:
                logging.error("Problem processing phenotype for '%s'", sentence.text_phenotype.text)
                raise
            known_sentences.append(sentence)

    return phenotype_description


def bulk_patient_phenotype_matching(patients=None):
    if patients is None:
        patients = Patient.objects.filter(phenotype__isnull=False)

    start = time.time()
    phenotype_matcher = PhenotypeMatcher()
    get_and_log_time_since(start, "load references")

    start = time.time()
    patients = list(patients)
    num_patients = len(patients)
    num_parsed_phenotypes = 0
    if num_patients:
        for i, p in enumerate(patients):

            parsed_phenotypes = p.process_phenotype_if_changed(phenotype_matcher=phenotype_matcher)
            num_parsed_phenotypes += parsed_phenotypes
            if not i % 50:
                perc_complete = 100.0 * i / num_patients
                logging.info("%d patients %0.2f%% complete", i, perc_complete)

        ts = get_and_log_time_since(start, "bulk_patient_phenotype_matching")
        if num_parsed_phenotypes:
            time_per_patient = ts / num_parsed_phenotypes
            pps = 1.0 / time_per_patient
            logging.info("%d parsed patient phenotypes - %.2f parsed per second/%.2f seconds per patient", num_parsed_phenotypes, pps, time_per_patient)

            total_matches = TextPhenotypeMatch.objects.count()
            tpm_qs = TextPhenotypeMatch.objects.values("match_type")
            tpm_qs = tpm_qs.annotate(count=Count("pk")).values_list("match_type", "count")
            logging.info("%d match types:", total_matches)
            for match_type, count in tpm_qs:
                logging.info("%s: %d", match_type, count)
    else:
        logging.info("No patients")
