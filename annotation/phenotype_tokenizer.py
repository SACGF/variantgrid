import nltk

from library.utils import invert_dict_of_lists
from library.utils.nltk_utils import ensure_nltk_data


class PhenotypeTokenizer:
    def __init__(self):
        # Make this a class - so we can ensure that we have NLTK resources loaded
        # for any dependent methods
        ensure_nltk_data('taggers/averaged_perceptron_tagger_eng')  # Needed for word_tokenize
        ensure_nltk_data('tokenizers/punkt_tab')

    @staticmethod
    def _initial_replace_characters(text):
        REPLACE_CHARS = {';': ' '}

        cleaned_text = text
        for from_char, to_char in REPLACE_CHARS.items():
            cleaned_text = cleaned_text.replace(from_char, to_char)

        return cleaned_text

    @staticmethod
    def _tokens_and_spans(txt, tokens):
        offset = 0

        tokenized_text_and_spans = []
        for token in tokens:
            offset = txt.find(token, offset)
            tokenized_text_and_spans.append((token, offset, offset + len(token)))
            offset += len(token)

        return tokenized_text_and_spans

    @staticmethod
    def _handle_unsplit(unsplit_digits_and_spans, split_words_and_spans):
        if unsplit_digits_and_spans:
            unsplit_word = '/'.join([d[0] for d in unsplit_digits_and_spans])
            my_offset_start = unsplit_digits_and_spans[0][1]
            my_offset_end = unsplit_digits_and_spans[-1][2]
            split_words_and_spans.append((unsplit_word, my_offset_start, my_offset_end))
            unsplit_digits_and_spans[:] = []

    def sentences_and_offsets(self, txt):
        """ Split by nltk.sent_tokenize and new line """
        tokens = nltk.sent_tokenize(txt)

        split_sentences_and_offsets = []
        for sentence, sentence_offset, _ in self._tokens_and_spans(txt, tokens):
            new_offset = sentence_offset
            for line in sentence.split('\n'):
                if line:
                    split_sentences_and_offsets.append((line, new_offset))
                new_offset += len(line) + 1

        return split_sentences_and_offsets

    def word_tokenise_and_spans(self, txt):
        """ Split by nltk.sent_tokenize and / (only if alphanumeric not for dates) """

        cleaned_text = self._initial_replace_characters(txt)
        tokens = nltk.word_tokenize(cleaned_text)
        split_words_and_spans = []

        # {TO: FROM} - These have to be the same size! Then switch them out later
        PRE_SPLIT_REPLACE = {"dv/p": "dv_p"}

        # replace word : list of alternate words to replace
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

        for word, offset_start, _ in self._tokens_and_spans(cleaned_text, tokens):
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
                        self._handle_unsplit(unsplit_digits_and_spans, split_words_and_spans)
                        split_words_and_spans.append(data)
                new_offset_start = new_offset_end
            self._handle_unsplit(unsplit_digits_and_spans, split_words_and_spans)

        return split_words_and_spans
