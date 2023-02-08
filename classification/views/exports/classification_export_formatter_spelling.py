import re
from collections import defaultdict
from functools import cached_property
from itertools import chain
from typing import Set, List, Dict

import nltk
from django.http import HttpRequest

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter2
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row

RE_HAS_BAD_CHAR = re.compile(r"[\d._]")


class ClassificationSpellingRow(ExportRow):

    def __init__(self, cm: ClassificationModification, spell: 'SpellChecker'):
        self.cm = cm
        self.spell = spell

    @export_column()
    def lab(self) -> str:
        return str(self.cm.classification.lab)

    @export_column("URL")
    def url(self) -> str:
        return get_url_from_view_path(self.cm.classification.get_absolute_url())

    @export_column("Interpretation Summary")
    def interpretation_summary(self) -> str:
        return self.cm.get(SpecialEKeys.INTERPRETATION_SUMMARY)

    @staticmethod
    def fix_word_token(word) -> List[str]:
        for bad_char in RE_HAS_BAD_CHAR.finditer(word):
            return []
        if "/" in word:
            words = word.split("/")
        else:
            words = [word]
        return [word for word in words if len(word) >= 4]

    @cached_property
    def suspect_words_set(self) -> Set[str]:
        if interpretation_summary := self.cm.get(SpecialEKeys.INTERPRETATION_SUMMARY):
            words = nltk.word_tokenize(interpretation_summary)
            words_2d = [ClassificationSpellingRow.fix_word_token(word) for word in words]
            words = list(chain.from_iterable(words_2d))

            return self.spell.unknown(words)
        return set()

    @export_column("Suspect Words")
    def spelling_warnings(self) -> str:
        if unknown := self.suspect_words_set:
            return "\n".join(sorted(unknown))


@register_classification_exporter("spelling")
class ClassificationExportFormatterSpelling(ClassificationExportFormatter2):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterSpelling':
        return ClassificationExportFormatterSpelling(
            classification_filter=ClassificationFilter.from_request(request)
        )

    def __init__(self, classification_filter: ClassificationFilter):
        from spellchecker import SpellChecker
        self.spell = SpellChecker()
        self.suspect_count: Dict[str, int] = defaultdict(int)
        super().__init__(classification_filter=classification_filter)

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
        return [delimited_row(ClassificationSpellingRow.csv_header(), ",")]

    def footer(self) -> List[str]:
        suspect_count_list = [(word, count) for word, count in self.suspect_count.items()]
        suspect_count_list.sort(key=lambda x: x[1], reverse=True)
        return [delimited_row([
            "", "", "Total Classifications Occurred In", "\n".join(f"{wc[0]}: {wc[1]}" for wc in suspect_count_list)
        ])]

    def row(self, allele_data: AlleleData) -> List[str]:
        rows: List[str] = []
        for cm in allele_data.cms:
            spelling_row = ClassificationSpellingRow(cm, self.spell)
            for suspect_word in spelling_row.suspect_words_set:
                self.suspect_count[suspect_word] += 1
            rows.append(
                delimited_row(spelling_row.to_csv())
            )
        return rows
