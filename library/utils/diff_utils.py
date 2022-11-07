import difflib
import re
from dataclasses import dataclass
from html import escape
from typing import List

from django.utils.safestring import SafeString


@dataclass(frozen=True)
class DiffTextSegment:
    operation: str
    text: str

    @property
    def effective_text(self):
        return SafeString(escape(self.text).replace("\n", "<span style='font-size:x-small;opacity:0.5'>&#8726;n</span><br/>"))

    @property
    def operation_name(self):
        if self.operation == ' ':
            return 'same'
        elif self.operation == '-':
            return 'subtract'
        elif self.operation == '+':
            return 'add'
        else:
            return self.operation

    def __str__(self):
        return f"{self.operation} {self.text}"


class DiffBuilder:
    # TODO, rather than trying to do a bunch of smarts in append() maybe all the smarts are best saved for optimize

    def __init__(self):
        self.diff_segments: List[DiffTextSegment] = []
        self.same_text = ''
        self.add_text = ''
        self.sub_text = ''

    def apply(self):
        if sub_text := self.sub_text:
            self.diff_segments.append(DiffTextSegment(operation='-', text=sub_text))
        if add_text := self.add_text:
            self.diff_segments.append(DiffTextSegment(operation='+', text=add_text))
        if same_text := self.same_text:
            self.diff_segments.append(DiffTextSegment(operation=' ', text=same_text))
        self.sub_text = ''
        self.add_text = ''
        self.same_text = ''

    def optimize_add(self, subtract: str, add: str) -> List[DiffTextSegment]:
        prefix_same: str = ''
        suffix_same: str = ''
        while subtract and add and subtract[0] == add[0]:
            prefix_same += add[0]
            subtract = subtract[1:]
            add = add[1:]
        # while subtract and add and subtract[-1] == add[-1]:
        #     suffix_same += subtract[-1]
        #     subtract = subtract[:-1]
        #     add = add[:-1]
        items: List[DiffTextSegment] = []
        if prefix_same:
            items.append(DiffTextSegment(operation=' ', text=prefix_same))
        if subtract:
            items.append(DiffTextSegment(operation='-', text=subtract))
        if add:
            items.append(DiffTextSegment(operation='+', text=add))
        if suffix_same:
            items.append(DiffTextSegment(operation=' ', text=suffix_same))
        return items

    def optimize(self):
        self.apply()
        optimized = []
        check_starts_with = None
        for element in self.diff_segments:
            if element.operation == '-':
                if check_starts_with:
                    optimized.append(check_starts_with)
                check_starts_with = element
            else:
                if element.operation == '+' and check_starts_with:
                    if items := self.optimize_add(check_starts_with.text, element.text):
                        optimized.extend(items)
                else:
                    optimized.append(check_starts_with)
                    optimized.append(element)
                check_starts_with = None
        if check_starts_with:
            optimized.append(check_starts_with)
        self.diff_segments = optimized

    def _has_pending_add_subtract(self):
        return bool(self.add_text) or bool(self.sub_text)

    def _has_pending_same(self):
        return bool(self.same_text)

    def append(self, diff_text: str):
        operation = diff_text[0]
        text = diff_text[2:]
        if operation == '?':
            # '? 'line not present in either input sequence
            # Lines beginning with ‘?’ attempt to guide the eye to intraline differences,
            # and were not present in either input sequence. These lines can be confusing if the sequences contain tab characters.
            return
        if operation == ' ':
            if self._has_pending_add_subtract():
                if not text.strip():
                    # if we have "bear down" -> "fat dog" it would result in -bear +fat (same space) -down +dog
                    # much cleaner to have that as -"bear down" +"fat dog"
                    # might need to look into having whitespace_text that we track separately
                    self.sub_text += text
                    self.add_text += text
                    return

                self.apply()
            self.same_text += text
        else:
            if self._has_pending_same():
                self.apply()
            if operation == '-':
                self.sub_text += text
            else:
                self.add_text += text

    def __iter__(self):
        return iter(self.diff_segments)

    def __bool__(self):
        return bool(self.diff_segments)


def diff_text(a: str, b: str) -> DiffBuilder:

    def _tokenize(text: str) -> List[str]:
        return re.split(r'(\s)', text)

    diff_builder = DiffBuilder()
    for diff_chars in difflib.Differ().compare(_tokenize(a), _tokenize(b)):
        diff_builder.append(diff_chars)
    diff_builder.optimize()
    return diff_builder
