import difflib
import itertools
import re
from dataclasses import dataclass, field
from html import escape
from typing import List, Optional, Any, Pattern

from django.utils.safestring import SafeString

from library.utils import first


@dataclass
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
        elif self.operation == 'd':  # only used by multi-diff
            return 'different'
        elif self.operation == '?':  # only used by multi-diff
            return 'unknown'
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


@dataclass(frozen=True)
class MultiDiffInput:
    identifier: Any
    text: str

@dataclass(frozen=True)
class MultiDiffOutput:
    input: MultiDiffInput
    parts: Optional[List[str]]
    matched: Optional[bool]
    diffs: List[DiffTextSegment] = field(default_factory=list)

    @property
    def identifier(self):
        return self.input.identifier

    def __str__(self):
        return " ".join(str(diff) for diff in self.diffs)

    def html(self):
        if len(self.diffs) == 1 and not self.diffs[0].text:
            return ""

        return SafeString("<span class='diff-text'>" + "".join(
            f"<span class='diff-text-{diff.operation_name}'>{escape(diff.text)}</span>" for diff in self.diffs
        ) + "</span>")

    def append(self, op: str, segment: str) -> 'MultiDiffOutput':
        if not self.diffs or self.diffs[-1].operation != op:
            self.diffs.append(DiffTextSegment(op, segment))
        else:
            self.diffs[-1].text += segment
        return self

    @staticmethod
    def from_input(input: MultiDiffInput, pattern: Pattern):
        text = input.text or ''
        parts: List[str]
        if match := pattern.match(text):
            parts = match.groups()
            return MultiDiffOutput(
                input=input,
                parts=parts,
                matched=True
            )
        else:
            return MultiDiffOutput(
                input=input,
                parts=None,
                matched=False,
            ).append('u', text)

class MultiDiff:

    @staticmethod
    def diff_index(parts: List[str]) -> Optional[int]:
        for index, by_the_letter in enumerate(itertools.zip_longest(*parts, fillvalue='$')):
            if any(x != first(by_the_letter) for x in by_the_letter):
                return index
        return None

    def __init__(self, re_parts: Pattern):
        self.re_parts = re_parts

    def diffs(self, compare: List[MultiDiffInput]) -> List[MultiDiffOutput]:
        length = len(compare)
        if length == 0:
            return []
        elif length == 1:
            first_compare = compare[0]
            return [MultiDiffOutput(input=first_compare, matched=None, parts=[]).append(' ', first_compare.text)]
        else:
            outputs = [MultiDiffOutput.from_input(input, self.re_parts) for input in compare]
            comparing = [output for output in outputs if output.matched]
            for index in range(0, self.re_parts.groups):
                compare_parts = [compare.parts[index] or '' for compare in comparing]
                diff_index = MultiDiff.diff_index(compare_parts)
                for diff_output, part in zip(comparing, compare_parts):
                    if diff_index is None:
                        diff_output.append(' ', part)
                        continue
                    if diff_index > 0:
                        diff_output.append(' ', part[:diff_index])
                    diff_output.append('d', part[diff_index:])

            return outputs
