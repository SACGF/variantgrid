"""
Adversarial unit tests for library utilities.
Pure-function tests — no DB needed (TestCase used for consistency).
"""
import os
import tempfile
import time
from datetime import UTC, date, datetime, timedelta
from decimal import Decimal

from django.test import TestCase

from library.cache import timed_cache
from library.log_utils import NotificationBuilder
from library.unit_percent import (
    convert_from_percent_to_unit,
    convert_from_unit_to_percent,
    server_side_format_percent,
)
from library.utils.color_utils import rgb_contrasting_text, rgb_relative_luminance
from library.utils.collection_utils import (
    IterableStitcher,
    LimitedCollection,
    batch_iterator,
    flatten_nested_lists,
    group_by_key,
)
from library.utils.date_utils import calculate_age, parse_yymm, utc_from_timestamp
from library.utils.file_utils import IteratorFile, file_to_array
from library.utils.hash_utils import string_deterministic_hash
from library.utils.json_utils import JsonDiffs, make_json_safe_in_place, strip_json
from library.utils.text_utils import (
    format_significant_digits,
    join_with_commas_and_ampersand,
    limit_str,
    pretty_label,
    split_dict_multi_values,
)

# ---------------------------------------------------------------------------
# unit_percent.py
# ---------------------------------------------------------------------------

class TestServerSideFormatPercent(TestCase):
    def test_zero_returns_empty_string(self):
        # Docstring says "Shows falsey values (eg 0.0) as blank" — 0.0 AF → blank
        self.assertEqual(server_side_format_percent(0.0), "")

    def test_nonzero_formats_with_percent(self):
        self.assertEqual(server_side_format_percent(0.5), "0.5%")

    def test_rare_frequency_is_not_scientific_notation(self):
        self.assertEqual(server_side_format_percent(0.0000123), "0.0000123%")

    def test_none_as_missing_value(self):
        self.assertEqual(server_side_format_percent(None, missing_value=None), "")

    def test_custom_missing_value_returns_empty(self):
        self.assertEqual(server_side_format_percent(-1, missing_value=-1), "")


class TestConvertPercentUnit(TestCase):
    def test_roundtrip(self):
        for pct in [1.0, 50.0, 100.0]:
            unit = convert_from_percent_to_unit(pct)
            self.assertAlmostEqual(convert_from_unit_to_percent(unit), pct)

    def test_zero_roundtrip(self):
        self.assertAlmostEqual(convert_from_percent_to_unit(0.0), 0.0)
        self.assertAlmostEqual(convert_from_unit_to_percent(0.0), 0.0)

    def test_missing_value_passthrough_none(self):
        self.assertIsNone(convert_from_percent_to_unit(None))
        self.assertIsNone(convert_from_unit_to_percent(None))


# ---------------------------------------------------------------------------
# text_utils.py
# ---------------------------------------------------------------------------

class TestFormatSignificantDigits(TestCase):
    def test_zero(self):
        self.assertEqual(format_significant_digits(0), "0")

    def test_integers_kept_whole(self):
        self.assertEqual(format_significant_digits(1), "1")
        self.assertEqual(format_significant_digits(10000), "10000")
        # 123456 → 3 sig figs → 123000
        self.assertEqual(format_significant_digits(123456), "123000")
        self.assertEqual(format_significant_digits(456.12), "456")

    def test_float_trailing_zeros_stripped(self):
        # Trailing zeros after decimal must be stripped, not left as "1.200000000000"
        self.assertEqual(format_significant_digits(1.2), "1.2")
        self.assertEqual(format_significant_digits(1.10004), "1.1")

    def test_float_exact_three_sig_figs(self):
        self.assertEqual(format_significant_digits(1.23456), "1.23")

    def test_rounds_half_up(self):
        self.assertEqual(format_significant_digits(1.114), "1.11")
        self.assertEqual(format_significant_digits(1.116), "1.12")

    def test_small_float(self):
        self.assertEqual(format_significant_digits(0.0001234), "0.000123")
        self.assertEqual(format_significant_digits(0.00000150002), "0.0000015")

    def test_negative_number(self):
        self.assertEqual(format_significant_digits(-1.234567), "-1.23")
        self.assertEqual(format_significant_digits(-0.00000150002), "-0.0000015")

    def test_custom_sig_digits(self):
        self.assertEqual(format_significant_digits(3.14159, sig_digits=2), "3.1")


class TestSlackMarkdownToHtml(TestCase):
    def test_emoji_stripped_and_bold_converted(self):
        self.assertEqual(NotificationBuilder.slack_markdown_to_html(":hospital: This is **bold**"),
                         "This is <strong>bold</strong>")


class TestUtcFromTimestamp(TestCase):
    def test_epoch_is_utc_aware(self):
        self.assertEqual(utc_from_timestamp(0), datetime(1970, 1, 1, tzinfo=UTC))


class TestSplitDictMultiValues(TestCase):
    def test_basic(self):
        data = {"a": "1|2|3", "b": "x|y|z"}
        result = split_dict_multi_values(data, "|")
        self.assertEqual(result, [{"a": "1", "b": "x"}, {"a": "2", "b": "y"}, {"a": "3", "b": "z"}])

    def test_unequal_splits_raise(self):
        # "a" splits into 3, "b" into 2 - mismatched counts must not silently drop keys
        data = {"a": "1|2|3", "b": "x|y"}
        with self.assertRaises(ValueError):
            split_dict_multi_values(data, "|")


class TestJoinWithCommasAndAmpersand(TestCase):
    def test_empty(self):
        self.assertEqual(join_with_commas_and_ampersand([]), "")

    def test_single(self):
        self.assertEqual(join_with_commas_and_ampersand(["Alice"]), "Alice")

    def test_two(self):
        self.assertEqual(join_with_commas_and_ampersand(["Alice", "Bob"]), "Alice & Bob")

    def test_three(self):
        self.assertEqual(join_with_commas_and_ampersand(["Alice", "Bob", "Carol"]), "Alice, Bob & Carol")


class TestLimitStr(TestCase):
    def test_exactly_at_limit_no_truncation(self):
        # Uses `>` not `>=` so exactly limit chars is NOT truncated
        self.assertEqual(limit_str("hello", 5), "hello")

    def test_one_over_limit(self):
        self.assertEqual(limit_str("hello!", 5), "hello...")

    def test_empty_string(self):
        self.assertEqual(limit_str("", 5), "")


class TestPrettyLabel(TestCase):
    def test_underscore_to_space_and_capitalize(self):
        self.assertEqual(pretty_label("hello_world"), "Hello World")

    def test_empty(self):
        self.assertEqual(pretty_label(""), "")


# ---------------------------------------------------------------------------
# date_utils.py
# ---------------------------------------------------------------------------

class TestCalculateAge(TestCase):
    def test_birthday_today(self):
        today = date.today()
        born = date(today.year - 30, today.month, today.day)
        self.assertEqual(calculate_age(born), 30)

    def test_birthday_is_tomorrow(self):
        today = date.today()
        tomorrow = today + timedelta(days=1)
        born = date(today.year - 30, tomorrow.month, tomorrow.day)
        self.assertEqual(calculate_age(born), 29)  # not yet had birthday

    def test_none_born_returns_none(self):
        self.assertIsNone(calculate_age(None))

    def test_with_died_day_before_birthday(self):
        born = date(1990, 6, 15)
        died = date(2020, 6, 14)
        self.assertEqual(calculate_age(born, died), 29)


class TestGetMonthAndYear(TestCase):
    def test_yymm_string(self):
        # "2201" → month=01, year=22
        month, year = parse_yymm("2201")
        self.assertEqual(month, 1)
        self.assertEqual(year, 22)

    def test_invalid_yyyymm_format_raises(self):
        with self.assertRaises(ValueError):
            parse_yymm("202201")


# ---------------------------------------------------------------------------
# file_utils.py
# ---------------------------------------------------------------------------

class TestFileToArray(TestCase):
    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
        self.tmp.write("line1\nline2\nline3\nline4\nline5\n")
        self.tmp.close()

    def tearDown(self):
        os.unlink(self.tmp.name)

    def test_no_max_lines_reads_all(self):
        self.assertEqual(len(file_to_array(self.tmp.name)), 5)

    def test_max_lines_is_a_count_not_an_index(self):
        result = file_to_array(self.tmp.name, max_lines=3)
        self.assertEqual(len(result), 3, f"max_lines=3 returned {len(result)} lines")

    def test_max_lines_zero(self):
        result = file_to_array(self.tmp.name, max_lines=0)
        self.assertEqual(len(result), 0, f"max_lines=0 should return 0 lines, got {len(result)}")

    def test_comment_skip(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("# comment\nreal_line\n# another comment\nmore_data\n")
            fname = f.name
        try:
            self.assertEqual(file_to_array(fname, comment="#"), ["real_line", "more_data"])
        finally:
            os.unlink(fname)

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            fname = f.name
        try:
            self.assertEqual(file_to_array(fname), [])
        finally:
            os.unlink(fname)


class TestIteratorFile(TestCase):
    def _make(self, chunks):
        return IteratorFile(iter(chunks))

    def test_read_across_chunks(self):
        f = self._make(["hel", "lo"])
        self.assertEqual(f.read(5), "hello")

    def test_read_exhausted_returns_empty_string(self):
        f = self._make(["hi"])
        f.read(100)
        self.assertEqual(f.read(1), '')

    def test_readline_across_chunks(self):
        f = self._make(["lin", "e1\n", "line2\n"])
        self.assertEqual(f.readline(), "line1\n")

    def test_readline_no_newline_at_end(self):
        f = self._make(["hello"])
        self.assertEqual(f.readline(), "hello")


# ---------------------------------------------------------------------------
# collection_utils.py
# ---------------------------------------------------------------------------

class TestBatchIterator(TestCase):
    def test_exact_multiple(self):
        self.assertEqual(list(batch_iterator([1, 2, 3, 4], batch_size=2)), [[1, 2], [3, 4]])

    def test_with_remainder(self):
        self.assertEqual(list(batch_iterator([1, 2, 3, 4, 5], batch_size=2)), [[1, 2], [3, 4], [5]])

    def test_empty(self):
        self.assertEqual(list(batch_iterator([], batch_size=2)), [])

    def test_batch_larger_than_input(self):
        self.assertEqual(list(batch_iterator([1, 2, 3], batch_size=10)), [[1, 2, 3]])


class TestGroupByKey(TestCase):
    def test_basic_sorted(self):
        data = [("a", 1), ("a", 2), ("b", 3)]
        result = list(group_by_key(data, key=lambda x: x[0]))
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0][0], "a")
        self.assertEqual(len(result[0][1]), 2)

    def test_none_key_raises(self):
        data = [("a", 1), (None, 2)]
        with self.assertRaises(ValueError):
            list(group_by_key(data, key=lambda x: x[0]))

    def test_empty_input(self):
        self.assertEqual(list(group_by_key([], key=lambda x: x)), [])


class TestIterableStitcher(TestCase):
    def test_merges_two_sorted_iterables(self):
        self.assertEqual(list(IterableStitcher([[1, 3, 5], [2, 4, 6]])), [1, 2, 3, 4, 5, 6])

    def test_empty_iterables(self):
        self.assertEqual(list(IterableStitcher([[], []])), [])

    def test_one_empty_one_not(self):
        self.assertEqual(list(IterableStitcher([[], [1, 2, 3]])), [1, 2, 3])

class TestLimitedCollection(TestCase):
    def test_over_limit(self):
        lc = LimitedCollection([1, 2, 3, 4, 5], limit=3)
        self.assertTrue(lc.is_limited)
        self.assertEqual(len(lc), 5)           # true_count
        self.assertEqual(list(lc), [1, 2, 3])  # only limit items iterated

    def test_bool_empty(self):
        self.assertFalse(bool(LimitedCollection([], limit=5)))

    def test_none_data(self):
        lc = LimitedCollection(None, limit=5)
        self.assertEqual(list(lc), [])
        self.assertFalse(bool(lc))

    def test_exactly_at_limit_not_limited(self):
        self.assertFalse(LimitedCollection([1, 2, 3], limit=3).is_limited)


class TestFlattenNestedLists(TestCase):
    def test_single_level_nested(self):
        self.assertEqual(flatten_nested_lists([[1, 2], [3, 4]]), [1, 2, 3, 4])

    def test_deeply_nested(self):
        self.assertEqual(flatten_nested_lists([[1, [2, [3]]]]), [1, 2, 3])

    def test_none_filtered_out(self):
        # None values are silently dropped during flattening
        self.assertEqual(flatten_nested_lists([1, None, 2]), [1, 2])

    def test_empty(self):
        self.assertEqual(flatten_nested_lists([]), [])


# ---------------------------------------------------------------------------
# json_utils.py
# ---------------------------------------------------------------------------

class TestStripJson(TestCase):
    def test_removes_none(self):
        self.assertEqual(strip_json({"a": None}), {})

    def test_removes_empty_string(self):
        self.assertEqual(strip_json({"a": ""}), {})

    def test_removes_false(self):
        self.assertEqual(strip_json({"a": False}), {})

    def test_removes_empty_list(self):
        self.assertEqual(strip_json({"a": []}), {})

    def test_keeps_zero_int(self):
        # 0 must NOT be stripped — it is a valid value (distinct from "absent")
        self.assertEqual(strip_json({"a": 0}), {"a": 0})

    def test_list_items_not_individually_stripped(self):
        # A non-empty list keeps the key even if all list items are falsey
        result = strip_json({"a": [None, False, ""]})
        self.assertIn("a", result)

    def test_nested_dict_stripped(self):
        result = strip_json({"outer": {"inner": None, "keep": 1}})
        self.assertEqual(result, {"outer": {"keep": 1}})

    def test_empty_dict_input(self):
        self.assertEqual(strip_json({}), {})

    def test_non_dict_passthrough(self):
        self.assertEqual(strip_json(42), 42)
        self.assertIsNone(strip_json(None))


class TestMakeJsonSafeInPlace(TestCase):
    def test_decimal_in_dict(self):
        d = {"val": Decimal("1.5")}
        make_json_safe_in_place(d)
        self.assertIsInstance(d["val"], float)
        self.assertAlmostEqual(d["val"], 1.5)

    def test_decimal_in_nested_dict(self):
        d = {"outer": {"inner": Decimal("2.0")}}
        make_json_safe_in_place(d)
        self.assertIsInstance(d["outer"]["inner"], float)

    def test_decimal_in_list(self):
        lst = [Decimal("3.0"), "hello"]
        make_json_safe_in_place(lst)
        self.assertIsInstance(lst[0], float)

    def test_no_modification_for_plain_types(self):
        d = {"a": 1, "b": "hello", "c": True}
        make_json_safe_in_place(d)
        self.assertEqual(d, {"a": 1, "b": "hello", "c": True})


class TestJsonDiffsDifferences(TestCase):
    def test_identical_dicts(self):
        self.assertFalse(bool(JsonDiffs.differences({"a": 1}, {"a": 1})))

    def test_simple_value_change(self):
        diffs = JsonDiffs.differences({"a": 1}, {"a": 2})
        self.assertEqual(len(diffs.json_diffs), 1)
        self.assertEqual(diffs.json_diffs[0].a, 1)
        self.assertEqual(diffs.json_diffs[0].b, 2)

    def test_missing_key(self):
        self.assertEqual(len(JsonDiffs.differences({"a": 1}, {}).json_diffs), 1)

    def test_nested_change(self):
        diffs = JsonDiffs.differences({"a": {"b": 1}}, {"a": {"b": 2}})
        self.assertEqual(len(diffs.json_diffs), 1)

    def test_list_with_id_uses_id_matching(self):
        a = [{"id": 1, "val": "x"}, {"id": 2, "val": "y"}]
        b = [{"id": 2, "val": "y"}, {"id": 1, "val": "changed"}]
        diffs = JsonDiffs.differences(a, b)
        # Should detect 1 change (id=1 val changed), not 2 spurious index diffs
        self.assertEqual(len(diffs.json_diffs), 1)

# ---------------------------------------------------------------------------
# hash_utils.py
# ---------------------------------------------------------------------------

class TestStringDeterministicHash(TestCase):
    def test_same_string_same_hash(self):
        self.assertEqual(string_deterministic_hash("hello"), string_deterministic_hash("hello"))


# ---------------------------------------------------------------------------
# cache.py
# ---------------------------------------------------------------------------

class TestTimedCache(TestCase):
    def test_basic_caching(self):
        call_count = 0

        @timed_cache()
        def fn(x):
            nonlocal call_count
            call_count += 1
            return x * 2

        fn(5)
        fn(5)
        self.assertEqual(call_count, 1)

    def test_ttl_expiry_recalculates(self):
        call_count = 0

        @timed_cache(ttl=0.05)
        def fn():
            nonlocal call_count
            call_count += 1
            return call_count

        first = fn()
        time.sleep(0.12)
        second = fn()
        self.assertNotEqual(first, second, "After TTL expiry, function should be re-called")

    def test_size_limit_evicts_oldest(self):
        call_tracker = []

        @timed_cache(size_limit=2, quick_key_access=True)
        def fn(x):
            call_tracker.append(x)
            return x

        fn(1)
        fn(2)
        fn(3)  # evicts key=1
        fn(1)  # 1 was evicted, must recalculate
        self.assertEqual(call_tracker.count(1), 2)


class TestRgbContrastingText(TestCase):
    """ Tag colors pick black/white text via luminance - see rgb_contrasting_text """

    @staticmethod
    def _contrast_ratio(rgb: str, other: str) -> float:
        luminances = sorted([rgb_relative_luminance(rgb), rgb_relative_luminance(other)])
        return (luminances[1] + 0.05) / (luminances[0] + 0.05)

    def test_readable_on_mid_tones(self):
        """ Mid-tone tag colors are where inverting the background used to fail """
        for rgb in ["#808080", "#4488cc", "#cc8844", "#996699", "#888844", "#556b7f"]:
            with self.subTest(rgb=rgb):
                ratio = self._contrast_ratio(rgb, rgb_contrasting_text(rgb))
                self.assertGreaterEqual(ratio, 4.5, f"{rgb} text fails WCAG AA")

    def test_picks_white_on_dark_black_on_light(self):
        self.assertEqual(rgb_contrasting_text("#000080"), "#ffffff")
        self.assertEqual(rgb_contrasting_text("#ffff00"), "#000000")
