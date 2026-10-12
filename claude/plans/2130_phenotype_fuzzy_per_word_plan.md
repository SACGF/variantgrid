# #2130 Phenotype fuzzy matching: one typo in one word, never a short word or a prefix

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-09
Status: in progress

Builds on `claude/plans/2131_phenotype_matcher_version_plan.md`: the change here alters stored matches, so it bumps
`PHENOTYPE_MATCHER_VERSION` to 2 and registers the rematch.

`PhenotypeMatcher.get_id_from_fuzzy_match` (`patients/phenotype_matcher.py`) accepted any term within Levenshtein
distance 1 of the whole joined text and returns the first one in dictionary order. One edit is enough to flip the
meaning of a clinical phrase when it lands in an abbreviation or a prefix:

| Input | Fuzzy match today | What went wrong |
|---|---|---|
| prolonged qt | Prolonged prothrombin time (alias "Prolonged PT") | QT → PT |
| recurrent urtis | Recurrent urinary tract infections (alias "Recurrent UTIs") | URTI → UTI |
| afebrile seizures | Febrile seizure | negation prefix dropped |
| decreased in | Decreased INR | 2-letter word edited (today a `COMMON_WORDS` entry) |

Only the special-case overrides keep these right; any phrase without one is exposed. A typo, by contrast, lives inside
a word that is long enough to have one: maplem → maple, xenterocyte → enterocyte, macrocephaky → macrocephaly.

## Rule

A fuzzy match is a term whose words line up one-to-one with the input, all identical except exactly one, where that
pair is a typo. `PhenotypeMatcher._is_typo_of(word, term_word)` is true when:

- both words are at least `MIN_LENGTH_SINGLE_WORD_FUZZY_MATCH` (5) characters, so abbreviations never fuzz
  (qt/pt, urtis/utis, in/inr);
- the input word is not a dictionary word (`_is_dictionary_word`, which already holds the NLTK corpus plus every
  ontology word) - a correctly spelled word isn't a typo; this is the single-word rule from #2125 applied per word,
  and it covers afebrile, asymptomatic and atypical;
- neither word is `"a" + the other` - the one negation prefix reachable in one edit (areflexia/reflexia is not in the
  dictionary);
- `Levenshtein.distance(word, term_word) == 1`.

That replaces `get_id_from_fuzzy_match`, `get_id_from_single_word_fuzzy_match`, `get_id_from_multi_word_fuzzy_match`
and `calculate_match_distance` (the "ae" allowance it gave is covered exactly by `_spelling_key` since #2125):

- single word: candidates are `single_words_by_length` at length ±1 (as now); the first candidate passing
  `_is_typo_of` wins;
- several words: candidates are the terms sharing a word with the input (`word_lookup`, as now, with the ±1 total
  length prefilter kept as the cheap cut); keep those with the same word count whose words differ at exactly one
  position, and the first whose differing pair passes `_is_typo_of` wins.

Candidates at distance 0 never reach here (exact matching, with spelling and plural variants, runs first), so every
accepted candidate is equally close and "first wins" stays deterministic.

`'decreased in'` leaves `COMMON_WORDS`: the rule now rejects it. The special-case overrides for the phrases above stay
(they still name the intended term); only the fallback changes.

## Version

`PHENOTYPE_MATCHER_VERSION = 2` with a changelog line `# 2: fuzzy matching is one typo in one word (#2130)`, and a new
`annotation/` migration registering `match_patient_phenotypes --stale` with `requires=["ontology-imported"]` and a
test of "any stale `TextPhenotype` exists" expressed on the historical model (`processed=True` and not
`match_version__matcher_version=2`).

## Tests (`patients/tests/test_phenotype_matching.py`)

Test terms added alongside the existing ones (the real HPO aliases): HP:0008151 Prolonged prothrombin time
["Prolonged PT"], HP:0000010 Recurrent urinary tract infections ["Recurrent UTIs"], HP:0002373 Febrile seizure
["Febrile seizures"], HP:0012514 Lower limb pain ["Leg pain"].

A phrase → expected-ids table run through `_match_ids` on a matcher whose three special-case dicts are emptied, so the
fallback is what is tested:

| phrase | expected |
|---|---|
| prolonged qt | nothing |
| recurrent urtis | nothing |
| afebrile seizures | nothing |
| decreased in | nothing |
| leg pains | HP:0012514 (plural, exact path) |
| febrile siezures | HP:0002373 (typo in a long word) |
| recurrent urinary tract infectons | HP:0000010 (typo in one word of a phrase) |

The existing `test_mispellings` (maplem, xenterocyte), `test_dictionary_word_not_fuzzy_matched` and
`test_exact_match_stops_fuzzy_match_in_other_ontologies` stay as the kept-behaviour guards; `test_skip_words` keeps
"decreased in" → nothing with its docstring updated to say why.
