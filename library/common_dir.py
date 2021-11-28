import os
from collections import Counter

_end = '_end_'

def make_trie(*words):
    root = dict()
    for word in words:
        current_dict = root
        for letter in word:
            current_dict = current_dict.setdefault(letter, {})
        current_dict = current_dict.setdefault(_end, _end)
    return root

def get_multi_entries(trie, prefix=''):
    counter = Counter()
    trie_keys = list(trie.keys())

    for k in trie_keys:
        if k == _end:
            continue
        counter.update(get_multi_entries(trie[k], prefix + k))

    if len(counter) < 1:
        nk = len(trie_keys)
        if nk > 1:
            dirname = os.path.dirname(prefix)
            counter[dirname] += nk

    return counter


def get_common_prefix_dirs(files):
    trie = make_trie(*tuple(files))
    most_common = get_multi_entries(trie).most_common(3)
    return [k for k, _ in most_common]
