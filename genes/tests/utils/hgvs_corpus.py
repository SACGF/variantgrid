import pathlib
from collections import defaultdict

HGVS_CORPUS_PATH = pathlib.Path(__file__).parent.parent / "test_data" / "hgvs_corpus.tsv"


def load_hgvs_corpus() -> dict[str, list[str]]:
    """ {category: [hgvs]} - see the file header for what each category is and where it came from """
    by_category = defaultdict(list)
    for line in HGVS_CORPUS_PATH.read_text().splitlines():
        if line.startswith("#"):
            continue
        category, _, hgvs = line.partition("\t")
        by_category[category].append(hgvs)
    return dict(by_category)


def all_hgvs() -> list[str]:
    return [hgvs for hgvs_list in load_hgvs_corpus().values() for hgvs in hgvs_list]
