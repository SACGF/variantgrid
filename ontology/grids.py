import abc
from collections import defaultdict

import pandas as pd

from library.pandas_jqgrid import DataFrameJqGrid
from ontology.models import OntologyTerm, OntologySnake


class AbstractOntologyGenesGrid(DataFrameJqGrid, abc.ABC):
    @abc.abstractmethod
    def _get_ontology_term_ids(self):
        pass

    def get_dataframe(self):
        # This uses the same method as gene filter (special_case_gene_symbols_for_hpo_and_omim) though with individual
        # calls per term so that it matches what gene filters is doing
        terms_dict = OntologyTerm.split_hpo_omim_mondo_as_dict(self._get_ontology_term_ids())
        gene_terms_set = defaultdict(lambda: defaultdict(set))

        for ontology_name, terms_qs in terms_dict.items():
            for ot in terms_qs:
                for gene in OntologySnake.cached_gene_symbols_for_terms_tuple((ot,)):
                    gene_terms_set[gene.symbol][ontology_name].add(str(ot))

        gene_dict = {k: {t: ", ".join(sorted(term_set)) for t, term_set in v.items()} for k, v in gene_terms_set.items()}
        df = pd.DataFrame.from_dict(gene_dict, orient='index')
        return df.sort_index()
