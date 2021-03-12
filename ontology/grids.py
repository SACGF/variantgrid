import abc
from collections import defaultdict

import pandas as pd
from library.pandas_jqgrid import DataFrameJqGrid
from ontology.models import OntologyTerm, OntologySnake


class AbstractOntologyGenesGrid(DataFrameJqGrid, abc.ABC):
    @abc.abstractmethod
    def _get_ontology_terms_ids(self):
        pass

    def get_dataframe(self):
        # This uses the same method as gene filter (special_case_gene_symbols_for_hpo_and_omim) though with individual
        # calls per term so that it matches what gene filters is doing
        hpo_qs, omim_qs = OntologyTerm.split_hpo_and_omim(self._get_ontology_terms_ids())
        gene_terms_set = defaultdict(lambda: defaultdict(set))
        for hpo in hpo_qs:
            for gene in OntologySnake.gene_symbols_for_terms([hpo]):
                gene_terms_set[gene.symbol]["hpo"].add(str(hpo))

        for omim in omim_qs:
            for gene in OntologySnake.gene_symbols_for_terms([omim]):
                gene_terms_set[gene.symbol]["omim"].add(str(omim))

        gene_dict = {k: {t: ", ".join(sorted(term_set)) for t, term_set in v.items()} for k, v in gene_terms_set.items()}
        df = pd.DataFrame.from_dict(gene_dict, orient='index')
        return df.sort_index()
