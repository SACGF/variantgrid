from genes.gene_matching import tokenize_gene_symbols, GeneSymbolMatcher
from genes.models import GeneList
from snpdb.models import ImportStatus
from upload.models import UploadedGeneList
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


def create_gene_list(user, category, name, gene_names_set, modification_info=None, gene_matcher=None):
    if gene_matcher is None:
        gene_matcher = GeneSymbolMatcher()

    gene_list = GeneList(category=category, name=name, user=user, import_status=ImportStatus.IMPORTING)
    gene_list.save()

    if gene_names_set:
        gene_matcher.create_gene_list_gene_symbols(gene_list, gene_names_set, modification_info)

    gene_list.import_status = ImportStatus.SUCCESS
    gene_list.save()
    return gene_list


class ImportGeneListTask(ImportTask):
    MIN_GENES_TO_USE_CACHING_GENE_MATCHER = 10

    def process_items(self, uploaded_file):

        uploaded_gene_list, _ = UploadedGeneList.objects.get_or_create(uploaded_file=uploaded_file)
        with open(uploaded_file.get_filename()) as f:
            gene_list_data = f.read()
        gene_names_set = tokenize_gene_symbols(gene_list_data)

        if len(gene_names_set) > self.MIN_GENES_TO_USE_CACHING_GENE_MATCHER:
            gene_matcher = GeneSymbolMatcher()
        else:
            gene_matcher = None

        modification_info = "From uploaded gene_list: %s" % uploaded_file.get_filename()
        gene_list = create_gene_list(uploaded_file.user, None, uploaded_file.name, gene_names_set, modification_info,
                                     gene_matcher=gene_matcher)
        uploaded_gene_list.gene_list = gene_list
        uploaded_gene_list.save()

        return gene_list.genelistgenesymbol_set.count()


ImportGeneListTask = app.register_task(ImportGeneListTask())  # @UndefinedVariable
