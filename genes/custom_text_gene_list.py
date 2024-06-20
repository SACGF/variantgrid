from django.contrib.auth.models import User

from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneListCategory, GeneList, CustomTextGeneList
from library.utils import sha256_str
from snpdb.models import ImportStatus


def create_custom_text_gene_list(custom_text_gene_list: CustomTextGeneList,
                                 username, gene_list_category_name=None, hidden=False, gene_matcher=None):
    """ Creates an associated gene_list for passed CustomTextGeneList object.
        Can pass gene_matcher instance to save reloading each call """

    if gene_matcher is None:
        gene_matcher = GeneSymbolMatcher()

    sha256_hash = sha256_str(custom_text_gene_list.text)
    if custom_text_gene_list.sha256_hash == sha256_hash:
        return  # No meaningful change

    user = User.objects.get(username=username)
    if custom_text_gene_list.gene_list:
        custom_text_gene_list.gene_list.delete()

    gene_list_kwargs = {}
    if gene_list_category_name:
        gene_list_kwargs["category"] = GeneListCategory.get_or_create_category(gene_list_category_name, hidden)

    gene_list = GeneList.objects.create(name=custom_text_gene_list.name,
                                        user=user,
                                        **gene_list_kwargs)
    try:
        gene_matcher.create_gene_list_gene_symbols_from_text(gene_list, custom_text_gene_list.text)
        gene_list.import_status = ImportStatus.SUCCESS
    except Exception as e:
        gene_list.error_message = str(e)
        gene_list.import_status = ImportStatus.ERROR

    gene_list.save()
    custom_text_gene_list.sha256_hash = sha256_hash
    custom_text_gene_list.gene_list = gene_list
    custom_text_gene_list.save()
