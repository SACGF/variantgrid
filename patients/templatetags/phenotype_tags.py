import uuid

from django.template import Library

from genes.forms import GeneSymbolForm
from ontology.forms import HPOForm, OMIMForm

register = Library()


@register.inclusion_tag("patients/tags/phenotype_entry_tag.html", takes_context=True)
def phenotype_entry(context, form_entry_field, phenotype_description,
                    edit=True, show_page_help=True,
                    show_hpo=True, show_omim=True, show_genes=True, show_grid=False):
    """ form_entry_field : form field to enter phenotypes """

    if phenotype_description:
        patient_results = phenotype_description.get_results()
    else:
        patient_results = None

    initial_phenotype_text = ''
    if phenotype_description:
        initial_phenotype_text = phenotype_description.original_text

    flattened_uuid = str(uuid.uuid4()).replace("-", "_")

    hpo_form = None
    omim_form = None
    gene_symbol_form = None
    if show_hpo:
        hpo_form = HPOForm()
    if show_omim:
        omim_form = OMIMForm()
    if show_genes:
        gene_symbol_form = GeneSymbolForm()

    return {'flattened_uuid': flattened_uuid,
            'user': context['user'],
            'initial_phenotype_text': initial_phenotype_text,
            'edit': edit,
            'hpo_form': hpo_form,
            'omim_form': omim_form,
            'gene_symbol_form': gene_symbol_form,
            "form_entry_field": form_entry_field,
            "patient_results": patient_results,
            "show_grid": show_grid,
            "has_patient_results": bool(patient_results),
            "show_page_help": show_page_help}
