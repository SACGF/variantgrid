from django.template import Library

from snpdb.forms import ProjectChoiceForm, VariantsTypeMultipleChoiceForm
from snpdb.models import VariantsType

register = Library()


@register.inclusion_tag("snpdb/tags/vcf_grid_filter_tags.html", takes_context=True)
def vcf_grid_filter(context, table_id, variants_type=False):
    # ProjectChoiceForm is shown twice (Samples and VCF) so we need to alter the IDs
    project_form = ProjectChoiceForm()
    project_form.fields["project"].widget.attrs["id"] = f"id_project_{table_id}"

    context = {
        'table_id': table_id,
        "project_form": project_form,
    }

    if variants_type:
        vt_initial = {"variants_type": [c[0] for c in VariantsType.choices]}
        context["variants_type_form"] = VariantsTypeMultipleChoiceForm(initial=vt_initial)
        context["variants_type_labels"] = dict(VariantsType.choices)

    return context
