from django.template import Library

from snpdb.forms import ProjectChoiceForm
from snpdb.models import GenomeBuild
from snpdb.templatetags.jqgrid_tags import get_grid_id

register = Library()


@register.inclusion_tag("snpdb/tags/vcf_grid_filter_tags.html", takes_context=True)
def vcf_grid_filter(context, grid_id_suffix):
    grid_id = get_grid_id(grid_id_suffix)
    # ProjectChoiceForm is shown twice (Samples and VCF) so we need to alter the IDs
    project_form = ProjectChoiceForm()
    project_form.fields["project"].widget.attrs["id"] = f"id_project_{grid_id}"

    context = {
        'grid_id': grid_id,
        "project_form": project_form,
        "genome_builds": GenomeBuild.builds_with_annotation(),
    }
    return context
