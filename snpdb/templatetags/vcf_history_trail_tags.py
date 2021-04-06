from django.template import Library

from seqauto.models import VCFFromSequencingRun

register = Library()


@register.inclusion_tag("snpdb/tags/vcf_history_trail.html", takes_context=True)
def vcf_history_trail(context, vcf):
    try:
        sequencing_run = vcf.vcffromsequencingrun.sequencing_run
    except VCFFromSequencingRun.DoesNotExist:
        sequencing_run = None

    user = context["user"]
    context.update({'sequencing_run': sequencing_run,
                    "vcf": vcf,
                    'can_view_vcf': user.has_perm('view_vcf', vcf)})
    return context
