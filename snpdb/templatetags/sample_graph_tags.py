from django.db.models.aggregates import Min, Max, Count
from django.template.library import Library

from library.date_utils import diff_month, date_to_month_year_string, month_year_string, month_range
from snpdb.models.models_vcf import VCF
import numpy as np

register = Library()


@register.inclusion_tag("snpdb/tags/samples_by_month_tag.html", takes_context=True)
def samples_by_month_graph(context):
    user = context["user"]

    vcf_qs = VCF.filter_for_user(user)
    data = vcf_qs.aggregate(earliest_date=Min("date"), latest_date=Max("date"))

    ed = data["earliest_date"]
    ld = data["latest_date"]

    tag_context = {}
    if ed and ld:
        num_months = diff_month(ld, ed)

        labels = month_range(ed.month, ed.year, 0, num_months, fmt=month_year_string)
        samples_by_month = {m: 0 for m in labels}

        for d, count in vcf_qs.annotate(count=Count("sample")).values_list("date", "count"):
            l = date_to_month_year_string(d)
            samples_by_month[l] += count

        cumulative_samples = np.array(list(samples_by_month.values())).cumsum().tolist()

        tag_context = {
            "labels": list(samples_by_month.keys()),
            "cumulative_samples": cumulative_samples
        }
    return tag_context
