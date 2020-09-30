from django.template.library import Library

from snpdb.models import Company

register = Library()

@register.inclusion_tag("snpdb/tags/company_logo.html", takes_context=True)
def company_logo(context):
    css_classes = None
    company = Company.get_our_company()
    if company:
        category = company.genelistcategory
        css_classes = ' '.join([category.icon_css_class, "pathology-test"])

    return {"company_css_class": css_classes}
