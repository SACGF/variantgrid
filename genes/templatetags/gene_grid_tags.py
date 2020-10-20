from django.conf import settings
from django.db.models.query_utils import Q
from django.template import Library

from annotation.forms import HPOSynonymForm, MIMAliasForm
from genes.forms import GeneListCategoryAutocompleteForm, NamedCustomGeneListForm, PanelAppPanelForm, GeneSymbolForm, \
    GeneAnnotationReleaseForm
from genes.models import GeneInfo, GeneListCategory
from pathtests.forms import ActivePathologyTestForm, SelectPathologyTestVersionForm
from pathtests.models import PathologyTest
from seqauto.forms import EnrichmentKitForm
from seqauto.models import EnrichmentKit
from snpdb.forms import LabSelectForm
from snpdb.models import Company, UserSettings

register = Library()


@register.inclusion_tag("genes/tags/gene_grid_tag.html", takes_context=True)
def gene_grid(context, columns_from_url=None,
              init_callback: str = None,
              update_url_callback: str = None,
              save_gene_list_callback: str = None,
              close_button_callback: str = None,
              data_columns_whitelist: set = None,
              gene_list_categories_whitelist: set = None,
              default_enrichment_kits=None,
              show_gene_annotation_release: bool = True,
              show_custom_gene_form: bool = True,
              show_help: bool = True):
    user = context["user"]
    user_settings = UserSettings.get_for_user(user)

    # Load system defaults
    if default_enrichment_kits is None:
        default_enrichment_kits = settings.GENE_GRID_DEFAULT_ENRICHMENT_KITS
    enrichment_kits = EnrichmentKit.get_enrichment_kits(default_enrichment_kits)
    initial_columns = [f"enrichment-kit-{ek.pk}" for ek in enrichment_kits]
    if columns_from_url:
        initial_columns.extend(columns_from_url.split("/"))

    if show_gene_annotation_release:
        gene_annotation_releases = user_settings.get_gene_annotation_releases()
        initial_columns.extend([f"gene-annotation-release-{release.pk}" for release in gene_annotation_releases])

    data_columns = []
    enrichment_kit_data = {"name": "Enrichment Kit",
                           "description": "Enrichment Kit",
                           "icon_css_class": "enrichment-kit-icon",
                           "form": EnrichmentKitForm()}
    data_columns.append(enrichment_kit_data)
    enrichment_kit_data = {"name": "Gene Annotation Release",
                           "description": "Gene Annotation Release (match genes in VEP annotation)",
                           #"icon_css_class": "",
                           "form": GeneAnnotationReleaseForm()}
    data_columns.append(enrichment_kit_data)

    categories = []

    gene_list_category_filter = None
    company = Company.get_our_company()
    if company:
        if PathologyTest.objects.exists():  # add test categories
            try:
                company_name = str(company).replace("_", " ")
                category = company.genelistcategory
                css_classes = ' '.join([category.icon_css_class, "pathology-test"])
                test_data = {"name": category.name,
                             "icon_css_class": css_classes,
                             "description": f"{company_name} current test",
                             "form": ActivePathologyTestForm(),
                             "form_css_class": "pathology-test-form"}
                categories.append(test_data)
                historical_test_data = {"name": f"historical-{category.name}",
                                        "icon_css_class": css_classes,
                                        "description": f"{company_name} historical test",
                                        "form": SelectPathologyTestVersionForm(),
                                        "form_css_class": "pathology-test-version-form"}
                categories.append(historical_test_data)
            except:
                pass

        gene_list_category_filter = ~Q(company=company)  # Skip as already added

    for category in GeneListCategory.get_gene_list_categories(gene_list_category_filter):
        initial = {"category": category.get("instance")}
        category_id = category.get("pk")
        prefix = f"category-{category_id}"
        category["form"] = GeneListCategoryAutocompleteForm(initial=initial, prefix=prefix)
        category["description"] = category.get("name")
        category["form_css_class"] = "category-gene-list-form"
        categories.append(category)

    panel_app_data = {"name": "Panel App Panel",
                      "description": "Panel App Panel",
                      "icon_css_class": "panel-app-icon",
                      "form": PanelAppPanelForm()}
    categories.append(panel_app_data)

    hpo_category_data = {"name": "Lab Gene Classification Counts",
                         "description": "Lab Gene Classification Counts",
                         "form": LabSelectForm()}
    categories.append(hpo_category_data)

    # HPO
    hpo_category_data = {"name": "HPO",
                         "icon_css_class": "hpo-icon",
                         "description": "Human Phenotype Ontology",
                         "form": HPOSynonymForm(),
                         "form_css_class": "hpo-gene-list-form"}
    categories.append(hpo_category_data)
    # OMIM
    omim_category_data = {"name": "OMIM",
                          "icon_css_class": "omim-icon",
                          "description": "OMIM",
                          "form": MIMAliasForm(),
                          "form_css_class": "omim-gene-list-form"}
    categories.append(omim_category_data)

    # Optionally filter based on whitelist given to tag
    if data_columns_whitelist is not None:
        def in_whitelist(c):
            return c["name"] in data_columns_whitelist
        data_columns = filter(in_whitelist, data_columns)

    if gene_list_categories_whitelist is not None:
        def in_whitelist(c):
            return c["name"] in gene_list_categories_whitelist
        categories = filter(in_whitelist, categories)

    has_gene_info = GeneInfo.objects.filter(gene_list__genelistgenesymbol__isnull=False).exists()
    if show_custom_gene_form:
        named_custom_gene_list_form = NamedCustomGeneListForm(username=user)
    else:
        named_custom_gene_list_form = None
    return {"ENRICHMENT_KIT_COLUMNS": settings.GENE_GRID_ENRICHMENT_KIT_COLUMNS,
            "ENRICHMENT_KIT_COLUMN_TOOL_TIPS": settings.GENE_GRID_ENRICHMENT_KIT_COLUMN_TOOL_TIPS,
            "ENRICHMENT_KIT_COLUMN_LABELS": settings.GENE_GRID_ENRICHMENT_KIT_COLUMN_LABELS,
            "ENRICHMENT_KIT_COLUMN_LABEL_TOOL_TIPS": settings.GENE_GRID_ENRICHMENT_KIT_COLUMN_LABEL_TOOL_TIPS,
            "initial_columns": initial_columns,
            "init_callback": init_callback,
            "update_url_callback": update_url_callback,
            "save_gene_list_callback": save_gene_list_callback,
            "close_button_callback": close_button_callback,
            "data_columns": data_columns,
            "categories": categories,
            "has_gene_info": has_gene_info,
            "gene_symbol_form": GeneSymbolForm(),
            "named_custom_gene_list_form": named_custom_gene_list_form,
            "user": user,
            "show_help": show_help}
