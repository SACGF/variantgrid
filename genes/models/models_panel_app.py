import logging
from datetime import timedelta

from django.conf import settings
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel

from genes.models.models_gene import HGNC
from genes.models.models_gene_list import GeneList, GeneListCategory
from library.guardian_utils import (
    add_public_group_read_permission,
    admin_bot,
)
from snpdb.models import Company
from snpdb.models.models_enums import ImportStatus


class PanelAppServer(models.Model):
    name = models.TextField(unique=True)
    url = models.TextField(unique=True)
    icon_css_class = models.TextField()

    def __str__(self):
        return self.name

    @staticmethod
    def australia_instance() -> 'PanelAppServer':
        return PanelAppServer.objects.get(name="PanelApp Australia")

    @staticmethod
    def england_instance() -> 'PanelAppServer':
        return PanelAppServer.objects.get(name="Genomics England PanelApp")


class PanelAppPanel(TimeStampedModel):
    """ Populated from PanelApp cached web resource task, updated to latest version.
        Soft-deleted when PanelApp returns 404/Not found for the panel — see issue #405. """
    server = models.ForeignKey(PanelAppServer, on_delete=CASCADE)
    panel_id = models.IntegerField()
    disease_group = models.TextField()
    disease_sub_group = models.TextField()
    name = models.TextField()
    # Mirrors PanelApp's own panel status: "public", "internal", "promoted", "retired"
    status = models.TextField()
    current_version = models.TextField()
    # PanelApp doesn't expose a "deleted" status — once a panel is removed upstream the API
    # just starts returning 404 (see issue #405), so we track that locally.
    deleted = models.BooleanField(default=False)

    class Meta:
        unique_together = ('server', 'panel_id')

    @property
    def url(self) -> str:
        return f"{self.server.url}/api/v1/panels/{self.panel_id}"

    @property
    def web_url(self) -> str:
        return f"{self.server.url}/panels/{self.panel_id}/"

    @property
    def cache_valid(self) -> bool:
        # Attempt to use cache if recent and present, otherwise fall through and do a query
        max_age = timedelta(days=settings.PANEL_APP_CACHE_DAYS)
        return timezone.now() < self.modified + max_age

    def __str__(self):
        return self.name


class PanelAppPanelRelevantDisorders(models.Model):
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    name = models.TextField()

    class Meta:
        unique_together = ('panel_app_panel', 'name')


class PanelAppPanelLocalCache(TimeStampedModel):
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    version = models.TextField()

    class Meta:
        unique_together = ("panel_app_panel", "version")

    def get_gene_list(self, panel_app_confidence) -> GeneList:
        from genes.gene_matching import GeneSymbolMatcher

        min_level = int(panel_app_confidence)
        name = f"{self.panel_app_panel.name} v.{self.version}.min_{min_level}"

        # We'll try to re-use gene lists - but it's possible due to race conditions we may occasionally make
        # a duplicate, but this should work fine and won't affect much
        category = GeneListCategory.get_or_create_category(GeneListCategory.PANEL_APP_CACHE, hidden=True)
        gene_list_kwargs = {
            "category": category,
            "name": name,
            "user": admin_bot(),
            "import_status": ImportStatus.SUCCESS,
            "url": self.panel_app_panel.url,
        }
        if gene_list := GeneList.objects.filter(**gene_list_kwargs).order_by("pk").first():
            logging.info("Reused existing gene list: %s", gene_list.pk)
        else:
            gene_list = GeneList.objects.create(**gene_list_kwargs)
            logging.info("Created gene list: %s", gene_list.pk)
            gene_names_list = []
            for pap_lc_gs in self.panelapppanellocalcachegenesymbol_set.select_related("hgnc"):
                confidence_level = int(pap_lc_gs.data["confidence_level"])
                if confidence_level >= min_level:
                    # Use our current approved symbol - PanelApp's own can be from an older HGNC snapshot
                    gene_names_list.append(pap_lc_gs.gene_symbol_str)

            logging.info("Creating symbols: %s", gene_names_list)
            gene_matcher = GeneSymbolMatcher()
            gene_matcher.create_gene_list_gene_symbols(gene_list, gene_names_list)
            gene_list.import_status = ImportStatus.SUCCESS
            gene_list.save()

            # PanelApp gene list should be public
            add_public_group_read_permission(gene_list)

        logging.info("Returning gene list: %s", gene_list.pk)
        return gene_list

    def __str__(self):
        return f"PanelApp cache for {self.panel_app_panel} v{self.version} (mod: {self.modified})"


class PanelAppPanelLocalCacheGeneSymbol(models.Model):
    """ A gene on a cached PanelApp panel.

        PanelApp identifies genes by HGNC ID, and separately reports a gene symbol taken from a dated
        Ensembl/HGNC snapshot - so their symbol can lag the current approved one (Genomics England is
        pinned to Ensembl release 90). PanelApp Australia asks integrators to key off HGNC IDs rather
        than symbols, so 'hgnc' is how we identify the gene and 'gene_symbol_reported' records what
        they called it.

        hgnc is nullable as PanelApp can reference an HGNC ID we have no record of. """
    panel_app_local_cache = models.ForeignKey(PanelAppPanelLocalCache, on_delete=CASCADE)
    # Django names the column 'hgnc_id' - that's HGNC's numeric pk, not the "HGNC:1234" string
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=SET_NULL)
    gene_symbol_reported = models.TextField()
    data = models.JSONField(null=False, blank=True, default=dict)  # API response

    @property
    def gene_symbol_str(self) -> str:
        """ Our current approved symbol where PanelApp's HGNC ID resolves, else the symbol they reported """
        if self.hgnc_id is not None:
            return str(self.hgnc.gene_symbol_id)
        return self.gene_symbol_reported

    def __str__(self):
        return self.gene_symbol_str


class CachedThirdPartyGeneList(models.Model):
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)
    company = models.ForeignKey(Company, on_delete=CASCADE)
