import types
from dataclasses import dataclass
from datetime import datetime
from typing import Optional, Union

from auditlog.models import LogEntry
from django.contrib.auth.models import Group, User
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models import StringAgg, TextField, Value
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.query_utils import Q
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django.utils import timezone
from django.utils.timezone import localtime
from guardian.shortcuts import get_objects_for_user
from model_utils.models import TimeStampedModel

from genes.models.models_gene import Gene, GeneSymbol, GeneSymbolAlias
from genes.models.models_gene_annotation_release import GeneAnnotationRelease, ReleaseGeneSymbol
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import (
    DjangoPermission,
    assign_permission_to_user_and_groups,
)
from snpdb.enums import ImportStatus
from snpdb.models import Company, Sample, Wiki


class GeneListCategory(models.Model):
    # This was a separate table rather than hard coded choices so that we can use dependent fields / chaining forms.
    # Used as CSS classes so have to have no spaces
    NODE_CUSTOM_TEXT = 'NodeCustomText'
    SAMPLE_GENE_LIST = 'SampleGeneList'
    QC_COVERAGE_CUSTOM_TEXT = 'QCCoverageCustomText'
    PATHOLOGY_TEST = 'PathologyTest'
    PATHOLOGY_TEST_ORDER = 'PathologyTestOrder'
    PANEL_APP_CACHE = 'PanelAppCache'

    name = models.TextField()
    company = models.OneToOneField(Company, null=True, blank=True, on_delete=CASCADE)
    icon_css_class = models.TextField(blank=True)
    hidden = models.BooleanField(default=False)  # for special ones
    public = models.BooleanField(default=False)
    description = models.TextField()

    @staticmethod
    def get_or_create_category(category_name, hidden=False):
        category, created = GeneListCategory.objects.get_or_create(name=category_name)
        if created:
            category.hidden = hidden
            category.save()
        return category

    @staticmethod
    def _get_pathology_test_gene_category(category_name, set_company=False):
        """ Requires settings.COMPANY to be set """
        category = None
        company = Company.get_our_company()
        if company:
            category, created = GeneListCategory.objects.get_or_create(name=category_name)
            if created:
                category.hidden = True
                if set_company:  # 1-to-1 field so can't always do it
                    category.company = company
                category.icon_css_class = "company-" + company.name.lower() + "-icon"
                category.save()
        return category

    @staticmethod
    def get_pathology_test_gene_category():
        return GeneListCategory._get_pathology_test_gene_category(GeneListCategory.PATHOLOGY_TEST, set_company=True)

    @staticmethod
    def get_pathology_test_order_gene_category():
        return GeneListCategory._get_pathology_test_gene_category(GeneListCategory.PATHOLOGY_TEST_ORDER)

    @staticmethod
    def get_gene_list_categories(extra_filters=None):
        categories = [{"name": "User", "pk": None}]  # Blank / no category
        qs = GeneListCategory.objects.exclude(genelist__isnull=True).exclude(hidden=True)
        if extra_filters:
            qs = qs.filter(extra_filters)
        for glc in qs.order_by("name"):
            categories.append({"name": glc.name,
                               "pk": glc.pk,
                               "icon_css_class": glc.icon_css_class,
                               "instance": glc})
        return categories

    def __str__(self):
        return self.name


class GeneList(GuardianPermissionsMixin, TimeStampedModel):
    """ Stores a gene/transcript list (to be used as a filter) """

    category = models.ForeignKey(GeneListCategory, null=True, blank=True, on_delete=CASCADE)
    name = models.TextField()
    user = models.ForeignKey(User, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    error_message = models.TextField(null=True, blank=True)
    locked = models.BooleanField(default=False)
    url = models.TextField(null=True, blank=True)
    # The list this one's genes were copied from - for a list created from PanelApp, the PanelApp cache's list
    # for that panel version. Edits since the copy are in the history (@see get_source) #2846
    source_gene_list = models.ForeignKey('self', null=True, blank=True, on_delete=SET_NULL,
                                         related_name="copied_gene_lists")

    @staticmethod
    def get_gene_ids_for_gene_lists(release: GeneAnnotationRelease, gene_lists: list['GeneList']):
        """ For GeneList node, we need to get query for multiple lists - and it's much faster to build the merged
            query here, than via joining separate queryset """
        # Route the gene-list symbols through a subquery rather than joining GeneListGeneSymbol via GeneSymbol.
        # The join-through-GeneSymbol path makes Postgres scan every ReleaseGeneSymbol for the release (~47k rows)
        # and hash-join; the subquery lets it use the (release_id, gene_symbol_id) unique index instead.
        symbols_qs = GeneListGeneSymbol.objects.filter(gene_list__in=gene_lists).values_list("gene_symbol_id", flat=True)
        rgs_qs = ReleaseGeneSymbol.objects.filter(release=release, gene_symbol__in=symbols_qs)
        return rgs_qs.values_list("releasegenesymbolgene__gene", flat=True)

    def get_genes(self, release: GeneAnnotationRelease):
        """ Get Genes (from a release) for symbols in this gene list """
        # Route symbols through a subquery rather than joining GeneListGeneSymbol via GeneSymbol, so Postgres
        # uses the ReleaseGeneSymbol (release_id, gene_symbol_id) unique index instead of scanning the release.
        symbols_qs = self.genelistgenesymbol_set.values_list("gene_symbol_id", flat=True)
        rgs_qs = ReleaseGeneSymbol.objects.filter(release=release, gene_symbol__in=symbols_qs)
        return Gene.objects.filter(releasegenesymbolgene__release_gene_symbol__in=rgs_qs)

    def get_gene_names(self):
        return self.genelistgenesymbol_set.filter(gene_symbol__isnull=False).values_list("gene_symbol", flat=True).distinct()

    def _copy_gene_symbols_to(self, gene_list: 'GeneList'):
        gene_list_gene_symbols = list(self.genelistgenesymbol_set.all())
        for glgs in gene_list_gene_symbols:
            glgs.pk = None
            glgs.gene_list = gene_list
        GeneListGeneSymbol.objects.bulk_create(gene_list_gene_symbols)

    def _copy_history_to(self, gene_list: 'GeneList'):
        log_entries = list(self._modification_log_entries())
        for log_entry in log_entries:
            log_entry.pk = None
            log_entry.object_pk = str(gene_list.pk)
            log_entry.object_id = gene_list.pk
        LogEntry.objects.bulk_create(log_entries)

    def clone(self):
        """ Needed to clone analysis node GeneListNode. Leaves self as the original """
        copy = GeneList.objects.get(pk=self.pk)
        copy.pk = None
        copy._state.adding = True
        copy.save()
        self._copy_gene_symbols_to(copy)
        self._copy_history_to(copy)
        return copy

    def copy_for_user(self, user, name: Optional[str] = None, category: Optional[GeneListCategory] = None) -> 'GeneList':
        """ An unlocked copy owned by user - unlike clone(), which keeps the owner, lock and category, and is
            used to snapshot a list for an analysis node or pathology test.
            The copy's source is where this list's genes came from: this list's own source along with the edits
            made to it since, otherwise this list """
        copy = GeneList.objects.create(name=name or f"Copy of {self.name}", user=user, category=category,
                                       import_status=self.import_status,
                                       source_gene_list_id=self.source_gene_list_id or self.pk)
        self._copy_gene_symbols_to(copy)
        if self.source_gene_list_id:
            self._copy_history_to(copy)
        return copy

    def _modification_log_entries(self):
        """ LogEntry is shared with everything auditlog tracks, so only the entries log_modification wrote """
        return LogEntry.objects.get_for_object(self).filter(additional_data__has_key="added").order_by("timestamp", "pk")

    def log_modification(self, user: Optional[User], added, removed):
        """ Record genes added to / removed from the list after it was created. Every edit goes through
            add_and_remove_gene_symbols or genes/gene_matching.py:GeneSymbolMatcher.set_gene_list_gene_symbols_from_text,
            which call this - so a list with no modifications is as it was created.
            GeneListGeneSymbol isn't registered with auditlog: bulk_create fires no signals for it to see, and
            an edit is one entry rather than one per gene """
        if added or removed:
            LogEntry.objects.log_create(self, force_log=True, action=LogEntry.Action.UPDATE, actor=user,
                                        additional_data={"added": sorted(added), "removed": sorted(removed)})

    def get_modifications(self) -> list['GeneListModification']:
        """ Edits since the list was created, oldest first """
        return [GeneListModification(user=log_entry.actor, date=log_entry.timestamp,
                                     added=log_entry.additional_data["added"],
                                     removed=log_entry.additional_data["removed"])
                for log_entry in self._modification_log_entries().select_related("actor")]

    def get_source(self) -> Optional['GeneListSource']:
        """ Where a copied list's genes came from, and how it has been edited since """
        if not self.source_gene_list_id:
            return None
        return GeneListSource(gene_list=self.source_gene_list, modifications=self.get_modifications())

    def save(self, *args, **kwargs):
        assign_permissions = kwargs.pop("assign_permissions", False)
        initial_save = not self.pk
        super().save(*args, **kwargs)
        if initial_save or assign_permissions:
            # logging.info("GeneList: assign_permission_to_user_and_groups")
            assign_permission_to_user_and_groups(self.user, self)

    # can_view() uses Guardian object-level perms (provided by GuardianPermissionsMixin).

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        # As mixin's can_write(), but a locked GeneList can't be written to by anyone
        write_perm = DjangoPermission.perm(self, DjangoPermission.WRITE)
        return user_or_group.has_perm(write_perm, self) and not self.locked

    @classmethod
    def filter_writable_for_user(cls, user, pks=None):
        """ Batch can_write - a locked GeneList is read only """
        return super().filter_writable_for_user(user, pks=pks).filter(locked=False)

    @classmethod
    def allow_group_permission_delete(cls) -> bool:
        return True  # User-created list; deletable via the group_permissions delete view

    def get_warnings(self, release: GeneAnnotationRelease) -> list[str]:
        counts = {"unmatched symbols": self.unmatched_gene_symbols.count(),
                  "aliased": self.aliased_genes.count(),
                  "unmatched genes": self.unmatched_genes(release).count()}

        warnings = []
        for name, num in counts.items():
            if num:
                warnings.append(f"{name} x {num}")
        return warnings

    @property
    def unmatched_gene_symbols(self):
        return self.genelistgenesymbol_set.filter(gene_symbol__isnull=True).order_by("original_name")

    @property
    def aliased_genes(self):
        return self.genelistgenesymbol_set.filter(gene_symbol_alias__isnull=False)

    def unmatched_genes(self, release: GeneAnnotationRelease):
        string_agg = GeneListGeneSymbol.get_joined_genes_qs_annotation_for_release(release)
        return self.genelistgenesymbol_set.annotate(matched_genes=string_agg).filter(matched_genes__isnull=True)

    def add_and_remove_gene_symbols(self, gene_symbol_additions, gene_symbol_deletions,
                                    gene_additions_modification_info=None, user: Optional[User] = None):
        """ Either works for all, or fails and makes no modification
            returns (num_added, num_deleted) """
        if gene_additions_modification_info is None:
            gene_additions_modification_info = {}

        if self.locked:
            msg = "Can't modify a locked GeneList!"
            raise PermissionDenied(msg)

        added = []
        for gene_symbol in gene_symbol_additions:
            modification_info = gene_additions_modification_info.get(gene_symbol)
            defaults = {"original_name": gene_symbol, "modification_info": modification_info}
            _, created = GeneListGeneSymbol.objects.get_or_create(gene_list=self, gene_symbol_id=gene_symbol,
                                                                  defaults=defaults)
            if created:
                added.append(gene_symbol)

        # Match original name and symbol so we can delete non-matched
        name_or_symbol_match = Q(original_name__in=gene_symbol_deletions) | Q(gene_symbol__in=gene_symbol_deletions)
        qs = self.genelistgenesymbol_set.filter(name_or_symbol_match)
        removed = [original_name or gene_symbol for original_name, gene_symbol
                   in qs.values_list("original_name", "gene_symbol_id")]
        qs.delete()

        if added or removed:
            self.set_modified_to_now()
            self.log_modification(user, added, removed)

        return len(added), len(removed)

    @staticmethod
    def get_for_user(user, gene_list_id, success_only=True):
        try:
            return GeneList.filter_for_user(user, success_only).get(pk=gene_list_id)
        except GeneList.DoesNotExist as exc:
            # Need to distinguish between does not exist and no permission
            get_object_or_404(GeneList, pk=gene_list_id)  # potentially throws GeneList.DoesNotExist
            # If we're here, object exists but we have a permission error
            msg = f"You don't have permission to access gene_list {gene_list_id}"
            raise PermissionDenied(msg) from exc

    @staticmethod
    def filter_for_user(user, success_only=True):
        perm = DjangoPermission.perm(GeneList, DjangoPermission.READ)
        user_qs = get_objects_for_user(user, perm, klass=GeneList, accept_global_perms=False)

        # Sample gene lists for samples we have permission to see
        samples_qs = Sample.filter_for_user(user)
        qc_sample_q = Q(customtextgenelist__qcgenelist__qc__sequencing_sample__samplefromsequencingsample__sample__in=samples_qs)
        sample_gene_list_q = Q(samplegenelist__sample__in=samples_qs)
        sample_gene_list_qs = GeneList.objects.filter(qc_sample_q | sample_gene_list_q,
                                                      category__name=GeneListCategory.SAMPLE_GENE_LIST)

        qs = user_qs | sample_gene_list_qs
        if success_only:
            qs = qs.filter(import_status=ImportStatus.SUCCESS)
        return qs.distinct()  # Sometimes gene list can be reached multiple ways via 'or' above

    @staticmethod
    def _visible_gene_lists(gene_lists_qs):
        hidden_categories_qs = GeneListCategory.objects.filter(hidden=True)
        return gene_lists_qs.exclude(category__in=hidden_categories_qs)  # Allow through no category

    @classmethod
    def visible_gene_lists_containing_gene_symbol(cls, gene_lists_qs, gene_symbol):
        qs = cls._visible_gene_lists(gene_lists_qs)
        gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
        return qs.filter(genelistgenesymbol__gene_symbol=gene_symbol)

    def set_modified_to_now(self):
        GeneList.objects.filter(pk=self.pk).update(modified=timezone.now())

    @property
    def name_and_size(self) -> str:
        return f"{self.name} ({self.genelistgenesymbol_set.count()} x genes)"

    def __str__(self):
        if source := self.get_source():
            return f"{self.name_and_size} - from {source}"
        return self.name_and_size

    def get_absolute_url(self):
        return reverse('view_gene_list', args=[str(self.id)])

    @classmethod
    def get_listing_url(cls):
        return reverse('gene_lists')


@dataclass(frozen=True)
class GeneListModification:
    """ One edit of a gene list @see GeneList.log_modification """
    user: Optional[User]
    date: datetime
    added: list[str]
    removed: list[str]


@dataclass(frozen=True)
class GeneListSource:
    """ @see GeneList.get_source """
    MAX_NAMED_CHANGES = 5

    gene_list: GeneList
    modifications: list[GeneListModification]

    @property
    def panel_app_version(self):
        """ The PanelApp panel version and confidence a PanelApp cache list holds, else None
            (PanelAppPanelLocalCacheGeneList - models_panel_app imports this module) """
        return getattr(self.gene_list, "panel_app_version", None)

    @property
    def description(self) -> str:
        if panel_app_version := self.panel_app_version:
            return str(panel_app_version)
        return f"gene list '{self.gene_list.name}'"

    @property
    def net_changes(self) -> tuple[list[str], list[str]]:
        """ (added, removed) against the source - a gene added then removed again is in neither """
        added = set()
        removed = set()
        for modification in self.modifications:
            for name in modification.added:
                if name in removed:
                    removed.discard(name)
                else:
                    added.add(name)
            for name in modification.removed:
                if name in added:
                    added.discard(name)
                else:
                    removed.add(name)
        return sorted(added), sorted(removed)

    @property
    def label(self) -> str:
        """ The source first, so the common "PanelApp panel plus a gene" reads as that """
        added, removed = self.net_changes
        if len(added) + len(removed) <= self.MAX_NAMED_CHANGES:
            changes = [f"{sign} {', '.join(names)}" for sign, names in (("+", added), ("-", removed)) if names]
        else:
            changes = [f"{sign} {len(names)} {verb}" for sign, names, verb
                       in (("+", added, "added"), ("-", removed, "removed")) if names]
        return " ".join([self.description, ", ".join(changes)]).strip()

    def __str__(self):
        return self.label


def create_fake_gene_list(*args, **kwargs):
    """ An unsaved GeneList with no URL. A patched instance rather than an abstract subclass, as Django won't
        instantiate an abstract model """

    def get_absolute_url(self):
        return None

    gene_list = GeneList(*args, **kwargs)
    gene_list.get_absolute_url = types.MethodType(get_absolute_url, gene_list)
    return gene_list


class GeneListGeneSymbol(models.Model):
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    original_name = models.TextField(null=True, blank=True)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    gene_symbol_alias = models.ForeignKey(GeneSymbolAlias, null=True, on_delete=CASCADE)
    modification_info = models.TextField(null=True)

    class Meta:
        unique_together = ('gene_list', 'original_name')

    @staticmethod
    def get_joined_genes_qs_annotation_for_release(release: GeneAnnotationRelease):
        """ Used to annotate GeneListGeneSymbol queryset """
        return StringAgg("gene_symbol__releasegenesymbol__releasegenesymbolgene__gene",
                         delimiter=Value(','), distinct=True, output_field=TextField(),
                         filter=Q(gene_symbol__releasegenesymbol__release=release))

    def __str__(self):
        if self.gene_symbol:
            description = f"{self.gene_symbol_id}"
            if self.gene_symbol_alias:
                description += f" (via {self.gene_symbol_alias})"
        else:
            description = f"'{self.original_name}' (unrecognised symbol)"
        return description


class CustomTextGeneList(models.Model):
    """' Some human entered text which gets pulled apart to create a gene list """
    sha256_hash = models.TextField()
    name = models.TextField()
    text = models.TextField()
    gene_list = models.OneToOneField(GeneList, null=True, on_delete=SET_NULL)

    def clone(self):
        gene_list = self.gene_list.clone() if self.gene_list else None
        return CustomTextGeneList.objects.create(sha256_hash=self.sha256_hash, name=self.name, text=self.text,
                                                 gene_list=gene_list)

    def __str__(self):
        text = self.text
        if len(text) > 50:
            text = text[:50] + " ..."
        return f"CustomTextGeneList for '{text}': {self.gene_list}"


class SampleGeneList(TimeStampedModel):
    """ A sample's genes of interest. There can be multiple SampleGeneLists per sample, but only 1 active one
        @see sample_gene_list_created for which becomes active """
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    visible = models.BooleanField(default=True, blank=False)

    class Meta:
        unique_together = ('sample', 'gene_list')
        ordering = ['created']

    @staticmethod
    def filter_for_user(user):
        return SampleGeneList.objects.filter(sample__in=Sample.filter_for_user(user))

    @property
    def from_pipeline(self) -> bool:
        """ Submitted through the API with the sequencing run (a QCGeneList), not created or applied by a user -
            who is gene_list.user """
        return GeneList.objects.filter(pk=self.gene_list_id, customtextgenelist__qcgenelist__isnull=False).exists()

    def __str__(self):
        s = f"{self.sample}: {localtime(self.modified)}"
        if not self.visible:
            s += " (hidden)"
        return s


@receiver(post_save, sender=SampleGeneList)
def sample_gene_list_created(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    """ A new GOI becomes the active one, except that a pipeline one doesn't replace a user's correction.
        Switching updates the existing ActiveSampleGeneList, so only the sample's first one is 'created' -
        the signal that auto-runs analysis templates requiring a sample gene list """
    if created:
        sample = instance.sample
        active = ActiveSampleGeneList.objects.filter(sample=sample).select_related("sample_gene_list").first()
        if active is None:
            ActiveSampleGeneList.objects.get_or_create(sample=sample, defaults={"sample_gene_list": instance})
        elif not instance.from_pipeline or active.sample_gene_list.from_pipeline:
            active.sample_gene_list = instance
            active.save()


class ActiveSampleGeneList(TimeStampedModel):
    """ Use 1-to-1 to enforce there's only 1 in DB
        (as compared to an "active" flag on SampleGeneList) """
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    sample_gene_list = models.ForeignKey(SampleGeneList, on_delete=CASCADE)


class GeneListWiki(Wiki):
    gene_list = models.OneToOneField(GeneList, on_delete=CASCADE)

    def _get_restricted_object(self):
        return self.gene_list
