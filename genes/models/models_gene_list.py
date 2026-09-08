import types
from typing import Union

from django.contrib.auth.models import Group, User
from django.core.exceptions import PermissionDenied
from django.db import IntegrityError, models, transaction
from django.db.models import StringAgg, TextField, Value
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.query_utils import Q
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django.utils import timezone
from django.utils.timezone import localtime
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user

from genes.models.models_gene import Gene, GeneSymbol, GeneSymbolAlias
from genes.models.models_gene_annotation_release import GeneAnnotationRelease, ReleaseGeneSymbol
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import (
    DjangoPermission,
    assign_permission_to_user_and_groups,
)
from snpdb.models import Company, Sample, Wiki
from snpdb.models.models_enums import ImportStatus


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

    def get_q(self, variant_annotation_version):
        """ For a Variant queryset """
        from annotation.models.models import VariantTranscriptAnnotation
        genes_qs = self.get_genes(variant_annotation_version.gene_annotation_release)
        return VariantTranscriptAnnotation.get_overlapping_genes_q(variant_annotation_version,
                                                                   genes_qs)

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

    def clone(self):
        """ Needed to clone analysis node GeneListNode """
        genelistgenesymbol_set = set(self.genelistgenesymbol_set.all())
        copy = self
        copy.pk = None
        copy.save()

        genes = 0
        for glgs in genelistgenesymbol_set:
            glgs.pk = None
            glgs.gene_list = copy
            glgs.save()
            genes += 1

        return copy

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
                                    gene_additions_modification_info=None):
        """ Either works for all, or fails and makes no modification
            returns (num_added, num_deleted) """
        if gene_additions_modification_info is None:
            gene_additions_modification_info = {}

        if self.locked:
            msg = "Can't modify a locked GeneList!"
            raise PermissionDenied(msg)

        num_added = 0

        for gene_symbol in gene_symbol_additions:
            modification_info = gene_additions_modification_info.get(gene_symbol)
            defaults = {"original_name": gene_symbol, "modification_info": modification_info}
            _, created = GeneListGeneSymbol.objects.get_or_create(gene_list=self, gene_symbol_id=gene_symbol,
                                                                  defaults=defaults)
            if created:
                num_added += 1

        # Match original name and symbol so we can delete non-matched
        name_or_symbol_match = Q(original_name__in=gene_symbol_deletions) | Q(gene_symbol__in=gene_symbol_deletions)
        qs = self.genelistgenesymbol_set.filter(name_or_symbol_match)
        num_deleted, _ = qs.delete()

        if num_added or num_deleted:
            self.set_modified_to_now()

        return num_added, num_deleted

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
        sample_gene_list_qs = GeneList.objects.filter(category__name=GeneListCategory.SAMPLE_GENE_LIST,
                                                      customtextgenelist__qcgenelist__qc__bam_file__sequencing_sample__samplefromsequencingsample__sample__in=samples_qs)

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

    def __str__(self):
        return f"{self.name} ({self.genelistgenesymbol_set.count()} x genes)"

    def get_absolute_url(self):
        url = self.url
        if url is None:
            url = reverse('view_gene_list', args=[str(self.id)])
        return url

    @classmethod
    def get_listing_url(cls):
        return reverse('gene_lists')


def create_fake_gene_list(*args, **kwargs):
    """  Originally FakeGeneList had abstract=True but with Django 3.2 got "Abstract models cannot be instantiated" """

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
        copy = self
        old_gene_list = self.gene_list
        if self.gene_list:
            copy.gene_list = self.gene_list.clone()
        copy.pk = None
        copy.save()
        self.gene_list = old_gene_list
        return copy

    def __str__(self):
        text = self.text
        if len(text) > 50:
            text = text[:50] + " ..."
        return f"CustomTextGeneList for '{text}': {self.gene_list}"


class SampleGeneList(TimeStampedModel):
    """ There can be multiple SampleGeneLists per sample, but only 1 active one.
        If multiple exist, the active one must be set manually """
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    visible = models.BooleanField(default=True, blank=False)

    class Meta:
        unique_together = ('sample', 'gene_list')
        ordering = ['created']

    def __str__(self):
        s = f"{self.sample}: {localtime(self.modified)}"
        if not self.visible:
            s += " (hidden)"
        return s


@receiver(post_save, sender=SampleGeneList)
def sample_gene_list_created(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    if created:
        sample = instance.sample
        if SampleGeneList.objects.filter(sample=sample).count() > 1:
            # Multiple exist, so need to set manually
            ActiveSampleGeneList.objects.filter(sample=sample).delete()
        else:
            try:
                with transaction.atomic():
                    # There can only be 1 - if this works it's active
                    # This will also trigger an event which is used to auto lauch analyses
                    ActiveSampleGeneList.objects.create(sample=sample, sample_gene_list=instance)
            except IntegrityError:
                ActiveSampleGeneList.objects.filter(sample=sample).delete()


class ActiveSampleGeneList(TimeStampedModel):
    """ Use 1-to-1 to enforce there's only 1 in DB
        (as compared to an "active" flag on SampleGeneList) """
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    sample_gene_list = models.ForeignKey(SampleGeneList, on_delete=CASCADE)


class GeneListWiki(Wiki):
    gene_list = models.OneToOneField(GeneList, on_delete=CASCADE)

    def _get_restricted_object(self):
        return self.gene_list
