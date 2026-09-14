"""
The named family structures over a Cohort. Trio (mother/father/proband), Quad (Trio plus a sibling)
and Duo (proband plus one relative - a parent or a sibling) each point at CohortSamples of one Cohort,
which is where their permissions and genome build come from; FamilyGroupMixin holds what they share.
Each provides get_cohort_samples() in pedigree order - that is what the analysis inheritance nodes and
the pedigree figures read.
"""
from typing import Optional

from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel

from library.django_utils import SortByPKMixin
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from library.preview_request import PreviewModelMixin, SvgSymbolPreviewIconMixin
from patients.models_enums import Sex
from snpdb.models.models_cohort import Cohort, CohortSample
from snpdb.models.models_enums import DuoRelationship
from snpdb.models.models_vcf import Sample


def _member_details(member, affected: bool) -> str:
    return f"{member} ({'affected' if affected else 'unaffected'})"


class FamilyGroupMixin:
    """ Shared by Duo, Trio and Quad - permissions and display that don't care how many members
        there are. Subclasses provide get_cohort_samples(), pedigree_icon_members and their own urls. """
    pk: int
    name: str
    cohort: Cohort
    proband: CohortSample
    proband_sex: Optional[str]
    pedigree_icon_members: tuple[str, ...]

    @classmethod
    def get_permission_class(cls):
        return Cohort

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-people-roof"

    def get_preview_icon_css_class(self) -> str:
        """ Blacken the affected members, the way the pedigree figure on the view page does """
        return " ".join(f"{member}-affected" for member in self.pedigree_icon_members
                        if getattr(self, f"{member}_affected"))

    @property
    def preview(self) -> 'PreviewData':
        return self.preview_with(identifier=str(self))

    def get_permission_object(self):
        # Permissions are based on the cohort
        return self.cohort

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(cohort__in=queryset)

    @property
    def genome_build(self):
        return self.cohort.genome_build

    @property
    def data_archived(self) -> bool:
        return self.cohort.data_archived

    def get_cohort_samples(self) -> list[CohortSample]:
        """ The members in pedigree order - what the inheritance nodes and pedigree figures read """
        raise NotImplementedError()

    def get_samples(self):
        return Sample.objects.filter(cohortsample__in=self.get_cohort_samples()).order_by("pk")

    @property
    def effective_proband_sex(self) -> Sex:
        """ What the scientist chose in the wizard, otherwise what the patient record says """
        if self.proband_sex:
            return Sex(self.proband_sex)
        return self.proband.sample.patient_sex

    def __str__(self):
        return self.name or f"{type(self).__name__} {self.pk}"


class ParentsMixin:
    """ Trio and Quad - both parents are present. A Duo has one relative that may be neither. """
    mother: CohortSample
    mother_affected: bool
    father: CohortSample
    father_affected: bool

    pedigree_icon_members = ("mother", "father")  # Quad adds the sibling

    @property
    def mother_details(self):
        return _member_details(self.mother, self.mother_affected)

    @property
    def father_details(self):
        return _member_details(self.father, self.father_affected)


class Trio(ParentsMixin, FamilyGroupMixin, GuardianPermissionsAutoInitialSaveMixin, SvgSymbolPreviewIconMixin, PreviewModelMixin,
           SortByPKMixin, TimeStampedModel):
    """ A simple pedigree used frequently for Mendellian disease (TrioNode in analysis)
        and karyomapping """
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    mother = models.ForeignKey(CohortSample, related_name='trio_mother', on_delete=CASCADE)
    mother_affected = models.BooleanField(default=False)
    father = models.ForeignKey(CohortSample, related_name='trio_father', on_delete=CASCADE)
    father_affected = models.BooleanField(default=False)
    proband = models.ForeignKey(CohortSample, related_name='trio_proband', on_delete=CASCADE)
    # Set in the trio wizard when the scientist resolves patient.sex vs sample.detected_sex, eg a male
    # fetus in a prenatal case entered under the mother's record. Null = go by the patient record
    proband_sex = models.CharField(max_length=1, choices=Sex.choices, null=True, blank=True)

    preview_icon_symbol = "node-icon-trio"  # TrioNode wears this too - see get_node_class_icon

    @classmethod
    def preview_if_url_visible(cls) -> str:
        return "trios"

    def get_cohort_samples(self):
        return [self.mother, self.father, self.proband]

    def get_absolute_url(self):
        return reverse('view_trio', kwargs={"pk": self.pk})

    def get_listing_url(self):
        return reverse('trios')


class Quad(ParentsMixin, FamilyGroupMixin, GuardianPermissionsAutoInitialSaveMixin, SvgSymbolPreviewIconMixin, PreviewModelMixin,
           SortByPKMixin, TimeStampedModel):
    """Mother + Father + Proband + Sibling.

    Extends the Trio concept to 4 family members. The sibling (typically
    unaffected) narrows down candidate variants because they share the same
    parental genome without sharing the proband's phenotype.
    """
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    mother = models.ForeignKey(CohortSample, related_name='quad_mother', on_delete=CASCADE)
    mother_affected = models.BooleanField(default=False)
    father = models.ForeignKey(CohortSample, related_name='quad_father', on_delete=CASCADE)
    father_affected = models.BooleanField(default=False)
    proband = models.ForeignKey(CohortSample, related_name='quad_proband', on_delete=CASCADE)
    sibling = models.ForeignKey(CohortSample, related_name='quad_sibling', on_delete=CASCADE)
    # Set in the quad wizard when the scientist resolves patient.sex vs sample.detected_sex
    proband_sex = models.CharField(max_length=1, choices=Sex.choices, null=True, blank=True)
    sibling_affected = models.BooleanField(default=False)

    preview_icon_symbol = "node-icon-quad"  # QuadNode wears this too - see get_node_class_icon
    pedigree_icon_members = ("mother", "father", "sibling")

    @classmethod
    def preview_if_url_visible(cls) -> str:
        return "quads"

    def get_cohort_samples(self):
        return [self.mother, self.father, self.proband, self.sibling]

    def get_absolute_url(self):
        return reverse('view_quad', kwargs={"pk": self.pk})

    def get_listing_url(self):
        return reverse('quads')

    @property
    def sibling_details(self):
        return _member_details(self.sibling, self.sibling_affected)


class Duo(FamilyGroupMixin, GuardianPermissionsAutoInitialSaveMixin, SvgSymbolPreviewIconMixin, PreviewModelMixin,
          SortByPKMixin, TimeStampedModel):
    """Proband + one relative - a parent, or a sibling (#1861).

    The rest of the family is unavailable (deceased, not consented, cost, or a prenatal case entered
    under the mother's record) so a Trio can't be made. One `relative` FK that is always set, plus
    the `relationship` the inheritance modes need - X-linked recessive is only meaningful through the
    mother, comp het is half phased on "one from the parent, one not", and a sibling pair carries no
    parental transmission at all, so the parent-only modes have nothing to filter on.
    """
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    proband = models.ForeignKey(CohortSample, related_name='duo_proband', on_delete=CASCADE)
    relative = models.ForeignKey(CohortSample, related_name='duo_relative', on_delete=CASCADE)
    relationship = models.CharField(max_length=1, choices=DuoRelationship.choices)
    relative_affected = models.BooleanField(default=False)
    # Set in the duo wizard when the scientist resolves patient.sex vs sample.detected_sex
    proband_sex = models.CharField(max_length=1, choices=Sex.choices, null=True, blank=True)

    preview_icon_symbol = "node-icon-duo"  # DuoNode wears this too - see get_node_class_icon
    pedigree_icon_members = ("relative",)

    @classmethod
    def preview_if_url_visible(cls) -> str:
        return "duos"

    def get_cohort_samples(self):
        return [self.relative, self.proband]

    def get_absolute_url(self):
        return reverse('view_duo', kwargs={"pk": self.pk})

    def get_listing_url(self):
        return reverse('duos')

    @property
    def relative_is_sibling(self) -> bool:
        return self.relationship == DuoRelationship.SIBLING

    @property
    def parent_is_mother(self) -> bool:
        return self.relationship == DuoRelationship.MOTHER

    @property
    def relationship_label(self) -> str:
        return DuoRelationship(self.relationship).label

    @property
    def missing_parent_label(self) -> str:
        """ The parent we don't have - named in the "absent in parent" warning. Only meaningful for a
            parent duo: a sibling duo is missing both, so callers ask relative_is_sibling first """
        if self.parent_is_mother:
            return DuoRelationship.FATHER.label
        return DuoRelationship.MOTHER.label

    def get_preview_icon_css_class(self) -> str:
        """ The symbol draws all three relative shapes - the relationship class picks which shows """
        css_classes = [f"duo-{self.relationship_label.lower()}"]
        if affected := super().get_preview_icon_css_class():
            css_classes.append(affected)
        return " ".join(css_classes)

    @property
    def relative_details(self):
        return _member_details(self.relative, self.relative_affected)
