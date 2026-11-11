""" The Duo/Trio/Quad wizards - turn a cohort's samples into a family group.

    Every wizard asks the same questions (which sample has which role, who's affected, and which sex
    to trust for the proband) so they share FamilyWizardView and analysis/family_wizard.html - a
    subclass names the family, its form and how to make the object from the roles.
"""
import re

from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.shortcuts import redirect
from django.views.generic import FormView

from analysis.enums import DuoSample, QuadSample, TrioSample
from analysis.forms import UserDuoWizardForm, UserQuadWizardForm, UserTrioWizardForm
from library.utils.collection_utils import remove_duplicates_from_list
from library.utils.text_utils import longest_common_substring
from patients.enums import Sex
from snpdb.enums import ImportStatus
from snpdb.models import Cohort, Duo, Quad, Sample, Trio

# Shorter than this, what the sample names share is a run or lane fragment ("_S", "L00") rather than a family
FAMILY_NAME_MIN_COMMON_LENGTH = 4
VCF_FILE_EXTENSIONS_RE = re.compile(r"(\.(vcf|bcf|gz|bgz))+$", re.IGNORECASE)
NAME_SEPARATORS = "_-. "


def family_name_suggestions(samples: list[Sample]) -> list[str]:
    """ Names for the family better than the roles-and-cohort default: the VCF's, when it was called on
        just these samples, then the whole words the sample names have in common (#2760) """
    suggestions = []
    vcfs = {sample.vcf for sample in samples}
    if len(vcfs) == 1:
        vcf = vcfs.pop()
        if vcf.name and vcf.sample_set.count() == len(samples):
            if vcf_name := VCF_FILE_EXTENSIONS_RE.sub("", vcf.name):
                suggestions.append(vcf_name)
    common = longest_common_substring([sample.name for sample in samples],
                                      separators=NAME_SEPARATORS).strip(NAME_SEPARATORS)
    if len(common) >= FAMILY_NAME_MIN_COMMON_LENGTH:
        suggestions.append(common)
    return remove_duplicates_from_list(suggestions)


def _confident_sex(sample: Sample) -> Sex:
    """ The sex the wizard works from: the patient record, or the chrX call when the record has none.
        The record wins a disagreement - it is the one a person can correct """
    patient_sex = sample.patient_sex
    if patient_sex != Sex.UNKNOWN:
        return patient_sex
    return sample.detected_sex


def _sample_sexes(samples: list[Sample]) -> list[dict]:
    """ The proband and the roles on offer are picked client side, so hand the JS every sample's sexes """
    return [{"patient_sex": s.patient_sex.value, "patient_sex_label": s.patient_sex.label,
             "detected_sex": s.detected_sex.value, "detected_sex_label": s.detected_sex.label,
             "sex": _confident_sex(s).value, "sex_label": _confident_sex(s).label}
            for s in samples]


def _patient_description_results(sample: Sample) -> list:
    """ The phenotype text and its ontology matches - the wizard lists the proband's on submit """
    description = ''
    results = []
    try:
        description = sample.patient.phenotype
        results = sample.patient.phenotype_description.get_results()
    except (AttributeError, ObjectDoesNotExist):
        pass
    return [description, results]


class FamilyWizardView(FormView):
    """ Give each of a cohort's samples a role, then make the family group out of the roles.
        Subclasses set the family_* attributes and implement get_or_create_family() """
    template_name = "analysis/family_wizard.html"
    family_name: str = None  # Heading, button and title - "Trio"
    family_icon: str = None  # Sprite symbol the family's node and preview icon wear too
    figure_width: str = None  # The pedigree figure is drawn from a badge-sized symbol - see the template
    heading_icon_width: str = None
    role_help: str = None  # What to do with the roles, under the heading
    role_classes: dict = {}  # Role -> the class that fills that member of the pedigree figure in
    role_shape_classes: dict = {}  # Role -> the class that picks which shape it's drawn as (Duo)
    auto_assign_roles: dict = {}  # Sex -> the role a sample of that sex takes once the proband is picked
    auto_affected_roles: list = []  # Roles that tick their own affected box when picked

    def dispatch(self, request, *args, **kwargs):
        self.cohort = Cohort.get_for_user(request.user, kwargs["cohort_id"])
        if self.cohort.import_status != ImportStatus.SUCCESS:
            import_status = self.cohort.get_import_status_display()
            msg = f"Can't create analysis for {self.cohort} of status {import_status}"
            raise PermissionDenied(msg)

        sample_ids = [kwargs[f"sample{i}_id"] for i in range(1, len(self.form_class.SAMPLE_FIELDS) + 1)]
        self.samples = [Sample.get_for_user(request.user, sample_id) for sample_id in sample_ids]
        self.sample_sexes = _sample_sexes(self.samples)
        self.name_suggestions = family_name_suggestions(self.samples)
        return super().dispatch(request, *args, **kwargs)

    def get_initial(self):
        initial = super().get_initial()
        if self.name_suggestions:
            initial["name"] = self.name_suggestions[0]
        return initial

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["patient_sexes"] = [sample.patient_sex for sample in self.samples]
        return kwargs

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        form = context["form"]
        sample_rows = []
        for i, sample in enumerate(self.samples, start=1):
            withheld_roles = form.roles_withheld_for_sex(sample.patient_sex)
            sample_rows.append({"sample": sample,
                                "index": i,
                                "withheld_roles": " or ".join(role.label for role in withheld_roles),
                                "role_field": form[f"sample_{i}"],
                                "affected_field": form[f"sample_{i}_affected"]})
        context.update({
            "cohort": self.cohort,
            "sample_rows": sample_rows,
            "sample_sexes": self.sample_sexes,
            "name_suggestions": self.name_suggestions,
            "patient_description_results": [_patient_description_results(s) for s in self.samples],
            "family_name": self.family_name,
            "family_icon": self.family_icon,
            "figure_width": self.figure_width,
            "heading_icon_width": self.heading_icon_width,
            "role_help": self.role_help,
            "role_options": {"roleClasses": self.role_classes,
                             "roleShapeClasses": self.role_shape_classes,
                             "parentRoleSexes": {role: sex for sex, role in form.PARENT_ROLE_BY_SEX.items()},
                             "autoAssignRoles": self.auto_assign_roles,
                             "autoAffectedRoles": self.auto_affected_roles},
        })
        return context

    def form_valid(self, form):
        cohort_samples = {role: self.cohort.cohortsample_set.get(sample=sample)
                          for sample, role in zip(self.samples, form.roles)}
        # Roles are declared parents first, so the name reads down the pedigree
        names = "/".join(cohort_samples[role].name for role in form.ROLE_ENUM if role in cohort_samples)
        name = form.cleaned_data["name"].strip()
        proband_sex = form.cleaned_data["proband_sex"] or None
        defaults = {"name": name or f"{names} from {self.cohort}", "proband_sex": proband_sex}

        family, created = self.get_or_create_family(form, cohort_samples, defaults)
        if not created:
            # The same roles find the family made before - what was asked for this time wins
            family.proband_sex = proband_sex
            if name:
                family.name = name
            family.save()
        return redirect(family)

    def get_or_create_family(self, form, cohort_samples: dict, defaults: dict) -> tuple:
        raise NotImplementedError()


class TrioWizardView(FamilyWizardView):
    form_class = UserTrioWizardForm
    family_name = "Trio"
    family_icon = "node-icon-trio"
    figure_width = "150px"
    heading_icon_width = "1.5em"
    role_help = ("give each sample a role in the family, and tick the parents that are affected. "
                 "Picking the proband fills in the parents where their sexes settle it.")
    role_classes = {TrioSample.MOTHER: "mother-affected", TrioSample.FATHER: "father-affected"}
    auto_assign_roles = {Sex.FEMALE: TrioSample.MOTHER, Sex.MALE: TrioSample.FATHER}

    def get_or_create_family(self, form, cohort_samples, defaults):
        affected = form.affected_by_role
        return Trio.objects.get_or_create(cohort=self.cohort,
                                          user=self.request.user,
                                          mother=cohort_samples[TrioSample.MOTHER],
                                          mother_affected=affected[TrioSample.MOTHER],
                                          father=cohort_samples[TrioSample.FATHER],
                                          father_affected=affected[TrioSample.FATHER],
                                          proband=cohort_samples[TrioSample.PROBAND],
                                          defaults=defaults)


class QuadWizardView(FamilyWizardView):
    form_class = UserQuadWizardForm
    family_name = "Quad"
    family_icon = "node-icon-quad"
    figure_width = "168px"
    heading_icon_width = "1.7em"
    role_help = "give each sample a role in the family, and tick the members that are affected."
    role_classes = {QuadSample.MOTHER: "mother-affected", QuadSample.FATHER: "father-affected",
                    QuadSample.SIBLING: "sibling-affected"}

    def get_or_create_family(self, form, cohort_samples, defaults):
        affected = form.affected_by_role
        return Quad.objects.get_or_create(cohort=self.cohort,
                                          user=self.request.user,
                                          mother=cohort_samples[QuadSample.MOTHER],
                                          mother_affected=affected[QuadSample.MOTHER],
                                          father=cohort_samples[QuadSample.FATHER],
                                          father_affected=affected[QuadSample.FATHER],
                                          proband=cohort_samples[QuadSample.PROBAND],
                                          sibling=cohort_samples[QuadSample.SIBLING],
                                          sibling_affected=affected[QuadSample.SIBLING],
                                          defaults=defaults)


class DuoWizardView(FamilyWizardView):
    form_class = UserDuoWizardForm
    family_name = "Duo"
    family_icon = "node-icon-duo"
    figure_width = "90px"
    heading_icon_width = "1.1em"
    role_help = ("say which relative you have - a parent, or a sibling - and which sample is the "
                 "proband, then tick the relative if they're affected. Picking the proband fills a "
                 "parent in where its sex settles it.")
    role_classes = {DuoSample.MOTHER: "relative-affected", DuoSample.FATHER: "relative-affected",
                    DuoSample.SIBLING: "relative-affected"}
    role_shape_classes = {DuoSample.MOTHER: "duo-mother", DuoSample.FATHER: "duo-father",
                          DuoSample.SIBLING: "duo-sibling"}
    auto_assign_roles = {Sex.FEMALE: DuoSample.MOTHER, Sex.MALE: DuoSample.FATHER}
    # A sibling duo is usually two affected siblings - #1861
    auto_affected_roles = [DuoSample.SIBLING]

    def get_or_create_family(self, form, cohort_samples, defaults):
        relative_role = form.relative_role
        return Duo.objects.get_or_create(cohort=self.cohort,
                                         user=self.request.user,
                                         proband=cohort_samples[DuoSample.PROBAND],
                                         relative=cohort_samples[relative_role],
                                         relationship=relative_role,
                                         relative_affected=form.affected_by_role[relative_role],
                                         defaults=defaults)
