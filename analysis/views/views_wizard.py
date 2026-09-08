""" The Duo/Trio/Quad wizards - turn a cohort's samples into a family group.

    Every wizard asks the same questions (which sample has which role, who's affected, and which sex
    to trust for the proband) so they share FamilyWizardView and analysis/family_wizard.html - a
    subclass names the family, its form and how to make the object from the roles.
"""
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.shortcuts import redirect
from django.views.generic import FormView

from analysis.forms import UserDuoWizardForm, UserQuadWizardForm, UserTrioWizardForm
from analysis.models.enums import DuoSample, QuadSample, TrioSample
from patients.models_enums import Sex
from snpdb.models import Cohort, Duo, ImportStatus, Quad, Sample, Trio


def _confident_sex(sample: Sample) -> Sex:
    """ The sex we'll hold the wizard's roles to: the patient record, or the chrX call when there is no
        record. Disagreement means we don't know - the proband sex mismatch warning asks about that """
    patient_sex = sample.patient_sex
    detected_sex = sample.detected_sex
    if patient_sex == Sex.UNKNOWN:
        return detected_sex
    if detected_sex in (Sex.UNKNOWN, patient_sex):
        return patient_sex
    return Sex.UNKNOWN


def _sample_sexes(samples: list[Sample]) -> list[dict]:
    """ The proband and the roles on offer are picked client side, so hand the JS every sample's sexes """
    return [{"patient_sex": s.patient_sex.value, "patient_sex_label": s.patient_sex.label,
             "detected_sex": s.detected_sex.value, "detected_sex_label": s.detected_sex.label,
             "sex": _confident_sex(s).value}
            for s in samples]


def _patient_description_results(sample: Sample) -> list:
    """ The phenotype text and its ontology matches - the wizard lists the proband's on submit """
    description = ''
    results = []
    try:
        description = sample.patient.phenotype
        results = sample.patient.patient_text_phenotype.phenotype_description.get_results()
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

    def dispatch(self, request, *args, **kwargs):
        self.cohort = Cohort.get_for_user(request.user, kwargs["cohort_id"])
        if self.cohort.import_status != ImportStatus.SUCCESS:
            import_status = self.cohort.get_import_status_display()
            msg = f"Can't create analysis for {self.cohort} of status {import_status}"
            raise PermissionDenied(msg)

        sample_ids = [kwargs[f"sample{i}_id"] for i in range(1, len(self.form_class.SAMPLE_FIELDS) + 1)]
        self.samples = [Sample.get_for_user(request.user, sample_id) for sample_id in sample_ids]
        self.sample_sexes = _sample_sexes(self.samples)
        return super().dispatch(request, *args, **kwargs)

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["sample_sexes"] = [Sex(ss["sex"]) for ss in self.sample_sexes]
        return kwargs

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        form = context["form"]
        sample_rows = []
        for i, sample in enumerate(self.samples, start=1):
            sample_rows.append({"sample": sample,
                                "index": i,
                                "role_field": form[f"sample_{i}"],
                                "affected_field": form[f"sample_{i}_affected"]})
        context.update({
            "cohort": self.cohort,
            "sample_rows": sample_rows,
            "sample_sexes": self.sample_sexes,
            "patient_description_results": [_patient_description_results(s) for s in self.samples],
            "family_name": self.family_name,
            "family_icon": self.family_icon,
            "figure_width": self.figure_width,
            "heading_icon_width": self.heading_icon_width,
            "role_help": self.role_help,
            "role_options": {"roleClasses": self.role_classes,
                             "roleShapeClasses": self.role_shape_classes,
                             "autoAssignRoles": self.auto_assign_roles},
        })
        return context

    def form_valid(self, form):
        cohort_samples = {role: self.cohort.cohortsample_set.get(sample=sample)
                          for sample, role in zip(self.samples, form.roles)}
        # Roles are declared parents first, so the name reads down the pedigree
        names = "/".join(cohort_samples[role].name for role in form.ROLE_ENUM if role in cohort_samples)
        proband_sex = form.cleaned_data["proband_sex"] or None
        defaults = {"name": f"{names} from {self.cohort}", "proband_sex": proband_sex}

        family, created = self.get_or_create_family(form, cohort_samples, defaults)
        if not created and family.proband_sex != proband_sex:
            family.proband_sex = proband_sex
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
    role_help = ("say which parent you have and which sample is the proband, then tick the parent if "
                 "they're affected. Picking the proband fills the parent in where its sex settles it.")
    role_classes = {DuoSample.MOTHER: "parent-affected", DuoSample.FATHER: "parent-affected"}
    role_shape_classes = {DuoSample.MOTHER: "duo-mother", DuoSample.FATHER: "duo-father"}
    auto_assign_roles = {Sex.FEMALE: DuoSample.MOTHER, Sex.MALE: DuoSample.FATHER}

    def get_or_create_family(self, form, cohort_samples, defaults):
        parent_role = form.parent_role
        return Duo.objects.get_or_create(cohort=self.cohort,
                                         user=self.request.user,
                                         proband=cohort_samples[DuoSample.PROBAND],
                                         parent=cohort_samples[parent_role],
                                         relationship=parent_role,
                                         parent_affected=form.affected_by_role[parent_role],
                                         defaults=defaults)
