from django.core.exceptions import PermissionDenied
from django.shortcuts import redirect, render

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


def trio_wizard(request, cohort_id, sample1_id, sample2_id, sample3_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    if cohort.import_status != ImportStatus.SUCCESS:
        import_status = cohort.get_import_status_display()
        msg = f"Can't create analysis for {cohort} of status {import_status}"
        raise PermissionDenied(msg)

    sample_1 = Sample.get_for_user(request.user, sample1_id)
    sample_2 = Sample.get_for_user(request.user, sample2_id)
    sample_3 = Sample.get_for_user(request.user, sample3_id)

    samples = [sample_1, sample_2, sample_3]
    patient_description_results = []
    for sample in samples:
        description = ''
        results = []
        try:
            description = sample.patient.phenotype
            results = sample.patient.patient_text_phenotype.phenotype_description.get_results()
        except:
            pass
        patient_description_results.append([description, results])

    sample_sexes = _sample_sexes(samples)

    form = UserTrioWizardForm(request.POST or None,
                              sample_sexes=[Sex(ss["sex"]) for ss in sample_sexes])
    if request.method == "POST":
        if form.is_valid():
            affected_by_role = form.affected_by_role
            mother_affected = affected_by_role[TrioSample.MOTHER]
            father_affected = affected_by_role[TrioSample.FATHER]
            proband_sex = form.cleaned_data['proband_sex'] or None

            mother_cs = None
            father_cs = None
            proband_cs = None

            def get_cohort_sample(sample):
                return cohort.cohortsample_set.get(sample=sample)

            for s, p in zip(samples, form.roles):
                if p == TrioSample.FATHER:
                    father_cs = get_cohort_sample(s)
                elif p == TrioSample.MOTHER:
                    mother_cs = get_cohort_sample(s)
                elif p == TrioSample.PROBAND:
                    proband_cs = get_cohort_sample(s)

            sample_names = "/".join((mother_cs.name, father_cs.name, proband_cs.name))
            trio_name = f"{sample_names} from {cohort}"
            trio, created = Trio.objects.get_or_create(cohort=cohort,
                                                       user=request.user,
                                                       mother=mother_cs,
                                                       mother_affected=mother_affected,
                                                       father=father_cs,
                                                       father_affected=father_affected,
                                                       proband=proband_cs,
                                                       defaults={"name": trio_name, "proband_sex": proband_sex})
            if not created and trio.proband_sex != proband_sex:
                trio.proband_sex = proband_sex
                trio.save()
            return redirect(trio)

    context = {"cohort": cohort,
               "sample_1": sample_1,
               "sample_2": sample_2,
               "sample_3": sample_3,
               "form": form,
               "sample_sexes": sample_sexes,
               "patient_description_results": patient_description_results}
    return render(request, 'analysis/trio_wizard.html', context)


def quad_wizard(request, cohort_id, sample1_id, sample2_id, sample3_id, sample4_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    if cohort.import_status != ImportStatus.SUCCESS:
        import_status = cohort.get_import_status_display()
        msg = f"Can't create analysis for {cohort} of status {import_status}"
        raise PermissionDenied(msg)

    samples = [Sample.get_for_user(request.user, sid)
               for sid in [sample1_id, sample2_id, sample3_id, sample4_id]]

    patient_description_results = []
    for sample in samples:
        description = ''
        results = []
        try:
            description = sample.patient.phenotype
            results = sample.patient.patient_text_phenotype.phenotype_description.get_results()
        except:
            pass
        patient_description_results.append([description, results])

    sample_sexes = _sample_sexes(samples)

    form = UserQuadWizardForm(request.POST or None,
                              sample_sexes=[Sex(ss["sex"]) for ss in sample_sexes])
    if request.method == "POST":
        if form.is_valid():
            affected_by_role = form.affected_by_role
            mother_affected = affected_by_role[QuadSample.MOTHER]
            father_affected = affected_by_role[QuadSample.FATHER]
            sibling_affected = affected_by_role[QuadSample.SIBLING]
            proband_sex = form.cleaned_data['proband_sex'] or None

            def get_cohort_sample(sample):
                return cohort.cohortsample_set.get(sample=sample)

            mother_cs = father_cs = proband_cs = sibling_cs = None
            for sample, role in zip(samples, form.roles):
                cs = get_cohort_sample(sample)
                if role == QuadSample.MOTHER:
                    mother_cs = cs
                elif role == QuadSample.FATHER:
                    father_cs = cs
                elif role == QuadSample.PROBAND:
                    proband_cs = cs
                elif role == QuadSample.SIBLING:
                    sibling_cs = cs

            quad_name = "/".join((mother_cs.name, father_cs.name, proband_cs.name, sibling_cs.name))
            quad_name += f" from {cohort}"
            quad, created = Quad.objects.get_or_create(
                cohort=cohort,
                user=request.user,
                mother=mother_cs,
                mother_affected=mother_affected,
                father=father_cs,
                father_affected=father_affected,
                proband=proband_cs,
                sibling=sibling_cs,
                sibling_affected=sibling_affected,
                defaults={"name": quad_name, "proband_sex": proband_sex},
            )
            if not created and quad.proband_sex != proband_sex:
                quad.proband_sex = proband_sex
                quad.save()
            return redirect(quad)

    context = {
        "cohort": cohort,
        "sample_1": samples[0],
        "sample_2": samples[1],
        "sample_3": samples[2],
        "sample_4": samples[3],
        "form": form,
        "sample_sexes": sample_sexes,
        "patient_description_results": patient_description_results,
    }
    return render(request, 'analysis/quad_wizard.html', context)


def duo_wizard(request, cohort_id, sample1_id, sample2_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    if cohort.import_status != ImportStatus.SUCCESS:
        import_status = cohort.get_import_status_display()
        msg = f"Can't create analysis for {cohort} of status {import_status}"
        raise PermissionDenied(msg)

    samples = [Sample.get_for_user(request.user, sid) for sid in [sample1_id, sample2_id]]

    patient_description_results = []
    for sample in samples:
        description = ''
        results = []
        try:
            description = sample.patient.phenotype
            results = sample.patient.patient_text_phenotype.phenotype_description.get_results()
        except:
            pass
        patient_description_results.append([description, results])

    sample_sexes = _sample_sexes(samples)

    form = UserDuoWizardForm(request.POST or None,
                             sample_sexes=[Sex(ss["sex"]) for ss in sample_sexes])
    if request.method == "POST":
        if form.is_valid():
            parent_affected = form.affected_by_role[DuoSample.PARENT]
            relationship = form.parent_relationship
            proband_sex = form.cleaned_data['proband_sex'] or None

            def get_cohort_sample(sample):
                return cohort.cohortsample_set.get(sample=sample)

            parent_cs = proband_cs = None
            for sample, role in zip(samples, form.roles):
                cs = get_cohort_sample(sample)
                if role == DuoSample.PARENT:
                    parent_cs = cs
                elif role == DuoSample.PROBAND:
                    proband_cs = cs

            duo_name = f"{parent_cs.name}/{proband_cs.name} from {cohort}"
            duo, created = Duo.objects.get_or_create(
                cohort=cohort,
                user=request.user,
                proband=proband_cs,
                parent=parent_cs,
                relationship=relationship,
                parent_affected=parent_affected,
                defaults={"name": duo_name, "proband_sex": proband_sex},
            )
            if not created and duo.proband_sex != proband_sex:
                duo.proband_sex = proband_sex
                duo.save()
            return redirect(duo)

    context = {
        "cohort": cohort,
        "sample_1": samples[0],
        "sample_2": samples[1],
        "form": form,
        "sample_sexes": sample_sexes,
        "patient_description_results": patient_description_results,
    }
    return render(request, 'analysis/duo_wizard.html', context)
