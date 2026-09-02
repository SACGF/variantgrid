from django.core.exceptions import PermissionDenied
from django.shortcuts import redirect, render

from analysis.forms import UserQuadWizardForm, UserTrioWizardForm
from analysis.models.enums import QuadSample, TrioSample
from snpdb.models import Cohort, ImportStatus, Quad, Sample, Trio


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

    # The proband is picked client side, so hand the JS every sample's sexes to compare
    sample_sexes = [{"patient_sex": s.patient_sex.value, "patient_sex_label": s.patient_sex.label,
                     "detected_sex": s.detected_sex.value, "detected_sex_label": s.detected_sex.label}
                    for s in samples]

    form = UserTrioWizardForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            mother_affected = form.cleaned_data['mother_affected']
            father_affected = form.cleaned_data['father_affected']
            proband_sex = form.cleaned_data['proband_sex'] or None
            sample_1_person = form.cleaned_data['sample_1']
            sample_2_person = form.cleaned_data['sample_2']
            sample_3_person = form.cleaned_data['sample_3']

            people = [sample_1_person, sample_2_person, sample_3_person]

            mother_cs = None
            father_cs = None
            proband_cs = None

            def get_cohort_sample(sample):
                return cohort.cohortsample_set.get(sample=sample)

            for s, p in zip(samples, people):
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

    # The proband is picked client side, so hand the JS every sample's sexes to compare
    sample_sexes = [{"patient_sex": s.patient_sex.value, "patient_sex_label": s.patient_sex.label,
                     "detected_sex": s.detected_sex.value, "detected_sex_label": s.detected_sex.label}
                    for s in samples]

    form = UserQuadWizardForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            mother_affected  = form.cleaned_data['mother_affected']
            father_affected  = form.cleaned_data['father_affected']
            sibling_affected = form.cleaned_data['sibling_affected']
            proband_sex = form.cleaned_data['proband_sex'] or None
            sample_roles = [form.cleaned_data[f'sample_{i}'] for i in range(1, 5)]

            def get_cohort_sample(sample):
                return cohort.cohortsample_set.get(sample=sample)

            mother_cs = father_cs = proband_cs = sibling_cs = None
            for sample, role in zip(samples, sample_roles):
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
