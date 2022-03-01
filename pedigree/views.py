from crispy_forms.helper import FormHelper
from django.conf import settings
from django.forms.forms import Form
from django.forms.formsets import formset_factory
from django.forms.models import ModelChoiceField
from django.http.response import HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse

from library.utils import full_class_name
from pedigree import forms
from pedigree.forms import BaseCohortSamplesForPedFileRecordsFormSet
from pedigree.graphs.pedigree_chart import PedigreeChart
from pedigree.models import PedFile, Pedigree, PedFileFamily, PedFileRecord, \
    CohortSamplePedFileRecord, create_automatch_pedigree
from snpdb.forms import UserCohortForm
from snpdb.graphs import graphcache
from snpdb.models import Cohort


def view_ped_file(request, ped_file_id):
    ped_file = PedFile.get_for_user(request.user, ped_file_id)
    ped_file_form = forms.ReadOnlyPedFileForm(instance=ped_file)
    draw_pedigree_chart = bool(getattr(settings, "PEDIGREE_MADELINE2_COMMAND", None))

    context = {'ped_file': ped_file,
               'ped_file_form': ped_file_form,
               'cohort_form': UserCohortForm(),
               'draw_pedigree_chart': draw_pedigree_chart}
    return render(request, 'pedigree/view_ped_file.html', context)


def pedigree_chart(request, ped_file_id):
    ped_file = PedFile.get_for_user(request.user, ped_file_id)  # Make sure we can access it
    graph_class_name = full_class_name(PedigreeChart)
    cached_graph = graphcache.async_graph(graph_class_name, ped_file.pk)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def ped_files(request):
    return render(request, 'pedigree/ped_files.html')


def pedigrees(request):
    form = forms.UserCohortandPedFileFamilyForm(user=request.user)
    context = {'form': form}

    return render(request, 'pedigree/pedigrees.html', context)


def view_pedigree(request, pedigree_id):
    pedigree = Pedigree.get_for_user(request.user, pedigree_id)

    cohort_samples_qs = pedigree.cohort.cohortsample_set.all()

    class CohortSamplesForm(Form):
        ped_file_record = ModelChoiceField(disabled=True, queryset=PedFileRecord.objects.all())
        cohort_sample = ModelChoiceField(required=False, queryset=cohort_samples_qs)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_show_labels = False

    CohortSamplesForPedFileRecordsFormSet = formset_factory(CohortSamplesForm, formset=BaseCohortSamplesForPedFileRecordsFormSet, extra=0)
    # TODO: Make this non read only and save name etc.
    pedigree_form = forms.PedigreeForm(request.POST or None, instance=pedigree)
    formset_initial = []
    existing_cspfr = pedigree.cohortsamplepedfilerecord_set.all()
    if existing_cspfr.exists():
        for cspfr in existing_cspfr:
            formset_initial.append({"ped_file_record": cspfr.ped_file_record.pk,
                                    "cohort_sample": cspfr.cohort_sample})
    else:
        pf_records_qs = pedigree.ped_file_family.pedfilerecord_set.all()
        for pfr in pf_records_qs:
            formset_initial.append({"ped_file_record": pfr.pk})

    formset = CohortSamplesForPedFileRecordsFormSet(request.POST or None, initial=formset_initial)
    if request.method == "POST":
        if pedigree_form.is_valid() and formset.is_valid():
            pedigree = pedigree_form.save()
            # Clear all existing records
            pedigree.cohortsamplepedfilerecord_set.all().delete()

            for form in formset:
                cohort_sample = form.cleaned_data['cohort_sample']
                if cohort_sample:
                    CohortSamplePedFileRecord.objects.create(pedigree=pedigree,
                                                             cohort_sample=cohort_sample,
                                                             ped_file_record=form.cleaned_data['ped_file_record'],)

    context = {'pedigree': pedigree,
               "cohort": pedigree.cohort,
               "ped_file_family": pedigree.ped_file_family,
               'pedigree_form': pedigree_form,
               "formset": formset,
               "has_write_permission": pedigree.can_write(request.user)}

    return render(request, 'pedigree/view_pedigree.html', context)


def create_pedigree_from_cohort_and_ped_file_family(request, cohort_id, ped_file_family_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    ped_file_family = get_object_or_404(PedFileFamily, pk=ped_file_family_id)

    pedigree = create_automatch_pedigree(request.user, ped_file_family, cohort)
    return HttpResponseRedirect(reverse('view_pedigree', kwargs={'pedigree_id': pedigree.pk}))
