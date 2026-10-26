"""
@see https://django-autocomplete-light.readthedocs.io/en/master/
"""

from django.db.models.query_utils import Q
from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from patients.enums import SampleSourceLevel
from patients.models import ExternalPK, Extraction, Patient, Specimen
from snpdb.views.views_autocomplete import GenomeBuildAutocompleteView, SampleAutocompleteView


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class PatientAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['patient_code', 'last_name', 'first_name']

    def get_user_queryset(self, user):
        qs = Patient.filter_for_user(user)
        # A patient's samples reach it either way round - the VCF import carries extraction down
        # without setting sample.patient, the patient CSV sets patient and may leave extraction null
        return self.filter_to_readable_samples(qs, ["sample", "specimen__extraction__sample"])


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class SpecimenAutocompleteView(GenomeBuildAutocompleteView):
    """ Narrows on whichever of Patient -> Specimen -> Extraction a form has already set """
    fields = ['reference_id']

    def get_user_queryset(self, user):
        qs = Specimen.filter_for_user(user)
        if patient := self.forwarded.get('patient'):
            qs = qs.filter(patient=patient)
        if extraction := self.forwarded.get('extraction'):
            qs = qs.filter(extraction=extraction)
        # An analysis is one genome build, so only offer specimens it can actually read
        return self.filter_to_readable_samples(qs, ["extraction__sample"])


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class ExtractionAutocompleteView(GenomeBuildAutocompleteView):
    """ Narrows on whichever of Patient -> Specimen a form has already set """
    # An extraction may be referred to by its own reference (eg a container suffix) or its specimen's
    fields = ['reference_id', 'specimen__reference_id']

    def get_user_queryset(self, user):
        qs = Extraction.filter_for_user(user)
        if patient := self.forwarded.get('patient'):
            qs = qs.filter(specimen__patient=patient)
        if specimen := self.forwarded.get('specimen'):
            qs = qs.filter(specimen=specimen)
        # An analysis is one genome build, so only offer extractions it can actually read
        return self.filter_to_readable_samples(qs, ["sample"])


@method_decorator([cache_page(30), vary_on_cookie], name='dispatch')
class ExternalPKAutocompleteView(AutocompleteView):
    """ A code can be an MRN or accession, so only those on a record the user can view """
    fields = ['code']

    def get_user_queryset(self, user):
        visible_q = Q()
        for rel in ExternalPK._meta.related_objects:
            if rel.one_to_one:  # each ExternallyManagedModel
                visible_q |= Q(pk__in=rel.related_model.filter_for_user(user).values(rel.field.name))

        qs = ExternalPK.objects.filter(visible_q)
        if external_type := self.forwarded.get('external_type'):
            qs = qs.filter(external_type=external_type)
        return qs


# Result text carries the parent, so two "DNA" extractions are told apart in one flat list
SAMPLE_SOURCE_TEXT = {
    SampleSourceLevel.SPECIMEN: lambda obj: f"{obj} \u2014 {obj.patient}",
    SampleSourceLevel.EXTRACTION: lambda obj: f"{obj} \u2014 {obj.specimen}",
    SampleSourceLevel.SAMPLE: lambda obj: f"{obj.name} ({obj.vcf})",
}


# str() on a result walks up the hierarchy, so pull the parent in with the row
SAMPLE_SOURCE_SELECT_RELATED = {
    SampleSourceLevel.SPECIMEN: ("patient",),
    SampleSourceLevel.EXTRACTION: ("specimen", "specimen__patient"),
    SampleSourceLevel.SAMPLE: ("vcf",),
}


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class SampleSourceAutocompleteView(AutocompleteView):
    """ One search box over all four levels of the hierarchy, so the SampleNode editor asks for the
        thing rather than for a level and then a thing.

        select2 renders a result carrying a `children` array as a group natively, so each group is
        just the matching model's own autocomplete queryset - permissions, genome build and
        exclude_archived forwards keep applying. Ids are "<level>:<pk>". """
    GROUPS = [
        (SampleSourceLevel.PATIENT, "Patients", PatientAutocompleteView),
        (SampleSourceLevel.SPECIMEN, "Specimens", SpecimenAutocompleteView),
        (SampleSourceLevel.EXTRACTION, "Extractions", ExtractionAutocompleteView),
        (SampleSourceLevel.SAMPLE, "Samples", SampleAutocompleteView),
    ]
    MAX_PER_GROUP = 8

    def get_queryset(self):
        return Patient.objects.none()  # Each group brings its own - see get()

    def get(self, request, *args, **kwargs):
        results = []
        for level, label, view_class in self.GROUPS:
            view = view_class()
            view.request = request
            view.q = self.q
            view.forwarded = self.forwarded
            get_text = SAMPLE_SOURCE_TEXT.get(level, str)
            qs = view.get_queryset().select_related(*SAMPLE_SOURCE_SELECT_RELATED.get(level, ()))
            children = [{"id": f"{level}:{obj.pk}", "text": get_text(obj)}
                        for obj in qs[:self.MAX_PER_GROUP]]
            if children:
                results.append({"text": label, "children": children})
        return JsonResponse({"results": results, "pagination": {"more": False}})
