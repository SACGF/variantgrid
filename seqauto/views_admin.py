from django.conf import settings
from django.shortcuts import render
from django.utils.decorators import method_decorator
from django.views.generic.edit import UpdateView

from library.django_utils import get_model_fields, staff_only
from seqauto.models import (
    Aligner,
    Assay,
    Library,
    Sequencer,
    VariantCaller,
    VariantCallingPipeline,
)


def get_sequencing_software_versions_template():
    if settings.SEQAUTO_ENABLED:
        base_template = "seqauto/menu_sequencing_data_base.html"
    else:
        base_template = "snpdb/menu/menu_settings_base.html"
    return base_template


def sequencing_software_versions(request):
    # TODO: Forms etc
    context = {"base_template": get_sequencing_software_versions_template()}
    return render(request, 'seqauto/sequencing_software_versions.html', context)


@method_decorator([staff_only], name='dispatch')
class SeqautoUpdateView(UpdateView):
    template_name = "snpdb/update_form.html"
    widgets = {}

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = get_sequencing_software_versions_template()
        context['title'] = str(self.model.__name__)
        return context

    def get_form(self, form_class=None):
        form = super().get_form(form_class)
        for f in form.fields:
            w = self.widgets.get(f)
            if w:
                form.fields[f].widget = w
        return form


class LibraryUpdate(SeqautoUpdateView):
    model = Library
    fields = get_model_fields(Library)
    #widgets = {"name" : TextInput()}


class SequencerUpdate(SeqautoUpdateView):
    model = Sequencer
    fields = get_model_fields(Sequencer)


class AssayUpdate(SeqautoUpdateView):
    model = Assay
    fields = get_model_fields(Assay)


class AlignerUpdate(SeqautoUpdateView):
    model = Aligner
    fields = get_model_fields(Aligner)


class VariantCallerUpdate(SeqautoUpdateView):
    model = VariantCaller
    fields = get_model_fields(VariantCaller)


class VariantCallingPipelineUpdate(SeqautoUpdateView):
    model = VariantCallingPipeline
    fields = get_model_fields(VariantCallingPipeline)
