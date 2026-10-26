"""
The "Classification" tab on the lab page (snpdb/templates/snpdb/settings/view_lab.html): the lab's
LabClassificationSettings as a form, saved by lab writers. Entry point: lab_classification_settings_tab.
"""
from django.core.exceptions import PermissionDenied
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

from classification.forms import LabClassificationSettingsForm
from library.django_utils import add_save_message, set_form_read_only
from snpdb.models.models import Lab


def lab_classification_settings_tab(request, pk: int):
    lab = get_object_or_404(Lab, pk=pk)
    if not (lab.is_member(request.user) or request.user.is_superuser):
        raise PermissionDenied(f"You are not a member of {lab}")
    form = LabClassificationSettingsForm(request.POST or None, instance=lab.classification_settings)
    has_write_permission = lab.can_write(request.user)

    if request.method == "POST":
        lab.check_can_write(request.user)
        valid = form.is_valid()
        if valid:
            form.save()
        add_save_message(request, valid, "Lab Classification Settings")
        return redirect(reverse('view_lab', kwargs={"lab_id": lab.pk}))

    if not has_write_permission:
        set_form_read_only(form)

    context = {
        "lab": lab,
        "form": form,
        "has_write_permission": has_write_permission,
    }
    return render(request, 'classification/lab_classification_settings_tab.html', context)
