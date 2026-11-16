"""
Lab request pages for snpdb/models/models_lab_request.py:LabRequest: a user with a private lab asks for a real one
(lab_request) and confirms their email from the link it sends (lab_request_confirm); admins review requests
(lab_requests, view_lab_request), approving through snpdb/lab_request_approval.py:approve_lab_request.
"""
from django.contrib import messages
from django.core import signing
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_POST

from library.django_utils import require_superuser
from snpdb.enums import LabRequestStatus
from snpdb.forms import LabRequestApproveForm, LabRequestForm
from snpdb.lab_request_approval import approve_lab_request
from snpdb.models.models_lab_request import LabRequest, get_user_private_lab


def lab_request(request):
    user = request.user
    private_lab = get_user_private_lab(user)
    pending = LabRequest.objects.filter(user=user, status=LabRequestStatus.PENDING).first()

    form = None
    if private_lab and not pending:
        form = LabRequestForm(request.POST or None, initial={"email": user.email})
        if form.is_valid():
            new_request = form.save(commit=False)
            new_request.user = user
            new_request.save()
            new_request.send_confirmation_email()
            messages.info(request, f"We've emailed a confirmation link to {new_request.email}")
            return redirect("lab_request")

    context = {
        "private_lab": private_lab,
        "pending": pending,
        "form": form,
        "lab_requests": LabRequest.objects.filter(user=user).exclude(status=LabRequestStatus.PENDING),
    }
    return render(request, "snpdb/lab_request/lab_request.html", context)


@require_POST
def lab_request_cancel(request):
    LabRequest.objects.filter(user=request.user, status=LabRequestStatus.PENDING).delete()
    return redirect("lab_request")


def lab_request_confirm(request, token):
    try:
        confirmed_request = LabRequest.from_confirmation_token(token)
    except signing.BadSignature:
        messages.error(request, "That confirmation link is invalid or has expired - cancel the request and make it again")
    else:
        if confirmed_request.status == LabRequestStatus.PENDING:
            confirmed_request.confirm_email()
            messages.success(request, "Email confirmed - an administrator will review your request")
    return redirect("lab_request")


@require_superuser
def lab_requests(request):
    qs = LabRequest.objects.select_related("user", "reviewed_by", "lab")
    context = {
        "awaiting_review": qs.filter(status=LabRequestStatus.PENDING, email_confirmed__isnull=False),
        "unconfirmed": qs.filter(status=LabRequestStatus.PENDING, email_confirmed__isnull=True),
        "reviewed": qs.exclude(status=LabRequestStatus.PENDING),
    }
    return render(request, "snpdb/lab_request/lab_requests.html", context)


@require_superuser
def view_lab_request(request, pk):
    reviewed_request = get_object_or_404(LabRequest, pk=pk)
    private_lab = get_user_private_lab(reviewed_request.user)
    approve_form = None
    if reviewed_request.can_review and private_lab:
        initial = {"organization_name": reviewed_request.organization_name, "lab_name": reviewed_request.lab_name}
        approve_form = LabRequestApproveForm(request.POST or None, initial=initial)
        if request.method == "POST":
            if request.POST.get("action") == "reject":
                reviewed_request.reject(request.user, notify=bool(request.POST.get("notify")))
                messages.info(request, "Lab request rejected")
                return redirect("lab_requests")
            if approve_form.is_valid():
                data = approve_form.cleaned_data
                lab = approve_lab_request(reviewed_request, request.user, data["organization"],
                                          data["organization_name"], data["lab_name"])
                messages.success(request, f"Created {lab}")
                return redirect(reverse("view_lab", kwargs={"lab_id": lab.pk}))

    context = {
        "lab_request": reviewed_request,
        "approve_form": approve_form,
        "private_lab": private_lab,
    }
    return render(request, "snpdb/lab_request/view_lab_request.html", context)
