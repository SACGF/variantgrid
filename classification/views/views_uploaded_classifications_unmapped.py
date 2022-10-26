import datetime
from dataclasses import dataclass
from typing import List, Iterable, Any, Dict, Optional, Set

import django
from django.contrib import messages
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import QuerySet
from django.http import HttpResponseRedirect, HttpRequest, StreamingHttpResponse
from django.shortcuts import render
from django.urls import reverse
from django.utils.timezone import now
from django.views import View

from classification.models import ClassificationImportRun, resolve_uploaded_url_to_handle, FileHandle
from classification.models.uploaded_classifications_unmapped import UploadedClassificationsUnmapped, \
    UploadedClassificationsUnmappedStatus
from classification.tasks.classification_import_map_and_insert_task import ClassificationImportMapInsertTask
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder, report_exc_info
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class UploadedClassificationsUnmappedColumns(DatatableConfig[UploadedClassificationsUnmapped]):

    def effective_date_set(self, row: Dict[str, Any]):
        return (row.get('effective_modified') or row.get('created')).timestamp()

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC, client_renderer="idRenderer"),
            RichColumn(key="filename", orderable=True, client_renderer='fileDownloaderRenderer'),
            RichColumn(key="effective_modified", label='Modified', orderable=True, client_renderer='TableFormat.timestamp', extra_columns=['created'], sort_keys=['created'], renderer=self.effective_date_set),
            RichColumn(key="file_size", label='Size', orderable=True, client_renderer='TableFormat.sizeBytes'),
            RichColumn(key="status", label="Status", orderable=True, renderer=lambda x: UploadedClassificationsUnmappedStatus(x["status"]).label),
            RichColumn(key="comment", client_renderer='TableFormat.text')
        ]

    def get_initial_queryset(self) -> QuerySet[UploadedClassificationsUnmapped]:
        if lab_id := self.get_query_param("lab_id"):
            return UploadedClassificationsUnmapped.objects.filter(lab=int(lab_id))
        else:
            raise ValueError("Must pass in lab_id")


def download_classification_unmapped_file(request: HttpRequest, uploaded_classification_unmapped_id: int):
    user = request.user
    record = UploadedClassificationsUnmapped.objects.filter(pk=uploaded_classification_unmapped_id).filter(lab__in=Lab.valid_labs_qs(user, admin_check=True)).first()
    if not record:
        raise PermissionDenied("You do not have access to this file")

    file_info = record.file_data
    response = StreamingHttpResponse(file_info.stream(), content_type='streaming_content')
    response['Content-Disposition'] = f'attachment; filename="{file_info.filename}"'
    return response


def view_uploaded_classification_unmapped(request: HttpRequest, uploaded_classification_unmapped_id: int):
    user = request.user
    record = UploadedClassificationsUnmapped.objects.filter(pk=uploaded_classification_unmapped_id).filter(
        lab__in=Lab.valid_labs_qs(user, admin_check=True)).first()
    if not record:
        raise PermissionDenied("You do not have access to this file")

    return render(request, 'classification/uploaded_classifications_unmapped.html', {
        "record": record,
    })


def view_uploaded_classification_unmapped_detail(request: HttpRequest, uploaded_classification_unmapped_id: int):
    user = request.user
    record = UploadedClassificationsUnmapped.objects.filter(pk=uploaded_classification_unmapped_id).filter(
        lab__in=Lab.valid_labs_qs(user, admin_check=True)).first()
    if not record:
        raise PermissionDenied("You do not have access to this file")
    in_progress = record.status not in {
        UploadedClassificationsUnmappedStatus.Manual,
        UploadedClassificationsUnmappedStatus.Error,
        UploadedClassificationsUnmappedStatus.Validated,
        UploadedClassificationsUnmappedStatus.Processed}

    import_runs: Iterable[ClassificationImportRun] = ClassificationImportRun.objects.filter(from_file=record).order_by('pk')

    http_response = render(request, 'classification/uploaded_classifications_unmapped_detail.html', {
        "record": record,
        "import_runs": import_runs,
        "now": now(),
        "in_progress": in_progress
    })

    if in_progress:
        # data is still being processed, we should continue to reload this page
        http_response['Auto-refresh'] = str(5000)  # auto-refresh in 5 seconds
    return http_response


@dataclass
class S3File:
    filename: str
    modified: Any
    size: Any

    def __str__(self):
        return f"{self.filename} {self.modified} {self.size}"


@dataclass(frozen=True)
class FileMeta:
    name: str
    modified: datetime
    size: int

    @staticmethod
    def from_unmapped(unmapped: UploadedClassificationsUnmapped) -> Optional['FileMeta']:
        if unmapped.effective_modified and unmapped.file_size:
            return FileMeta(
                name=unmapped.filename,
                modified=unmapped.effective_modified,
                size=unmapped.file_size
            )

    @staticmethod
    def from_file_handle(handle: FileHandle) -> 'FileMeta':
        return FileMeta(
            name=handle.filename,
            modified=handle.modified,
            size=handle.size
        )


class UploadedClassificationsUnmappedView(View):

    def get(self, request, **kwargs):
        lab_picker = LabPickerData.from_request(
            request=request,
            selection=kwargs.get('lab_id'),
            view_name='classification_upload_unmapped_lab',
            multi_select=False)

        if redirect_response := lab_picker.check_redirect():
            return redirect_response

        context = {"lab_picker": lab_picker}
        selected_lab = lab_picker.selected_lab
        if upload_location := selected_lab.upload_location:
            if server_address := resolve_uploaded_url_to_handle(selected_lab.upload_location):
                existing: Set[FileMeta] = set()
                for unmapped in UploadedClassificationsUnmapped.objects.filter(lab=selected_lab):
                    if meta := FileMeta.from_unmapped(unmapped):
                        existing.add(meta)

                server_files: List[FileHandle] = list()
                for file in server_address.list():
                    meta = FileMeta.from_file_handle(file)
                    if meta not in existing:
                        server_files.append(file)

                context["server_files"] = server_files

        return render(request, 'classification/classification_upload_file_unmapped.html', context)

    def post(self, requests, **kwargs):
        # lazily have s3boto3 requirements

        django.utils.encoding.force_text1 = django.utils.encoding.force_str
        django.utils.encoding.smart_text = django.utils.encoding.smart_str

        lab_picker = LabPickerData.from_request(request=requests, selection=kwargs.get('lab_id'))
        lab = lab_picker.selected_lab
        if not lab:
            raise ValueError("Single Lab Required")

        try:
            if server_address := resolve_uploaded_url_to_handle(lab.upload_location):
                sub_file: FileHandle

                if existing_file := requests.POST.get('import-existing'):
                    sub_file = server_address.sub_file(existing_file)
                    if not sub_file.exists:
                        raise ValueError(f"File {existing_file} does not appear to exist")

                else:
                    file_obj = requests.FILES.get('file')
                    if not file_obj:
                        raise ValueError("No File Provided")
                    sub_file = server_address.save(file_obj)

                status = UploadedClassificationsUnmappedStatus.Pending if lab.upload_automatic else UploadedClassificationsUnmappedStatus.Manual
                user: User = requests.user
                uploaded_file = UploadedClassificationsUnmapped.objects.create(
                    url=sub_file.clean_url,
                    filename=sub_file.filename,
                    user=requests.user,
                    lab=lab,
                    status=status,
                    effective_modified=sub_file.modified,
                    file_size=sub_file.size
                )

                if lab.upload_automatic:
                    task = ClassificationImportMapInsertTask.si(uploaded_file.pk)
                    task.apply_async()
                    return HttpResponseRedirect(reverse("classification_upload_unmapped_status", kwargs={
                        "uploaded_classification_unmapped_id": uploaded_file.pk}))
                else:
                    # if no automatic process of automated file,
                    admin_url = get_url_from_view_path(
                        f"/admin/classification/uploadedfilelab/{uploaded_file.pk}/change")
                    notifier = NotificationBuilder(
                        message="File Uploaded"
                    ).add_header(":file_folder: File Uploaded"). \
                        add_field("For Lab", lab.name). \
                        add_field("By User", user.username). \
                        add_field("Path", sub_file.clean_url). \
                        add_markdown(f"*URL:* <{admin_url}>")

                    # Soon going to replace this with automated import

                    notifier.add_markdown("This file will need to be handled manually!")
                    messages.add_message(requests, messages.INFO,
                                         f"File {sub_file.filename} uploaded for {lab.name}. This file will be uploaded after manual review.")

                    notifier.send()
            else:
                raise ValueError("Only s3 storage is currently supported for lab uploads")

        except ValueError as ve:
            report_exc_info()
            messages.add_message(requests, messages.ERROR, str(ve))

        if lab:
            return HttpResponseRedirect(reverse("classification_upload_unmapped_lab", kwargs={"lab_id": lab.pk}))
        else:
            return HttpResponseRedirect(reverse("classification_upload_unmapped"))
