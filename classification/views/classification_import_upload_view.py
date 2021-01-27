import os

from django.contrib.auth.models import User
from django.contrib import messages
from django.shortcuts import render
from django.views import View

from library.log_utils import report_event
from snpdb.models import Lab


class FileUploadView(View):

    def get(self, request, **kwargs):
        user: User = request.user
        labs = Lab.valid_labs_qs(user=user, admin_check=True).filter(upload_location__isnull=False).exclude(upload_location__iexact='')
        context = {"labs": labs}
        return render(request, 'classification/import_upload.html', context)

    def post(self, requests, **kwargs):
        # lazily have s3boto3 requirements
        from storages.backends.s3boto3 import S3Boto3Storage

        lab_str = requests.POST.get('lab')
        if not lab_str:
            raise ValueError("lab required")

        lab = Lab.objects.get(group_name=lab_str)
        if not lab.is_member(user=requests.user, admin_check=True):
            raise PermissionError("User does not have access to lab")

        protocol, path = lab.upload_location.split(":", maxsplit=1)

        if protocol == "s3":
            bucket, directory = path.split("/", maxsplit=1)

            file_obj = requests.FILES.get('file', '')

            # do your validation here e.g. file size/type check

            # organize a path for the file in bucket
            file_directory_within_bucket = directory

            # synthesize a full file path; note that we included the filename
            file_path_within_bucket = os.path.join(
                file_directory_within_bucket,
                file_obj.name
            )

            media_storage = S3Boto3Storage(bucket_name=bucket)

            media_storage.save(file_path_within_bucket, file_obj)
            file_url = media_storage.url(file_path_within_bucket)

            report_event(name="s3 upload", request=requests, extra_data={
                "filename": file_obj.name,
                "destination": file_url
            })
            messages.add_message(requests, messages.INFO, f"File {file_obj.name} uploaded to {lab.name}, importing should begin shortly. See the classifications sub-menu.")
            return self.get(requests)
        else:
            raise ValueError("Only s3 storage is currently supported for lab uploads")
