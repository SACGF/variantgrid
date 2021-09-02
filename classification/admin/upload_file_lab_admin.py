from classification.models.upload_file_lab import UploadedFileLab
from snpdb.admin_utils import ModelAdminBasics
from django.contrib import admin


@admin.register(UploadedFileLab)
class ConditionTextAdmin(ModelAdminBasics):
    pass