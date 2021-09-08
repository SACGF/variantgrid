from django.contrib import admin

from classification.models.upload_file_lab import UploadedFileLab
from snpdb.admin_utils import ModelAdminBasics


@admin.register(UploadedFileLab)
class ConditionTextAdmin(ModelAdminBasics):
    pass