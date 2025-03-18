from assay_detailed.models import AssayDetailedRNA
from snpdb.admin_utils import ModelAdminBasics
from django.contrib import admin, messages


@admin.register(AssayDetailedRNA)
class AssayAdmin(ModelAdminBasics):
    pass