from django.contrib import admin
from manual import models
from manual.models import ManualMigrationAttempt, ManualMigrationRequired
from snpdb.admin_utils import ModelAdminBasics


@admin.register(ManualMigrationAttempt)
class ManualMigrationAttemptAdmin(ModelAdminBasics):
    list_display = ('id', 'created', 'task', 'note')

    def has_add_permission(self, request):
        return False


@admin.register(ManualMigrationRequired)
class ManualMigrationRequiredAdmin(ModelAdminBasics):
    list_display = ('id', 'created', 'task', 'note')

    def has_add_permission(self, request):
        return False
