from django.contrib import admin

from library.django_utils.admin_utils import ModelAdminBasics
from manual.models import ManualMigrationAttempt, ManualMigrationRequired


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
