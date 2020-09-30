from django.contrib import admin
from manual import models

class ManualMigrationTaskAdmin(admin.ModelAdmin):
    list_display = ('id', 'created', 'task', 'note')

class ManualMigrationRequiredAdmin(admin.ModelAdmin):
    list_display = ('id', 'created', 'task', 'note')

admin.site.register(models.ManualMigrationTask)
admin.site.register(models.ManualMigrationAttempt, ManualMigrationTaskAdmin)
admin.site.register(models.ManualMigrationRequired, ManualMigrationRequiredAdmin)
