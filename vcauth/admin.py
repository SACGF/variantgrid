from django.contrib import admin
from django.contrib.auth.models import User
from vcauth.user_admin import CustomUserAdmin

admin.site.unregister(User)
admin.site.register(User, CustomUserAdmin)
