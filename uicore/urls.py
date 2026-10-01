from django.urls import path

from uicore.views.chrome_view import page_chrome

# Django's path rather than variantgrid.perm_path.path: every page needs its chrome, whatever URLS_NAME_REGISTER hides
urlpatterns = [
    path('chrome', page_chrome, name='page_chrome'),
]
