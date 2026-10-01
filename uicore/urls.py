from django.urls import path

from uicore.views.page_frame_view import page_frame

# Django's path rather than variantgrid.perm_path.path: every page needs its page frame, whatever URLS_NAME_REGISTER hides
urlpatterns = [
    path('page_frame', page_frame, name='page_frame'),
]
