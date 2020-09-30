from rest_framework.urlpatterns import format_suffix_patterns

from flags.views.views import FlagsView, FlagView
from variantgrid.perm_path import perm_path

rest_urlpatterns = [
    perm_path('api/flags/<flag_collection_id>', FlagsView.as_view(), name='flags_api'),
    perm_path('api/flag/<int:flag_id>', FlagView.as_view(), name='flag_api')
]

urlpatterns = format_suffix_patterns(rest_urlpatterns)
