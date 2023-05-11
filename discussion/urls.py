from discussion import views
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('discuss/<int:discussed_object_pk>/<str:topic_pk>/', views.start_discussion, name='start_discussion')
]