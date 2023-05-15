from discussion import views
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('discuss/<int:discussed_object_pk>/<str:topic_pk>/', views.new_discussion, name='start_discussion'),
    perm_path('discuss/<int:answer_group_pk>/', views.edit_discussion, name='edit_discussion')
]