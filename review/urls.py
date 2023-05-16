from review import views
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('new/<int:reviewed_object_pk>/<str:topic_pk>/', views.new_review, name='start_review'),
    perm_path('<int:answer_group_pk>/', views.edit_review, name='edit_review')
]