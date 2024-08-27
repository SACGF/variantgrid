from review import views
from variantgrid.perm_path import path

urlpatterns = [
    path('new/<int:reviewed_object_id>/<str:topic_id>/', views.new_review, name='start_review'),
    path('detail/<int:review_id>', views.view_discussion_detail, name='review_detail'),
    path('<int:review_id>/', views.edit_review, name='edit_review')
]
