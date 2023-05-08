from discussion import views
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('discussion/<int:answer_group>', views.view_discussion, name='discussion_answer')
]