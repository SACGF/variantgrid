from django.urls import re_path
from django.views.generic import RedirectView

from user_messages import views

urlpatterns = [
    re_path(r'^$', RedirectView.as_view(permanent=True, url='inbox/'), name='messages_redirect'),
    re_path(r'^inbox/$', views.inbox, name='messages_inbox'),
    re_path(r'^outbox/$', views.outbox, name='messages_outbox'),
    re_path(r'^compose/$', views.compose, name='messages_compose'),
    re_path(r'^compose/(?P<recipient>[\w.@+-]+)/$', views.compose, name='messages_compose_to'),
    re_path(r'^view/(?P<message_id>[\d]+)/$', views.view, name='messages_detail'),
    re_path(r'^delete/(?P<message_id>[\d]+)/$', views.delete, name='messages_delete'),
    re_path(r'^undelete/(?P<message_id>[\d]+)/$', views.undelete, name='messages_undelete'),
    re_path(r'^trash/$', views.trash, name='messages_trash'),
]
