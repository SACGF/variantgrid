from typing import List
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import QuerySet, Q
from django.http import HttpRequest
from django.shortcuts import render
from email_manager.models import EmailLog
from library.django_utils import require_superuser
from library.utils.django_utils import render_ajax_view
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, CellData


@require_superuser
def email_manager_view(request):
    return render(request, 'email_log.html', context={})


@require_superuser
def email_detail(request, email_id: int):
    email_log = EmailLog.objects.get(pk=email_id)
    recipients = [address.strip() for address in email_log.recipient_list.split(";")]
    recipients = [address for address in recipients if address]
    users: List[User] = list(User.objects.filter(email__in=recipients).order_by('email').all())
    unrecognised_emails = set(recipients)
    check_user: User
    for check_user in users:
        unrecognised_emails.discard(check_user.email)
    unrecognised_email_list = list(sorted(unrecognised_emails))

    return render_ajax_view(
        request,
        'email_detail.html',
        context={"email": email_log, "users": users, "unrecognised": unrecognised_email_list},
        menubar='settings'
    )


class EmailColumns(DatatableConfig[EmailLog]):

    def power_search(self, qs: QuerySet[EmailLog], search_string: str) -> QuerySet[EmailLog]:
        if search_string:
            qs = qs.filter(Q(recipient_list__icontains=search_string) | Q(subject__icontains=search_string))
        return qs

    def recipient_renderer(self, row: CellData):
        if filename := row.get('filename'):
            return filename
        elif detail := row.get('details'):
            return detail.split('\n', 1)[0]

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        if not self.user.is_superuser:
            raise PermissionDenied("Email log requires admin privileges")

        self.search_box_enabled = True

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('email_detail', expected_height=300)
        self.rich_columns = [
            RichColumn('id', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True, default_sort=SortOrder.DESC),
            RichColumn('subject', orderable=True),
            RichColumn('recipient_list', label='Recipients', orderable=False),
            RichColumn('probably_sent', client_renderer='TableFormat.boolean.bind(null, "false_is_error")'),
        ]

    def get_initial_queryset(self) -> QuerySet[EmailLog]:
        return EmailLog.objects.all()
