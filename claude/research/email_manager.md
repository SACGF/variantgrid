# email_manager — research notes

Verified against 7c4408c62 on 2026-09-06

The app exists so that every email VariantGrid sends leaves a row behind. `email_manager/models/email_model.py:EmailLog`
is the only model (`claude/maps/models.md#email_manager`) and `EmailLog.send_mail` is the only entry point: callers hand it
a subject, an HTML and a text body, a from-address and a recipient list, and get back whether the message probably went.
The pages under `claude/maps/urls.md#email_manager` are a superuser-only log with a detail row and a bare-HTML view. The
readme is `email_manager/__email_manager_readme.md`.

## Flows

Sending is `email_manager/models/email_model.py:EmailLog.send_mail`. It strips the bodies, then only attempts delivery
when all three of `from_email`, `recipient_list` and `settings.SEND_EMAILS` are truthy. With `allow_users_to_see_others=True`
it is one `django.core.mail.send_mail` with every address in To; otherwise it opens one SMTP connection and sends a separate
`EmailMultiAlternatives` per recipient (HTML attached as the `text/html` alternative) through `connection.send_messages`, so
no recipient sees another. Either way the `EmailLog` row is written afterwards with `recipient_list` joined by `; `,
`probably_sent` set to the backend's return value and `single_email` recording which mode was used. Both sends use
`fail_silently=True`: a dead mail server produces `probably_sent=False`, never an exception in the caller.

The callers are the notification builders and the account flows. `library/log_utils.py:AdminNotificationBuilder.send`
emails every superuser (batch mode) when `is_communication` and `settings.ADMIN_EMAIL_NOTIFICATION` are set;
`snpdb/utils.py:LabNotificationBuilder` emails a lab's address, or the lab's members whose `UserSettings` opted into
discordance or weekly updates, from `settings.DISCORDANCE_EMAIL`; the Monday `discordance-emails-weekly` beat entry in
`variantgrid/celery.py` runs `classification/views/classification_email_view.py:send_summary_emails`, which calls
`send_summary_email_to_user` per opted-in user and skips anyone without an active lab; `user_messages/signal_handlers.py:email_new_message_handler`
emails a new inbox message; `variantgrid/views.py:keycloak_admin` sends the welcome email as the `pre_password_reset`
callback of `library/keycloak.py:Keycloak.add_user`; and MME match notifications go through the same path
(`mme/apps.py:MMEConfig.ready` refuses to start with `MME_ENABLED` unless `MME_FROM_EMAIL` and `SEND_EMAILS` are both set,
precisely because `send_mail` is silently a no-op otherwise). `library/email.py:Email` is the one sender that bypasses the
log - it wraps Django's `send_mail` with `fail_silently=False` and is only built (not sent) in `keycloak_admin`.

Reading the log: `email_manager/views/email_manager_views.py:email_manager_view` mounts a DataTable whose config
`EmailColumns` raises `PermissionDenied` for non-superusers, searches subject and recipients in `power_search`, sorts
newest first and expands each row through `email_detail` (`DatatableConfig._row_expand_ajax`). `email_detail` splits
`recipient_list` back into addresses, matches them to `User.email` and lists the unrecognised remainder; `email_pure` renders
only the body. Both detail templates put the stored HTML in a sandboxed `<iframe srcdoc>` rather than inlining it
(`email_manager/templates/email_detail.html`, `email_manager/templates/email_pure.html`). `EmailLog` is a
`library/preview_request.py:PreviewModelMixin`, so the nightly health check can hover-preview each email:
`eventlog/signals/active_users_health_check.py:email_health_check` counts the window's rows and lists subjects, collapsing to a
count once there are more than ten distinct subjects so Slack is not flooded.

## Why it is shaped this way

The log is written whether or not the send happened, and `settings.SEND_EMAILS` defaults to `False` in
`variantgrid/settings/components/default_settings.py` (only `vgaws.py`, `shariant.py` and `runx1db2.py` turn it on;
`shariantcommon.py` leaves it off so test and demo Shariant never mail real labs). That gives every non-production box a
record of what it would have sent, which is how notification templates are checked without a mail server. Per-recipient
sending is the default because most traffic is discordance and lab notifications where the recipient list is itself
sensitive; batch mode is opt-in and used only for admin broadcasts. The admin `email_manager/admin.py:EmailLogAdmin` denies
add, change and delete so the log stays an audit trail.

## History

The app began as the discordance-notification sender ("Switching over discordance notifications to new framework") and
gained its admin, ID column and the "pure" view as people needed to inspect rendered HTML. The detail and pure views were
moved into a sandboxed iframe (#1521) and then restricted to superusers ("Require super_user to view content of automated
emails"), closing the hole where any logged-in user could read another lab's notifications. `mme/apps.py:MMEConfig.ready`
was added with MatchMaker Exchange, when a silent no-op send became a compliance problem rather than a convenience.

## Traps

`EmailLog.send_mail` returns falsy without raising when `from_email` is `None` - a deployment with `DISCORDANCE_EMAIL` or
`ADMIN_EMAIL_NOTIFICATION` unset logs rows with `probably_sent=False` and nobody is told. `probably_sent` is the backend's
count of messages accepted by the SMTP server, not delivery. `recipient_list` is stored as one `; `-joined text field, so
searching for an address is `icontains`, and `email_detail` re-splits it. The app has no tests directory of its own; the
callers' tests patch `EmailLog.send_mail` (see `mme/tests/test_inbound_notification.py`).
