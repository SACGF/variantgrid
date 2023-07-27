## Notification Builder

NotificationBuilder can be used to send notifications to HTML Email and/or Slack.
Generally provide it wirth Markdown (which can then be converted to basic HTML if needed).

Instead of using NotificationBuilder directly, consider:

AdminNotificationBuilder: Use this to notify administrators of things such as health checks, important events etc.

LabNotificationBuilder: (in the snpdb project) Use this to notify a lab via their preferred notification method of
discordances or other lab specific events.


## jQGrid

This is a way of providing paginated data using jQuery's jQGrid. This has been deprecated in favour of using
DataTables, but many jQGrid's still exist.

## Utils directory
class_utils
collection_utils
color_utils
etc

Just handy utils for common functions e.g. color_utils.rgb_invert, 