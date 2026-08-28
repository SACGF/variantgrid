## Notification Builder

NotificationBuilder can be used to send notifications to HTML Email and/or Slack.
Generally provide it wirth Markdown (which can then be converted to basic HTML if needed).

Instead of using NotificationBuilder directly, consider:

AdminNotificationBuilder: Use this to notify administrators of things such as health checks, important events etc.

LabNotificationBuilder: (in the snpdb project) Use this to notify a lab via their preferred notification method of
discordances or other lab specific events.


## Grids

Every grid on the site renders with DataTables. There are two server side implementations:

* `snpdb.views.datatable_view.DatatableConfig` - use this for new grids.
* `library.jqgrid.JqGrid` - the older engine, still driving the variant grids (which build their
  columns dynamically from `CustomColumn`) and the CSV/VCF export. `JqGridDatatableView`
  (`library/django_utils/jqgrid_datatable_adapter.py`) serves it to the DataTables client.

## Utils directory
class_utils
collection_utils
color_utils
etc

Just handy utils for common functions e.g. color_utils.rgb_invert, 