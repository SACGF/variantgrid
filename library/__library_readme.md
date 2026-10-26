## Grids

Every grid on the site renders with DataTables off one server side engine,
`snpdb.views.datatable_view.DatatableConfig`, served by `DatabaseTableView`.

* `library.django_utils.filter_rules` - the column filter rule vocabulary the grids' filter builder
  and FilterNode share, and the rule -> Q conversion both go through.
* `library.django_utils.grid_export` - streaming CSV for the server side download.

## Utils directory
class_utils
collection_utils
color_utils
etc

Just handy utils for common functions e.g. color_utils.rgb_invert, 