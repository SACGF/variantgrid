# pathtests — research notes

Verified against a96540a68 on 2026-09-25

pathtests holds two loosely related things: **pathology tests** - a named, versioned gene list that a lab curator
maintains and that analyses filter on - and **cases / test orders**, the clinical-request side that an external LIMS
populates. It is off by default (`PATHOLOGY_TESTS_ENABLED = False`), disabled on Shariant, and has no tasks, commands
or signals of its own; the only writers of cases and orders live in the private SA Pathology `sapath` app. Fields and
URLs: [models](../maps/models.md#pathtests), [urls](../maps/urls.md#pathtests).

## Flows

### Creating and curating a test

`pathtests/views.py:manage_pathology_tests` creates a `pathtests/models.py:PathologyTest` (the name is the primary key)
with the requesting user as curator, and a v1 `pathtests/models.py:PathologyTestVersion` whose gene list is a clone of a
chosen gene list or of another test version's list, or a new empty list, put in the `PathologyTest` gene list category
(`genes/models/models_gene_list.py:GeneListCategory.get_pathology_test_gene_category`). A version is a draft until
`confirmed_date` is set; `can_confirm` also needs an enrichment kit. Confirming happens in
`pathtests/views.py:view_pathology_test_version` (curator-only POST), which stamps the date, and
`PathologyTestVersion.save` then locks the gene list and `set_as_active_test` repoints
`pathtests/models.py:ActivePathologyTestVersion`. "Active" therefore means "most recently confirmed", and the REST view
`pathtests/views_rest.py:PathologyTestLatestVersionView` ("latest confirmed version") simply returns the active one.
`PathologyTestVersion.save` also refuses (after the row is written) a gene list whose import did not succeed, and touches
the parent test so its `modified` reflects version changes.

### Gene modification requests

Anyone can file a request to add or remove a gene: the gene grid JS (`variantgrid/static_files/default_static/js/gene_grid.js`)
POSTs to `pathtests/views.py:modify_pathology_test_version`, creating a PENDING
`pathtests/models.py:PathologyTestGeneModificationRequest`. The curator reviews them on the version page:
`pathtests/views.py:get_gene_modification_request` splits pending requests into additions and deletions **by whether
the gene is currently in the list**, not by the request's `operation`, and the form posts `add-<gene>` / `del-<gene>`
radios (ignore / reject / accept). If the version is already confirmed and anything is accepted,
`PathologyTestVersion.next_version` clones it (gene list cloned and unlocked, version +1, unconfirmed), all reviewed
requests move to the new draft, and the genes are changed there via
`genes/models/models_gene_list.py:GeneList.add_and_remove_gene_symbols` with a provenance note per added gene
(`pathtests/views.py:handle_modification_requests`). The old version stays active until the new one is confirmed. This is
the only UI route to a new version.

### Where tests are consumed

`analysis/models/nodes/filters/gene_list_node.py:GeneListNode` has a `pathology_test_version` FK and warns when the
version is older than `PATHOLOGY_TEST_STALE_WARNING_DAYS`. The gene grid adds a column per test version through the two
REST endpoints. `pathtests/templatetags/enrichment_kit_coverage_tags.py:enrichment_kit_coverage` shows, per version,
which genes each kit in `PATHOLOGY_TEST_SORTED_ENRICHMENT_KITS` misses.

### Cases and orders

`pathtests/models.py:Case` and `pathtests/models.py:PathologyTestOrder` are
`patients/models.py:ExternallyManagedModel`s: each carries a one-to-one `external_pk`, and
`ExternallyManagedModel.can_write` is false when that key's external manager says so - `pathtests/forms.py:CaseForm`
disables every field in that case. There is no create or edit view in this app: SA Pathology's Helix import makes a Case
per Helix accession and a PathologyTestOrder per SAP order number, and links samples to patients through those cases.
The pages here are read-only views, the three datatables in `pathtests/grids.py`, lookup-by-LIMS-id redirects
(`view_external_case`, `view_external_pathology_test_order`) and "my cases" by lead scientist, which includes the
scientists a user follows (`patients/models.py:get_lead_scientist_users_for_user`, toggled by
`pathtests/views.py:follow_scientist`). `pathtests/models.py:get_external_order_system_last_checked` reaches into
`sapath` for the "last checked" time, and returns None when that app is absent.

## Why it is shaped this way

- **Name as primary key** lets the gene grid and API address a test by the name clinicians use
  (`api/view_latest_pathology_test_version/<name>`), and lets `sapath`'s importer `get_or_create` by name.
- **Confirmed versions are immutable** (locked gene list) because orders and analyses reference a specific version; a
  change after confirmation must be a new version so historical results keep the list they were run against.
- **ActivePathologyTestVersion is a separate row** so switching the active version is one update, and deleting a
  test can drop "active" without touching versions.

## History

The app dates from 2020 (migrations stop in 2021). Clinician pages and SA Path test-form generation were removed in
2026-08 (SACGF/variantgrid_sapath#434), the grids moved from jqGrid to datatables (#1785), and the 2026-09 dead code sweep
(#1795) trimmed it further. `PathologyTestOrderSample`, `PathologyTestOrderPopulation`, `CaseClinician` and
`Case.workflow_status` are schema only - nothing (sapath included) writes them. The order's library/sequencing
timestamps are copied from Helix by sapath.

## Traps

- `pathtests/models.py:PathologyTestGeneModificationRequest.__str__` does `" ".join(...)` over model instances and raises
  `TypeError`, so the admin changelist for that model (registered in `pathtests/admin.py`) errors.
- `pathtests/views.py:manage_pathology_tests` clones a source gene list without resetting `locked`. A test created from
  a *confirmed* version gets a locked v1 list, and accepting any request on that draft fails in
  `GeneList.add_and_remove_gene_symbols` with PermissionDenied. `next_version` unlocks; this path does not.
- Deleting a test (`pathtests/views.py:view_pathology_test`, type "delete") removes its ActivePathologyTestVersion;
  restoring does not put it back, and a confirmed version cannot be re-confirmed, so a restored test has no active
  version (API 404) until a gene request is accepted to make a new one, or an admin fixes it.
- Deletion vs addition is inferred from list membership, not `operation`: a REMOVE request for a gene not in the list
  is shown and applied as an addition.
- A request filed after the curator loaded the page has no radio in the POST, so `handle_modification_requests` saves it
  with `outcome=None` and hits the NOT NULL constraint.
- `PathologyTestVersion.next_version` and `GeneList.clone` turn `self` into the copy (`copy = self; copy.pk = None`).
  Re-fetch if you still need the original.
- `PathologyTestVersion.replace_gene_list` has no callers; it overwrites the new list's row with the old list's fields.
- Curator checks are an exact user match (no superuser override) and raise `PermissionError`, which is a 500, not a 403.
- `pathtests/models.py:cases_for_user` always returns an empty queryset (disabled TODO), so the variantopedia dashboard's
  `user_has_cases` is always false.
