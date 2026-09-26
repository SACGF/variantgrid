# pedigree — research notes

Verified against a96540a68 on 2026-09-25

A Pedigree is the general-shape family: a PED file family laid over a Cohort, which the analysis `PedigreeNode` filters
by an inheritance model. Trio / Duo / Quad (`snpdb/models/models_family.py`) are the fixed-shape alternatives and get
far more use; this app is small (under 500 lines) and quiet. Fields, URLs and settings are in the maps
([models](../maps/models.md#pedigree), [urls](../maps/urls.md#pedigree)); there are no tasks or commands of its own,
and the one signal receiver is the search hook `pedigree/signals/pedigree_search.py:search_pedigree` (see
[signals](../maps/signals.md)).

## Flows

### Upload a PED file

A `.ped` upload goes through the upload app: `upload/tasks/import_ped_task.py:ImportPedTask` calls
`pedigree/ped/import_ped.py:import_ped`, which reads the first six whitespace-separated columns with pandas (indexed on
family + individual), creates the `PedFile` and grants `view_pedfile` to the uploader and each of their groups, then per
family builds a dependency graph of child -> parents and saves `PedFileRecord`s in `toposort` order, so a parent row
always exists before a child's `father` / `mother` FK points at it. Value parsing (`0`/`.`/`Unknown` parents, `1/2/M/F`
sex, `1/2/U/A/0/X/-9` affection) is in `pedigree/ped/ped_file_utils.py`. Each family is then checked with
`pedigree/models.py:validate` (at least one record, fathers male, mothers female, at least one affected) and a failure
raises, marking the `PedFile` ERROR. Only after a successful import does the task create the `UploadedPedFile` link back
to the original file.

### Link a family to samples

`pedigree/models.py:CohortSamplePedFileRecord` maps each PED record to a `CohortSample`. Three ways to get them:

- **On upload**, `pedigree/ped/import_ped.py:automatch_pedigree_samples` runs when
  `PEDIGREE_MIN_COHORT_SAMPLE_MATCHES_FOR_AUTO_MATCH` is set (default 3): every cohort the user can see with that many
  sample names equal to PED individual IDs gets a Pedigree via `pedigree/models.py:create_automatch_pedigree`.
- **By hand**, the pedigrees page / PED file page pick a cohort and family and hit
  `pedigree/views.py:create_pedigree_from_cohort_and_ped_file_family`, which also auto-matches by name.
- **Then editing** on `pedigree/views.py:view_pedigree`: a formset of record -> cohort sample; saving deletes and
  recreates every mapping. `pedigree/forms.py:BaseCohortSamplesForPedFileRecordsFormSet` stops one sample being used twice.

Records left unmapped are allowed; nothing enforces "one mapping per record".

### Filter in an analysis

`analysis/models/nodes/sources/pedigree_node.py:PedigreeNode` is an `AbstractCohortBasedNode` over `pedigree.cohort`.
With `inheritance_model` set it builds a per-sample zygosity regex over the cohort genotype packed string
(`snpdb/models/models_cohort.py:CohortGenotypeCollection.get_zygosity_q`) from the mapped records only: recessive =
affected HOM_ALT, unaffected HET; dominant = affected HET/HOM_ALT, unaffected MISSING. `require_zygosity=False` also lets
a no-call through. Any inheritance model defeats the label-count cache (`_has_filters_that_affect_label_counts`).

### Chart

`pedigree/views.py:pedigree_chart` renders through `snpdb.graphs.graphcache` to
`pedigree/graphs/pedigree_chart.py:PedigreeChart`, which converts the *original uploaded file* with `ped_parser
--to_madeline` then runs Madeline2 (`PEDIGREE_MADELINE2_COMMAND`, off by default - the view hides the chart when unset).
`get_ped_parser_command` runs the `ped_parser` script with the current interpreter because it is not on `PATH` for a
uv venv under systemd (#1572). `variantgrid/deployment_validation/tool_version_checks.py` checks both tools.

## Why it is shaped this way

- **PedFile is kept separate from Pedigree** so one uploaded family can be laid over several cohorts (the same family
  sequenced twice, or re-called after a VCF replace - `snpdb/models/vcf_replace_data.py` repoints pedigrees and their
  mappings at the new cohort).
- **Affection is a nullable boolean**, collapsing PED's several encodings (and Phenotips' `-9`) to affected /
  unaffected / unknown.
- **Export lives here but is used by snpdb**: `pedigree/ped/export_ped.py` writes the trio and unrelated PED files that
  `snpdb/models/models_somalier.py` passes to `somalier relate --ped`.

## History

- 2020-12 (#46): Somalier work added `SomalierPedigreeRelate`, which was never wired up (below).
- 2026-03 (#1460, 879657742): pedigree grids moved from jqGrid to DataTables; a unit-testing pass fixed parser bugs
  (the trio export wrote affection as 1/0 rather than PED's 2/1, so somalier read the
  affected proband as unaffected; blank parent cells and lower-case sex now parse).
- 2026-06 (#1572): `ped_parser` invoked through the interpreter.
- 2026-09 (#183): somalier 0.3.5 support in snpdb; the pedigree subclass was not touched.

## Traps

- **Viewers can edit a pedigree.** `pedigree/views.py:view_pedigree` saves the POSTed form and formset without checking
  `can_write` - the template only hides the Save button. `pedigree/forms.py:PedigreeForm` is `ALL_FIELDS` with
  unfiltered querysets, so the POST can also repoint `cohort` / `ped_file_family` at objects the user cannot see.
- **No permission check on the family.** `pedigree/views.py:create_pedigree_from_cohort_and_ped_file_family` fetches the
  `PedFileFamily` with `get_object_or_404`, so any user can build a pedigree over another user's PED family by id and then
  read its records on the pedigree page. It is also a state-changing GET.
- **Unknown affection is "unaffected" in the node but not in `get_samples`.**
  `PedigreeNode.get_affected_unaffected_sample_zygosities_dict` tests `if ped_file_record.affection`, so `None` gets the
  unaffected zygosities; `pedigree/models.py:Pedigree.get_samples` with `affected=False` excludes `None`.
- **Dominant unaffected only matches MISSING (`.`), not HOM_REF (`R`).** An unaffected sample with an explicit 0/0 call
  (a gVCF-style or joint-called VCF writing `R`) excludes the variant. The trio code uses
  `analysis/models/nodes/family_inheritance.py` `NO_VARIANT = {MISSING, HOM_REF}`.
- **Recessive treats every unaffected member as a carrier** (HET required), so an unaffected HOM_REF sibling drops a real
  recessive hit unless `require_zygosity` is off - and even that only admits no-calls.
- **One bad family aborts the rest.** `import_ped` raises on the first invalid family: earlier families are already
  saved, later ones never are, and a family with no affected member (or a parent ID not in the family - a `KeyError`,
  or a parent cycle - `toposort` error) fails the whole upload.
- **Auto-match is not idempotent.** Re-uploading the same PED file creates another PedFile and another auto Pedigree per
  matching cohort.
- **The chart needs an upload.** `PedigreeChart.save` reads `ped_file.uploadedpedfile`, so a PedFile made any other way
  (tests, shell) cannot be charted; it charts every family in the file, and its cache key is the PedFile id only.
- **`SomalierPedigreeRelate` is dead.** `pedigree/models.py:SomalierPedigreeRelate` is never created,
  its `write_ped_file` is `pass`, and it inherits `has_ped_file() -> False`; it is also missing from the pre_delete cleanup
  receiver in `snpdb/models/models_somalier.py`. Wiring it up needs a real `write_ped_file` from the mapped records.
- **`PedFileFamily.errors` re-queries** its records on every access and `is_valid` calls it again.
