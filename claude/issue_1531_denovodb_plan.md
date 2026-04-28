# Issue #1531 — Add denovo-db annotation

## Goal

Annotate variants with [denovo-db](https://denovo-db.gs.washington.edu/denovo-db/Download.jsp)
records via VEP `--custom`, exposing three fields on `VariantAnnotation` so users
can see whether (and where) a variant has been reported as a confirmed/observed
de novo in published cohorts.

## Intended workflow

The fields are evidence, not a primary filter:

1. User filters for rare variants (gnomAD AF < threshold).
2. User filters by trio/quad inheritance (variant absent in parents).
3. On the remaining shortlist, user reads denovo-db fields as supporting evidence
   (was this variant seen de novo in a published cohort? in what phenotype?).
4. Optional first-pass coarse filter on the analysis grid: "denovo-db field is
   not null" — i.e. variant present in denovo-db at all.

We deliberately do **not** add a per-variant occurrence count. Recurrence in the
same phenotype shows up directly when users read the `&`-joined values
(e.g. `autism&autism&epilepsy`), and a count divorced from phenotype is more
misleading than useful.

## Why denovo-db (not Gene4Denovo2 or alternatives)

Surveyed alternatives in 2026-04 before committing:

- **[Gene4Denovo2](https://academic.oup.com/nar/article/54/D1/D1069/8268403)** —
  much larger (1.6M DNMs vs ~600k, 130k probands vs ~25k, 96 phenotypes vs
  small handful) and ships GRCh37 + GRCh38 + T2T natively. **Rejected: no
  bulk download** — interactive search and user-data analysis only, so
  unusable as a VEP `--custom` annotation source. Revisit if they ever
  publish a bulk VCF.
- **[Kaplanis 2020 / DDD](https://www.nature.com/articles/s41586-020-2832-5)**
  — 31k DD trios, excellent quality. **Rejected: data only available as
  supplementary TSV/Excel tables to a one-off paper, GRCh37 only, not a
  maintained resource.** Already incorporated into denovo-db / Gene4Denovo2.
- **Genomics England 100kGP de novo dataset.** **Rejected: only accessible
  inside the GE Research Environment**, no redistribution.
- **DECIPHER research track, SFARI Gene / GPF.** **Rejected: no clean bulk
  VCF**; UI / interactive access only.

denovo-db is stale (frozen at v1.6.1 / 2019, no DNMs published after 2018) but
it's the only maintained public DNM database with a bulk VCF download, which
makes it the only viable choice for VEP `--custom` annotation today.

## Source data

- Source: denovo-db v1.6.1, **SSC + non-SSC merged** (full dataset).
- denovo-db only publishes a GRCh37 / hg19 release. No GRCh38 build.

**The upstream VCF is not directly usable.** Inspecting the actual VCF INFO
fields (`SAMPLE_CT`, `DBSNP`, `AA`, `ESP_AA_AF`, `ESP_EA_AF`, `EXAC_AF`,
`G1K_ALLELE_CT`, `CADD`, `LOF`, `LRT`, `STUDIES`, `TRANSCRIPTS`) — the
clinically useful fields **`StudyName`, `PubmedID`, `PrimaryPhenotype` are
only in the TSV, not the VCF**. Most VCF INFO fields duplicate annotation we
already have from VEP / dbNSFP / gnomAD. So we build our own VCF from the
upstream TSV.

**Pipeline (run once per denovo-db release, by deployer):**

1. Download SSC and non-SSC TSVs from the denovo-db site.
2. Run `annotation/annotation_data/generate_annotation/denovo_db/denovo_db_generate_grch37.sh`,
   which calls `denovo_db_tsv_to_vcf.py` to convert + aggregate per-variant
   (`&`-joined parallel arrays of StudyName / PubmedID / PrimaryPhenotype +
   `SAMPLE_CT` count), then `bcftools sort -Oz` to produce
   `denovo-db.variants.v.1.6.1.vcf.gz` in `annotation_data/GRCh37/`.
3. Run `denovo_db_liftover_grch37_to_grch38.sh` (uses `bcftools +liftover` with
   `GRCh37_to_GRCh38.chain.gz`) to produce
   `denovo-db.variants.v.1.6.1.GRCh38.vcf.gz` in `annotation_data/GRCh38/`.

Resulting deployed files:
- GRCh37: `annotation_data/GRCh37/denovo-db.variants.v.1.6.1.vcf.gz`
- GRCh38: `annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.vcf.gz`

Both pipeline scripts live in `annotation/annotation_data/generate_annotation/denovo_db/`
and are checked in alongside the existing gnomAD / VEP fasta generators.
- T2T-CHM13v2.0: not supported initially (settings `None`). Could be added
  later by running `bcftools +liftover` GRCh38 → T2T with the existing
  `hg38ToHs1` chain.
- **SSC licence note.** The denovo-db Usage page states: *"The use of Simons
  Simplex Collection (SSC) and Simons VIP data sets is limited to projects
  related to advancing the field of autism and related developmental disorder
  research."* The decision is to ship SSC by default — primary deployments
  (SA Pathology autism diagnostics, neurogenetics) fall within that scope.
  Other deployments that don't want SSC content can swap to the non-SSC file
  via env-specific settings (the `denovo_db` settings path is per-deployment).
  Citation requirement (always): *"denovo-db, Seattle, WA
  (URL: denovo-db.gs.washington.edu) [date (month, yr) accessed]."*

## Fields to annotate

Three TextFields on `VariantAnnotation`, each storing the raw VEP `&`-joined
values when a variant matches multiple denovo-db records:

| VariantAnnotation column        | denovo-db INFO field | Notes                                                  |
| ------------------------------- | -------------------- | ------------------------------------------------------ |
| `denovodb_study_name`           | `StudyName`          | Identifies cohort (e.g. `AutismKaplanis`).             |
| `denovodb_pubmed_id`            | `PubmedID`           | Kept distinct from existing VEP-derived `pubmed`.      |
| `denovodb_primary_phenotype`    | `PrimaryPhenotype`   | e.g. `autism`, `epilepsy`, `intellectualDisability`.   |

Fields **deliberately skipped** (and why):

- `NumProbands` / `NumControls` — study-level totals, not per-variant counts;
  redundant under the workflow above. Users can look up cohort sizes from the
  StudyName / PubmedID if needed.
- `SampleID`, `SequenceType`, `Validation` — low filter and evidence value.
- Gene / transcript / functional annotations / prediction scores — VEP + dbNSFP
  already cover these better.
- A derived per-variant occurrence count — see "Intended workflow" above.

## Implementation

The custom-track machinery is fully data-driven: adding a `VEPCustom` enum
entry, three rows in `vep_columns.py`, three model fields, a settings path, and
a migration is enough. `annotation/vep_annotation.py` already iterates over
`VEPCustom` to build `--custom` arguments and will pick up the new track
automatically (see `vep_annotation.py:251-277`).

### 1. `annotation/models/models_enums.py` — `VEPCustom`

Add to the `VEPCustom(models.TextChoices)` enum:

```python
DENOVO_DB = 'D', 'denovo_db'
```

The `label` (`'denovo_db'`) is used both as the VEP `--custom` prefix and,
lower-cased, as the settings key (see `BulkVEPVCFAnnotationInserter` and
`vep_annotation.py:269`).

### 2. `annotation/vep_columns.py` — column definitions

Add a thin helper alongside the existing `_gnomad3` / `_uk10k`-style helpers:

```python
def _denovo_db(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ denovo-db custom track, GRCh37 + GRCh38 (locally lifted-over),
        standard pipeline. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.LITERATURE,
        "vep_custom": VEPCustom.DENOVO_DB,
        "source_field_has_custom_prefix": True,
        "genome_builds": GRCH37_38,
        "pipeline_types": STANDARD,
        **overrides,
    })
```

Add three entries inside the `VEP_COLUMNS` tuple, in a new
`# ---------- denovo-db (GRCh37) -------------------------` section:

```python
_denovo_db('StudyName',        'denovodb_study_name'),
_denovo_db('PubmedID',         'denovodb_pubmed_id'),
_denovo_db('PrimaryPhenotype', 'denovodb_primary_phenotype'),
```

Should we want to gate this behind a columns_version bump (recommended — see
section 7), add `min_columns_version=N` to the helper defaults.

### 3. `VariantAnnotation` model — three new fields

In `annotation/models/models.py`, alongside the existing literature-style
TextFields (e.g. near `repeat_masker`):

```python
denovodb_study_name = models.TextField(null=True, blank=True)
denovodb_pubmed_id = models.TextField(null=True, blank=True)
denovodb_primary_phenotype = models.TextField(null=True, blank=True)
```

These store raw `&`-joined values straight from VEP. No formatter is needed —
the existing fall-through in `BulkVEPVCFAnnotationInserter.vep_to_db_dict()`
copies the source value verbatim when no entry in `self.field_formatters`
matches.

These belong on the **representative** record (`VariantAnnotation`), not per
transcript (`VariantTranscriptAnnotation`) — denovo-db is a variant-level
attribute. Confirm by leaving them off `VariantTranscriptAnnotation`; the
`variant_only_columns` derivation will then exclude them from the per-transcript
table automatically.

### 4. `BulkVEPVCFAnnotationInserter` — no code changes

The inserter is fully data-driven from `vep_columns.filter_for(...)`. With no
matching entry in `field_formatters`, the raw `&`-joined string is copied
straight into the TextField column. The custom-track skip path in
`_add_vep_field_handlers` (lines ~253-263) silently drops the track on builds
where settings have it as `None`, so GRCh38 / T2T need no special handling.

### 5. `variantgrid/settings/components/annotation_settings.py`

Under `vep_config` for each build:

```python
# GRCh37  (combined SSC + non-SSC file, see "Source data" above)
"denovo_db": "annotation_data/GRCh37/denovo-db.variants.v.1.6.1.vcf.gz",

# GRCh38  (locally lifted-over with `bcftools +liftover`)
"denovo_db": "annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.vcf.gz",

# T2T-CHM13v2.0
"denovo_db": None,   # not supported initially; can be added via GRCh38 → T2T liftover
```

The settings key must be the lower-cased `VEPCustom.DENOVO_DB.label`
(`"denovo_db"`).

### 6. Track the source filename / version on `VariantAnnotationVersion`

`VariantAnnotationVersion` records the version string for every annotation
data source it pulled in (see `models.py:611-644`):

```python
gnomad   = models.TextField(blank=True, null=True)   # e.g. "4.1"
cosmic   = models.IntegerField(blank=True, null=True)
dbnsfp   = models.TextField(blank=True, null=True)
hgmd     = models.TextField(blank=True, null=True)
...
```

denovo-db needs the same treatment so we can tell, for any annotated variant,
*which* denovo-db release the `denovodb_*` fields came from (and so we can
distinguish SSC-included vs non-SSC-only deployments at the data layer rather
than relying on settings hygiene).

Add to `VariantAnnotationVersion`:

```python
denovo_db = models.TextField(blank=True, null=True)   # e.g. "1.6.1"  or  "1.6.1-ssc"
```

Convention: store the upstream version, with a `-ssc` suffix if the SSC file
was included (default for our deployments) and no suffix for non-SSC-only.
This makes the choice queryable (e.g. `qs.filter(denovo_db__contains='ssc')`)
without inventing a separate boolean.

Population: `annotation_settings.py` will need a parallel
`"denovo_db_version": "1.6.1-ssc"` (or similar) entry, copied into the new
`VariantAnnotationVersion` row when annotation runs are created — same pattern
as `gnomad`, `dbnsfp`, etc. Find the existing version-string population point
(grep for where `vav.gnomad` / `vav.dbnsfp` get set during annotation-run
creation) and extend it.

### 7. Migration

New migration in `annotation/migrations/` (next number: `0131_…`). Three
changes:

1. `AddField` × 3 for the new `VariantAnnotation` columns (all
   `TextField(null=True, blank=True)`).
2. `AddField` × 1 for `VariantAnnotationVersion.denovo_db` (TextField, null,
   blank).
3. Likely a `VariantAnnotationVersion.columns_version` bump if we want to mark
   data as needing re-annotation — see section 8.

Because `VariantAnnotation` is partitioned (`psqlextra`), follow the same
migration patterns used by recent column additions (see e.g.
`0127_manualvariantentry_warning_message.py` neighbours and the dbNSFP /
gnomAD-era migrations for VariantAnnotation column adds).

### 8. `columns_version` bump (decision needed)

`VariantAnnotationVersion.columns_version` currently caps at `3`. Two options:

- **A. Bump to 4 and gate the new fields on `min_columns_version=4`.** Cleaner
  semantically: existing v3 annotation data is left untouched, and new
  deployments / re-runs at v4 pick up denovo-db. Requires updating
  `VariantAnnotationVersion.get_variant_annotation_columns()` (see
  `models.py:673-690`) to include the three new fields when
  `columns_version == 4`. Also touch `annotation_settings.py` defaults if we
  want new deployments at v4 by default.
- **B. Add fields without bumping.** Existing v3 annotations would silently get
  `NULL` for these fields until re-annotated. Simpler, but obscures which
  annotation runs actually have denovo-db data.

Recommend **A**. Confirm with @davmlaw before implementing.

### 9. VEP custom-track wiring — automatic

`vep_annotation.get_vep_command()` already loops `for vep_custom in VEPCustom`
and emits `--custom file=…,short_name=denovo_db,format=vcf,fields=StudyName%PubmedID%PrimaryPhenotype`
based on the entries in `vep_columns.py`. No code change needed here. Verify
behaviour by running `python manage.py shell` and inspecting
`get_vep_command(...)` output for a GRCh37 standard pipeline annotation run.

### 10. Variant details page

Out of scope for the inserter change, but called out so it's not forgotten:

- Split each of the three fields on `&`.
- Render as a small table: one row per denovo-db record, columns = StudyName,
  PrimaryPhenotype, PubMed link.
- Cardinality of the three split lists is always equal (one per denovo-db
  record) — they're parallel arrays produced by VEP.
- Locate alongside other "external evidence" sections on the variantopedia
  variant detail template.

Track as a follow-up — separate PR is fine.

### 11. Analysis grid / column registry

Three new columns auto-flow through the existing `VariantGridColumn`
infrastructure once they're in `VEP_COLUMNS` and on the model. Confirm they
appear in the column picker under the **Literature** category
(`ColumnAnnotationCategory.LITERATURE`), and that "is not null" filtering works
out of the box (it should — they're plain TextFields).

## Testing

- **Unit test in `annotation/tests/`**: feed a synthetic VEP-annotated VCF
  through `BulkVEPVCFAnnotationInserter` with a CSQ entry containing
  `denovo_db_StudyName=AutismFoo&AutismBar`,
  `denovo_db_PubmedID=12345&67890`,
  `denovo_db_PrimaryPhenotype=autism&autism`,
  and assert the three TextFields land on the resulting `VariantAnnotation`
  unchanged (still `&`-joined).
- **T2T smoke test**: with `denovo_db=None` in settings, ensure the inserter
  logs the "Skipping custom denovo_db due to missing settings" warning
  (existing path) and produces no errors.
- **VEP command test**: assert `--custom` for `denovo_db` is present in
  GRCh37 and GRCh38 commands and absent from T2T. Pattern: see existing tests
  covering `repeat_masker` / `uk10k`.

## Out of scope (note for future work)

- T2T-CHM13v2.0 support (would need a second `bcftools +liftover` step
  GRCh38 → T2T using the existing `hg38ToHs1` chain).
- Generating the lifted-over GRCh38 VCF — performed by the user/deployer with
  `bcftools +liftover`; this codebase only consumes the resulting file.
- Per-variant occurrence count (`denovodb_proband_count`) — only add if a real
  recurrence-ranking use case shows up.
- Variant details page rendering (separate follow-up).

## Open questions

1. **columns_version bump**: option A (bump to 4) or B (no bump)? Recommend A.
2. Exact format of the `VariantAnnotationVersion.denovo_db` string —
   `"1.6.1-ssc"` vs `"1.6.1 (SSC)"` vs separate `denovo_db_includes_ssc`
   boolean. Match whatever stylistic convention is closest to existing
   `gnomad` / `dbnsfp` strings on `VariantAnnotationVersion`.
3. Whether to merge SSC + non-SSC into a single VCF at deployment time, or
   ship two `--custom` tracks. Single merged file is simpler (one VEP custom
   track, one settings entry) and is the recommended approach.
