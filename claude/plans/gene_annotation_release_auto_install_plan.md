# Gene Annotation Release auto-install plan

## Goal

Remove the manual work in getting a `GeneAnnotationRelease` that matches the VEP a
deployment is actually running.

Today, after `create_new_variant_annotation_version` makes a new
`VariantAnnotationVersion` (VAV), an admin has to:

1. Read the VEP version strings off the annotation build detail page.
2. Work out which cdot per-GFF file corresponds to them.
3. Find and download that file by hand.
4. Run `import_gene_annotation --json-file=... --release=...`, inventing a release name.
5. Go into Django Admin and set `VariantAnnotationVersion.gene_annotation_release`,
   because the auto-link never fires.

Two changes fix this:

- **Auto-link** an existing `GeneAnnotationRelease` when one already matches the new
  VAV's config. This is the common case on GRCh37, whose RefSeq annotation does not
  change between VEP releases.
- **Auto-install** when one does not exist: a single command that resolves the right
  cdot GitHub release asset, downloads it, imports it as a `GeneAnnotationRelease`,
  links it to the VAV, and runs gene annotation.
  `create_new_variant_annotation_version` prints that command when it cannot link.

## Background: the facts this design rests on

### cdot publishes per-GFF assets with a parseable naming scheme

The current cdot data release (`data_v0.2.33`) has, alongside the combo files:

```
cdot-0.2.33.Homo_sapiens_GRCh37_Ensembl_87.gtf.json.gz
cdot-0.2.33.Homo_sapiens_GRCh37_RefSeq_105.20201022.gff.json.gz
cdot-0.2.33.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_112.gtf.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_113.gtf.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_114.gtf.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_115.gtf.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_116.gtf.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_110.gff.json.gz
cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz
```

Shape:

```
cdot-{data_version}.Homo_sapiens_{genome_build}_{RefSeq|Ensembl}_{release_token}.{gff|gtf}.json.gz
```

The genome build segment matches `GenomeBuild.name` verbatim. The release token can
contain underscores (`RS_2025_08`), so parsing anchors on the consortium segment.

T2T-CHM13v2.0 has combo files only — no per-GFF asset. Auto-install cannot serve
T2T; auto-link still can.

### The release token is derivable from the VAV

Verified against the live VAV rows:

| build | consortium | source field | value | token |
|---|---|---|---|---|
| GRCh37 | RefSeq | `refseq` | `105.20220307 - GCF_000001405.25_GRCh37.p13_genomic.gff` | `105.20220307` |
| GRCh38 | RefSeq | `refseq` | `GCF_000001405.40-RS_2023_10 - GCF_000001405.40_GRCh38.p14_genomic.gff` | `RS_2023_10` |
| GRCh38 | RefSeq | `refseq` | `GCF_000001405.40-RS_2025_08 - GCF_000001405.40_GRCh38.p14_genomic.gff` | `RS_2025_08` |
| GRCh38 | RefSeq | `refseq` | `110 - GCF_000001405.40_GRCh38.p14_genomic.gff` | `110` |
| GRCh37 | Ensembl | fixed | last GRCh37 Ensembl release | `87` |
| GRCh38 | Ensembl | `vep_cache` or `vep` | `116` | `116` |
| T2T-CHM13v2.0 | Ensembl | `genebuild` | `2022-06` | `2022_06` |

RefSeq rule: take the part of `refseq` before `" - "`, then strip a leading
`GCF_\d+\.\d+-` prefix.

One legacy GRCh37 RefSeq value is a timestamp rather than a release
(`2020-10-26 17:03:42 - GCF_000001405.25_GRCh37.p13_genomic.gff`); the existing code
in `_gene_annotation_release_and_gff_url` maps it to `105.20201022`, so keep that
special case.

### The existing GFF URL derivation is wrong for every current build

`VariantAnnotationVersion._gene_annotation_release_and_gff_url`
(`annotation/models/models.py:885`) guesses a URL, and `import_gene_annotation.py:176`
string-compares it against the URL cdot recorded, to decide whether to auto-link.
Neither side agrees, which is why the auto-link never fires:

| | derived by our code | actually in the cdot file |
|---|---|---|
| RefSeq RS_2025_08 | `http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz` | `https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz` |
| RefSeq 110 | as above pattern | `https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz` |
| RefSeq GRCh37 105.20220307 | as above pattern | `https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz` |
| Ensembl 116 | `ftp://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz` | `https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` |

cdot moved to `genomes/all/annotation_releases/9606/...` for RefSeq and uses GTF over
HTTPS for Ensembl. The NCBI layout also differs by release generation: a `GCF_`-prefixed
release directory has the GFF directly inside it, while a bare numeric release
(`110`, `105.20220307`) has an extra assembly directory whose name is the GFF filename
with `_genomic.gff` removed.

### `GeneAnnotationRelease.version` is not a usable match key

Live values are ad-hoc and inconsistent — hand-typed as `--release=` arguments over
several years:

```
refseq_grch38_110   refseq_grch37_105   GRCh38_RS_2023   GRCh37_105
T2Tv2_Ensembl_2022_06   grch38_2024_08   GRCh38_RefSeq_2025_08
```

`gene_annotation_import.url` is the reliable identity: cdot writes it into every
transcript record, so it is consistent across deployments and release generations.
Match on it.

Corollary worth noting: GRCh37 VAVs 4 and 6 report `refseq=105.20220307` but are
linked to a release whose URL is `105.20201022`. The current hand-linking has already
gone wrong; a URL-keyed match surfaces that rather than hiding it.

### `import_cdot_latest` is a superset of the per-GFF files

The combo file for a build contains the latest copy of every transcript, including
those in the individual per-GFF files. So creating a release does not need a targeted
backfill of missing transcript versions — the combo import covers it.

`import_cdot_latest` also short-circuits when the installed cdot data version already
matches GitHub's latest, so chaining it as a prerequisite step is cheap. It currently
does that via `exit(0)` (`genes/management/commands/import_cdot_latest.py:68`), a hard
process exit that would kill a chained run before the release import.

## Design

### 1. `genes/cdot_data_release.py` (new)

Thin layer over `cdot.data_release` for per-GFF assets, kept separate from the
management commands so it is unit-testable without a DB.

```python
CDOT_RELEASE_ASSET_PATTERN = re.compile(
    r"^cdot-(?P<cdot_version>\d+\.\d+\.\d+)\.Homo_sapiens_"
    r"(?P<genome_build>.+?)_(?P<annotation_consortium>RefSeq|Ensembl)_"
    r"(?P<release>.+)\.(?P<file_type>gff|gtf)\.json\.gz$"
)


@dataclass(frozen=True)
class CdotReleaseAsset:
    cdot_version: str
    genome_build: str
    annotation_consortium: str
    release: str
    file_type: str
    url: str


def get_latest_release_assets() -> list[CdotReleaseAsset]
def find_latest_release_asset(genome_build_name, annotation_consortium_label, release) -> Optional[CdotReleaseAsset]
```

`get_latest_release_assets` calls `cdot.data_release.get_latest_data_release()` (which
already filters GitHub releases to ones whose data schema matches the installed cdot),
then keeps the assets matching the pattern. Build and consortium comparisons are
case-insensitive, matching the style of `get_latest_combo_file_urls`.

### 2. `VariantAnnotationVersion` — replace URL guessing with token derivation

In `annotation/models/models.py`:

- Add `cdot_gene_release_token` — the release token per the table above. Returns None
  when it cannot be determined.
- Rewrite `_gene_annotation_release_and_gff_url` to emit the URLs cdot actually
  records, so URL comparison can work:
  - RefSeq: `https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/{release_dir}/{assembly_dir/}{gff_filename}.gz`,
    where `release_dir` and `gff_filename` come straight out of the `refseq` field's two
    halves, and `assembly_dir` is present only for a bare numeric release.
  - Ensembl GRCh38: `https://ftp.ensembl.org/pub/release-{token}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{token}.gtf.gz`
  - Ensembl GRCh37: `https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz`
  - T2T: no derivation; rely on the token fallback in the matcher below.
- Add `find_matching_gene_annotation_release()` returning an existing
  `GeneAnnotationRelease` for this VAV's build and consortium, or None. Matching, in
  order:
  1. `gene_annotation_import__url` equal to the derived URL, comparing with the scheme
     normalised so an `http`/`https` difference does not miss.
  2. Failing that, and when the token is known, a release whose
     `gene_annotation_import.url` contains the token as a path-ish segment. This is what
     lets T2T (`2022_06`) and any future NCBI layout change still match.

  A unique match links; multiple candidates log them and link nothing.
- Replace `cdot_gene_release_filename` with the per-GFF cdot asset name derived from the
  token, so the UI can name the file that actually exists.
- Keep `suggested_gene_annotation_release_name` as the `--release=` default, now built
  from the token (e.g. `RefSeq_RS_2025_08`).

### 3. `create_new_variant_annotation_version` — auto-link, or print the command

After `get_or_create_variant_annotation_version_from_current_vep`, for each VAV with no
`gene_annotation_release`:

- Call `find_matching_gene_annotation_release()`. On a hit, set and save it, and log
  which release and why it matched.
- On a miss, log the release token and the exact command to run:

  ```
  python3 manage.py import_gene_annotation_release --genome-build=GRCh38
  ```

The command stays advisory here rather than running the download inline: creating a VAV
is quick and routine, while installing a release is a large download plus a long
gene-annotation build.

### 4. `genes/management/commands/import_cdot_gene_annotation_release.py` (new)

The one command that does everything, and the only place release creation lives.

All release handling moves here out of `import_gene_annotation`, which keeps only the
bulk-insert machinery (`import_cdot_data`, `import_cdot_data_file`, `read_cdot_version`,
`GeneAnnotationImportManager`) and loses its `--release` argument. `create_release`
becomes a module-level `create_gene_annotation_release`, and the auto-link block it used
to carry goes away — the command knows which VariantAnnotationVersion it is installing
for, so it links directly.

`--json-file` (with an optional `--release-name`) covers what
`import_gene_annotation --json-file --release=NAME` used to: creating a release from a
locally produced cdot file rather than a published one. That path skips the GitHub lookup
and download.

Arguments:

- `--genome-build` — optional, defaults to every `GenomeBuild.builds_with_annotation()`.
- `--status` — which VAV to target, defaulting to NEW then falling back to ACTIVE, matching
  how `gene_annotation --new-releases` / `--latest-releases` pick theirs.
- `--force` — reinstall even when the VAV already has a linked release.
- `--gene-annotation` / `--no-gene-annotation` — run the gene annotation step, default on.

Per build:

1. Resolve the target VAV. Skip when it already has a linked release, unless `--force`.
2. Auto-link first, via `find_matching_gene_annotation_release()`, and skip the download
   when that succeeds.
3. Derive the release token. When it cannot be derived, report the VAV's version strings
   and move on to the next build.
4. Find the cdot asset. When absent, report the token alongside the assets that *are*
   published for that build and consortium, so the mismatch is obvious.
5. Ensure transcripts are present by running the `import_cdot_latest` check for this
   build and consortium, importing the combo file when the installed cdot data version is
   behind.
6. Stream the per-GFF asset to a `NamedTemporaryFile` and import it via
   `import_gene_annotation`'s `_create_release`, using
   `suggested_gene_annotation_release_name` as the release name. Reuse the streaming
   approach already in `import_cdot_latest`.
7. Set `vav.gene_annotation_release` from the release just created, rather than inferring
   it from a URL comparison.
8. When `--gene-annotation` is on, invoke the `gene_annotation` command for the release.

### 5. `import_cdot_latest` — make the up-to-date check reusable

Extract the "is our cdot data current?" comparison into a function returning the
installed and latest versions, so both commands share it, and replace the `exit(0)` with
a return so a chained caller continues. The combo import becomes
`import_latest_combo_file`, callable as a prerequisite step.

The streaming download (`requests` → temp file → gzip handle) becomes
`download_cdot_json` in `genes/cdot_data_release.py`, shared by the combo import and the
release install rather than written out in both. On disk rather than in memory because the
importer reads the file twice, so it has to be seekable but needn't live in RAM.

### 6. `annotation/templates/annotation/annotation_build_detail.html`

Replace the stale two-line `wget` + `import_gene_annotation` instructions
(currently pointing at a download host that no longer serves these files) with the single
new command. Keep the Django Admin fallback link, since a manual override is still
occasionally needed.

### 7. `genes/gene_annotation_import_urls.py` (new)

`GeneAnnotationImport.url` is the dedup key in `GeneAnnotationImportManager`, so every
time the scheme or file format cdot recorded changed, a duplicate row appeared for the
same underlying annotation. Live counts across 117 rows:

| family | count | data on the legacy rows |
|---|---|---|
| scheme-only duplicates (`ftp://` or `http://` alongside an `https://` twin) | 32 groups | one row (18,000 gene versions) |
| `…/gff3/Homo_sapiens.GRCh38.N.gff3.gz`, each with an equivalent `…/gtf/…N.gtf.gz` row | 25 | 75,984 gene versions, no transcript versions |
| `postgresql://uta.biocommons.org/uta_…` | 2 | 19,607 transcript versions |

A single pure function, importable from both the model layer and a migration:

```python
def canonical_import_url(url: str) -> str
```

- `postgresql://` URLs (UTA imports, not GFF files) pass through untouched.
- Other schemes become `https`.
- Ensembl gff3 paths become their gtf equivalent: `/gff3/` → `/gtf/` and
  `.gff3.gz` → `.gtf.gz`.

Used by `find_matching_gene_annotation_release` so URL comparison tolerates both
variations, and by the migration below.

### 8. `genes/migrations/…_one_off_canonical_gene_annotation_import_urls.py` (new)

Collapses each set of `GeneAnnotationImport` rows sharing a `canonical_import_url` into
one row, keeping the lowest pk among those whose URL is already canonical (falling back
to the lowest pk overall). For every other row in the set:

1. Repoint `GeneVersion.import_source` and `TranscriptVersion.import_source`, and
   `GeneAnnotationRelease.gene_annotation_import`, to the keeper.
2. Delete the row.

Then set the keeper's `url` to the canonical form. The two UTA rows are untouched.
Written with `apps.get_model` and bulk `update()` calls, reversible as a no-op.

### No backfill of `VariantAnnotationVersion.gene_annotation_release`

Existing versions are already linked — the work was done by hand in Django Admin, which is
what the missing auto-link forced. Running the derivation over the live data found nothing
to fill in. New versions link themselves at creation, so a backfill migration would only
ever be a no-op.

Two GRCh37 versions do report `refseq=105.20220307` while linked to a `105.20201022`
release. Re-pointing one changes which transcript versions existing analyses resolve to,
so that stays a human decision rather than something a migration does silently.

## Testing

- `genes/tests/test_cdot_data_release.py` — asset-name parsing over the real
  `data_v0.2.33` names, including a release token containing underscores
  (`RS_2025_08`), the T2T combo-only case, and rejection of combo file names.
- `genes/tests/test_gene_annotation_import_urls.py` — canonicalisation: scheme, Ensembl
  gff3 → gtf, RefSeq archive `.gff3.gz` left alone, UTA connection strings untouched,
  idempotence.
- `annotation/tests/test_gene_annotation_release_matching.py` — token and URL derivation
  for each row in the table above, asserting the URLs equal the ones cdot actually
  records; the legacy timestamp `refseq` special case; and matching by exact URL, across
  scheme, across gff3/gtf, via the T2T token fallback, plus the ambiguous and no-match
  cases.
- `annotation/tests/test_link_gene_annotation_release.py` — `link_gene_annotation_release`
  against real model rows: links the release for our VEP rather than an earlier one for
  the same build, ignores `GeneAnnotationRelease.version`, and declines to link across
  builds, across consortia, and when candidates are ambiguous.

These use `annotation.fake_annotation.get_fake_vep_version` and run with `--keepdb`.

## Notes

- Asset resolution and the combo import both hit the network. Tests stub
  `get_latest_release_assets` rather than calling GitHub.
- T2T-CHM13v2.0 auto-links but cannot auto-install while cdot publishes no per-GFF asset
  for it. The command reports that rather than failing opaquely.
