# genes — research notes

Verified against 7c4408c62 on 2026-09-06

The genes app is the reference layer everything else resolves against: gene symbols and their aliases, the stable
RefSeq/Ensembl gene and transcript ids, the per-build, per-version transcript geometry that turns an HGVS string into a
`VariantCoordinate` and back, and the release bookkeeping that ties a VEP build to the exact gene and transcript
versions it reported. On top of that sit the things users make from symbols - gene lists (typed, PanelApp, per-sample),
canonical transcript collections for coverage, gene fusions - and the reference downloads (HGNC, MANE, LRG, gnomAD
constraint, Pfam). This document is the story behind `genes/CLAUDE.md`: how the data gets in, why symbol matching is
per release, how HGVS resolution falls through biocommons to ClinGen, and what has bitten before. Fields, URLs,
commands, tasks and signals are in the generated maps ([models](../maps/models.md#genes),
[urls](../maps/urls.md#genes), [commands](../maps/commands.md), [tasks](../maps/tasks.md#genes),
[signals](../maps/signals.md), [settings](../maps/settings.md)); `claude/domain.md` has the vocabulary.

## Flows

### Genes and transcripts arrive from cdot

Every Gene, GeneVersion, Transcript and TranscriptVersion row comes from a cdot JSON file
(https://github.com/SACGF/cdot, the spin-off that parses RefSeq/Ensembl GFFs into one JSON per build). cdot publishes
two kinds of file per data release, described in `genes/cdot_data_release.py`: a *combo* file per build and consortium
holding the latest copy of every transcript, and a *per-GFF* file per gene set VEP can be built against.
`genes/management/commands/import_cdot_latest.py:import_latest_combo_file` downloads the combo file to disk
(`genes/cdot_data_release.py:download_cdot_json` - a temp file, because the importer reads it twice and the larger
builds OOMed in memory) and `genes/management/commands/import_gene_annotation.py:Command.import_cdot_data_file`
streams genes then transcripts out of it with ijson. `genes/management/commands/import_cdot_latest.py:cdot_data_needs_update` compares the `cdot`
key stamped into the most recently inserted TranscriptVersion's `data` with the latest GitHub release tag, so a
re-run is a no-op unless cdot moved.

`genes/management/commands/import_gene_annotation.py:Command._import_cdot_data` is an upsert keyed on accession: it
loads every known symbol, gene id, transcript id and `GeneVersion.id_by_accession` /
`TranscriptVersion.id_by_accession` map up front, then bulk-creates what is new and bulk-updates the rest. RefSeq genes have no version and are stored as version 0; a cdot
gene accession starting with `_` (a fake from UTA data) is renamed with `Gene.FAKE_GENE_ID_PREFIX`. Each distinct GFF
URL in the file becomes a `genes/models/models_gene.py:GeneAnnotationImport`
(`genes/management/commands/import_gene_annotation.py:GeneAnnotationImportManager`), and a version row keeps the
`import_source` that first created it - many versions are shared between GFF releases, so a version is created once
and re-linked afterwards. The `hgnc` FK is only set when the HGNC row exists, so HGNC should be loaded first. A
transcript's whole geometry lives in `TranscriptVersion.data["genome_builds"][build]`; pre-cdot rows lack
`genome_builds`, which `genes/models/models_gene.py:TranscriptVersion.data_is_current_cdot_format` detects.

### A GeneAnnotationRelease per VEP build

VEP reports genes and transcripts from one specific GFF, and the annotation import links every row to a
TranscriptVersion by FK, so a `genes/models/models_gene_annotation_release.py:GeneAnnotationRelease` records which
gene and transcript versions that GFF contained (`ReleaseGeneVersion`, `ReleaseTranscriptVersion`). The release
stores nothing beyond those links and a `gene_annotation_import`, which is why the combo file must be in first:
`genes/management/commands/import_cdot_gene_annotation_release.py:Command._ensure_transcripts_installed` runs the
combo import if cdot has moved, then `create_gene_annotation_release` builds the link rows from the per-GFF file,
creating any GeneVersion the combo file skipped because a later transcript version pointed at a different gene.

The command is driven from the VariantAnnotationVersion, not from a filename.
`annotation/gene_annotation_release_matching.py:VepGeneSetVersions.release_token` reads the release identifier out of
the VAV's `refseq` / `vep_cache` / `genebuild` pins (`RS_2025_08`, `110`, `2022_06`),
`genes/cdot_data_release.py:find_latest_release_asset` finds the per-GFF asset with that token, and
`genes/management/commands/import_cdot_gene_annotation_release.py:Command._ensure_release` downloads, creates and
links it. Before downloading anything it tries `annotation/models/models.py:VariantAnnotationVersion.link_gene_annotation_release`,
which matches existing releases by the GFF URL in their import
(`annotation/gene_annotation_release_matching.py:match_gene_annotation_release`) - a build's gene set usually does
not change between VEP versions, so the previous VAV's release is normally the right one. An already-linked release
that disagrees with the VEP pins is reported and only re-pointed with `--relink`
(`genes/management/commands/import_cdot_gene_annotation_release.py:Command._check_existing_release`), because
changing a release changes which transcript versions existing analyses resolve to. The command finishes by running
`gene_annotation` for the release and bumping the AnnotationVersion (`Command._ensure_gene_annotation`), since a
release without a GeneAnnotationVersion makes `AnnotationVersion.validate` raise and the variant page 500.

### Symbols, aliases and per-release matching

`genes/models/models_gene.py:GeneSymbol` is keyed on the symbol string under a case-insensitive collation, with
`genes/models/models_gene.py:GeneSymbolAlias` rows from HGNC, NCBI gene_info, legacy UCSC and manual entry. Most
come from `genes/cached_web_resource/hgnc.py:save_hgnc_records`, which upserts a GeneSymbol per approved symbol and
an alias per `prev_symbol` / `alias_symbol`; `genes/cached_web_resource/hgnc.py:store_hgnc_from_web` then re-matches
every release. Resolving one string is
`genes/models/models_gene.py:GeneSymbol.cast` (exact, cached); resolving text or aliases is
`genes/gene_matching.py:GeneSymbolMatcher`, which loads the whole symbol and alias tables into upper-cased dicts on
first use and records which alias matched on the `GeneListGeneSymbol`. `genes/gene_matching.py:HGNCMatcher.match_hgnc`
is the HGNC-id variant: an approved symbol always outranks an alias redirect, because HGNC's `alias_symbol` synonyms
are frequently another gene's approved symbol (AURKAIP1 lists AIP).

A symbol does not name a gene; it names a gene *in a release*. The same symbol has pointed at different Ensembl and
RefSeq ids over time (and TAZ became TAFAZZIN between builds), so `genes/gene_matching.py:ReleaseGeneMatcher` writes
`genes/models/models_gene_annotation_release.py:ReleaseGeneSymbol` / `ReleaseGeneSymbolGene` rows per release: first
a direct hit on the release's own GeneVersion symbols, then `genes/gene_matching.py:ReleaseGeneMatcher.aliases_dict`,
which walks the alias graph (`ReleaseGeneMatcher._aliases`, loop-guarded) and the symbols other builds' GeneVersions
gave the same gene, recording the path as `match_info` so the gene list grid can show why. Readers never touch the
matcher: `genes/models/models_gene_annotation_release.py:GeneAnnotationRelease.genes_for_symbols` and
`genes/models/models_gene_list.py:GeneList.get_genes` read the cached rows, and the release always comes from the
VAV (`GeneAnnotationRelease.get_for_latest_annotation_versions_for_builds`). Matching is triggered whenever symbols
appear: `genes/gene_matching.py:GeneSymbolMatcher.create_gene_list_gene_symbols` runs
`ReleaseGeneMatcher.match_unmatched_in_hgnc_and_gene_lists` for every release after a bulk create, and the HGNC
download and release creation do the same, so every symbol in a gene list or HGNC has rows in every release.

### HGVS resolution

`genes/hgvs/hgvs_matcher.py:HGVSMatcher` is the one entry point for HGVS in both directions, built per build through
`genes/hgvs/hgvs_matcher.py:HGVSMatcher.instance` (lru-cached, since construction opens the genome fasta). It wraps a
`genes/hgvs/biocommons_hgvs/hgvs_converter_biocommons.py:BioCommonsHGVSConverter`: biocommons `hgvs` driven by
`genes/hgvs/biocommons_hgvs/data_provider.py:DjangoTranscriptDataProvider`, a cdot `LocalDataProvider` whose
transcripts are our TranscriptVersion rows and whose sequence fetcher chains the transcript sequences we hold
(`DBTranscriptSeqFetcher`, backed by `genes/models/models_gene.py:TranscriptVersionSequenceInfo.get` and, when
`HGVS_RETRIEVE_TRANSCRIPT_SEQUENCE` allows, `genes/transcript_sequence_retrieval.py:TranscriptSequenceFetcher`) ahead
of the build's fasta (`SingleBuildFastaSeqFetcher`, which refuses contigs from other builds). The ClinGen Allele
Registry is the second converter type (`genes/hgvs/hgvs_converter.py:HGVSConverterType`), enabled by
`clingen_resolution` and only for builds ClinGen supports (`genes/hgvs/hgvs_matcher.py:HGVSConverterFactory.factory`).

For a c./n. string, `genes/hgvs/hgvs_matcher.py:HGVSMatcher.get_variant_coordinate_and_details` asks
`HGVSMatcher.filter_best_transcripts_and_converter_type_by_accession` for an ordered list of (transcript version,
converter) candidates and takes the first that succeeds. The candidate versions span the lowest to highest we have
seen in any build or sequence-info row; versions we lack become a `genes/hgvs/hgvs_matcher.py:FakeTranscriptVersion`
that only ClinGen can answer, and a local candidate is added only when `TranscriptVersion.hgvs_ok`. Ordering is
`HGVSMatcher._sort_transcript_converter_types`: version distance comes from cdot's `rank_transcript_versions` (ported
out of this file, so change it there) with local-before-ClinGen layered on top. ClinGen calls are guarded by
`HGVSMatcher._clingen_allele_registry_ok` - an "unknown reference" answer is cached in Redis for a week per transcript
so the next candidate is tried without another request, and any other server error turns ClinGen off for the life of
the matcher. The result carries the transcript actually used, a `HGVSConverterInfo` (converter, method string, cdot
data version or the ClinGen call date) and two flags the UI turns into warnings:
`genes/hgvs/hgvs_converter.py:HgvsMatchRefAllele` (the reference base given differs from the genome) and
`HgvsOriginallyNormalized` (the nomen was not the normalised form). LRG accessions are first rewritten to the mapped RefSeq
transcript via `genes/models/models_gene.py:LRGRefSeqGene.get_transcript_version`
(`HGVSMatcher._lrg_get_variant_coordinate_used_transcript_method_and_matches_reference`), ClinGen taking over when
the mapping has no local transcript. g. and m. strings skip the loop; `HGVSMatcher._validate_genomic_kind` rejects an
m. on a nuclear contig, which biocommons would otherwise resolve silently as g. (#1632).

The reverse direction, `HGVSMatcher.variant_coordinate_to_hgvs_used_converter_type_and_method`, walks the same
candidate list, checks the variant's contig against the transcript's, and on the ClinGen path overwrites the gene
symbol with ours so the string is the same whichever converter produced it. g.HGVS for a SNV is string-formatted
without biocommons (`HGVSMatcher._fast_variant_coordinate_to_g_hgvs`). Symbolic DEL/DUP/INV go through a second
`AssemblyMapper` with reference replacement and normalisation off, so a 2 Mb deletion never reads 2 Mb of reference
(#1571). Parsing without the database is `genes/hgvs/hgvs.py:HGVSComponents` (with `HGVSComponents.diff` for comparing
two) and `genes/hgvs/hgvs.py:HGVSDisplay` for the view layer; the matcher is only for coordinates. A c.HGVS with a gene symbol and no transcript raises in the matcher - the search app owns that
feature: `snpdb/signals/variant_search.py:_search_hgvs` catches the error, ranks the symbol's transcripts through
`HGVSMatcher.rank_gene_symbol_transcripts` (cdot, MANE first) and searches each under the `SEARCH_HGVS_GENE_SYMBOL*`
settings (`snpdb/signals/variant_search.py:_search_hgvs_using_gene_symbol`).

### Canonical transcripts and MANE

"Canonical" is three things. cdot tags (`genes/models/models_gene.py:TranscriptVersion.tags`, read from the build's
data with `basic` stripped because Ensembl puts it on nearly everything) give `TranscriptVersion.canonical_score`
from `TranscriptVersion.CANONICAL_SCORES` - MANE Select 2, RefSeq Select 1, anything else 0 - which sorts a gene's
transcripts on the gene page and in cdot's symbol ranking. `genes/models/models_gene.py:MANE` is the MANE summary
file pinned to one version in `genes/cached_web_resource/mane.py:store_mane_from_web` (GRCh38 only, replaced
wholesale on each load), linking a symbol, HGNC and both consortia's gene and transcript versions; search uses
`MANE.get_mane_and_aliases_list_from_symbol`. VEP's own `canonical` flag is a third opinion
(`genes/models/models_gene.py:Gene.get_vep_canonical_transcript`). None of these is the
`genes/models/models_gene_coverage.py:CanonicalTranscriptCollection`, a lab's TSV of one transcript per gene
(`genes/canonical_transcripts/create_canonical_transcripts.py:create_canonical_transcript_collection` reads it;
`ChosenTranscript` overrides `CanonicalTranscript`), attached to an enrichment kit and defaulted by
`settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID` through
`genes/canonical_transcripts/canonical_transcript_manager.py:CanonicalTranscriptManager.get_canonical_collection_for_enrichment_kit`.

### Gene lists and PanelApp

A `genes/models/models_gene_list.py:GeneList` is symbols, not genes: `GeneListGeneSymbol` keeps the text the user
typed (`original_name`), the symbol it matched and the alias it matched through, so an unmatched entry is visible
rather than dropped. Typed text goes through `genes/custom_text_gene_list.py:create_custom_text_gene_list`, which
hashes the text so an unchanged edit is a no-op and otherwise replaces the list, tokenising with
`genes/gene_matching.py:tokenize_gene_symbols`; tokens over `genes/gene_matching.py:MAX_GENE_SYMBOL_LENGTH` are
skipped with a recorded warning because they cannot fit the unique index on `(gene_list, original_name)`. Analyses
read a list per release: `genes/models/models_gene_list.py:GeneList.get_q` resolves symbols to that release's genes
and delegates to `annotation/models/models.py:VariantTranscriptAnnotation.get_overlapping_genes_q`, and the gene
list node merges several lists with `GeneList.get_gene_ids_for_gene_lists`, routing symbols through a subquery so
Postgres uses the `(release, gene_symbol)` unique index instead of hash-joining every ReleaseGeneSymbol. Categories
(`genes/models/models_gene_list.py:GeneListCategory`) are a table rather than choices so forms can chain on them,
and the hidden ones mark node text, sample, coverage and PanelApp-cache lists. A sample's gene list is a
`SampleGeneList` plus an `ActiveSampleGeneList` one-to-one; `genes/models/models_gene_list.py:sample_gene_list_created`
makes the first list active and, on a second, deletes the active row so a human has to choose.

PanelApp is cached in two layers. The panel *listing* is a CachedWebResource task per server
(`genes/panel_app.py:store_panel_app_panels_from_web`, paging the API into `PanelAppPanel` rows for autocomplete);
the *genes* of a panel are fetched on demand by `genes/panel_app.py:get_panel_app_local_cache`, which reuses a
`PanelAppPanelLocalCache` for the current panel version while `genes/models/models_panel_app.py:PanelAppPanel.cache_valid`
(younger than `PANEL_APP_CACHE_DAYS`) and otherwise calls the API. Genes are identified by HGNC id rather than
PanelApp's symbol (`genes/panel_app.py:resolve_panel_app_genes`), because their symbols come from a dated Ensembl
snapshot; `genes/models/models_panel_app.py:PanelAppPanelLocalCache.get_gene_list` then builds one public GeneList
per (panel, version, minimum confidence) under `admin_bot`, naming each gene by our current approved symbol
(`PanelAppPanelLocalCacheGeneSymbol.gene_symbol_str`). A panel PanelApp no longer serves is soft-deleted rather than
removed (`genes/panel_app.py:_mark_panel_deleted`, #405), since gene lists and analyses still reference it.

### Gene coverage and reference data

`genes/models/models_gene_coverage.py:GeneCoverageCollection` is one coverage file (from SeqAuto QC or a standalone
upload) loaded by `GeneCoverageCollection.load_from_file` into two partitioned tables under the same collection:
every row (`GeneCoverage`, when `SEQAUTO_QC_GENE_COVERAGE_STORE_ALL`) and the rows whose transcript is in the kit's
canonical collection (`GeneCoverageCanonicalTranscript`, what the per-gene metrics and
`GeneCoverageCollection.get_uncovered_gene_symbols` read). Rows are `COPY`ed from CSV straight into the partition,
with one `GeneSymbolMatcher` and one `transcript_versions_by_id` map per load
(`genes/tasks/gene_coverage_tasks.py:reload_gene_coverage_collection`). The loader raises when more genes have no
canonical match than do - the symptom of the wrong collection on a kit. Collections can be archived
(`DataArchiveMixin`); `genes/models/models_gene_coverage.py:gene_coverage_collection_pre_delete_handler` drops the
partitions.

HGNC, MANE, LRG, gnomAD constraint, Pfam, UniProt, RefSeq gene summaries and the PanelApp listings all load through
`annotation.CachedWebResource` tasks (`genes/tasks/cached_web_resource_tasks.py`, wired by
`genes/signals/manual_signals.py:hgnc_post_save_handler` and siblings): saving the named resource row runs the
download on a worker and the row's description records what arrived. Pfam domains are the exception:
`genes/interpro.py:store_domains_for_transcripts` fetches them per gene from InterPro the first time a gene page
needs them (#1554). Gene fusions live here too: `genes/gene_fusions.py:GeneFusionResolver` turns `BCR::ABL1` into
`genes/models/models_gene_fusion.py:FusionGeneId` pairs that become gene-level Variants (#1506).

## Why it is shaped this way

Gene and Transcript are versionless, build-independent ids and GeneVersion / TranscriptVersion carry the build so
that a symbol, a coordinate or an HGVS can be asked about in any build without going through liftover; the price is
that the symbol lives on GeneVersion and `Gene.get_symbols` legitimately returns several. cdot replaced the in-house
GTF parsing (#566, 2022) and `TranscriptVersion.data` is the cdot record verbatim rather than columns, so the HGVS
library reads the shape it was developed against and a cdot schema change is a data reimport, not a migration. The
per-release matching cache exists because an analysis pinned to an old AnnotationVersion must keep resolving its
gene list to the genes VEP actually reported then; populating `ReleaseGeneSymbol` eagerly whenever symbols appear
keeps the read path to two indexed lookups.

HGVS moved from PyHGVS to biocommons (#839) because biocommons validates against real transcript sequence and handles
normalisation, and the local-then-ClinGen candidate list exists because neither source is complete: we lack old
transcript versions ClinGen has, and ClinGen lacks Ensembl and anything recent. `allow_alternative_transcript_version`
is a flag rather than always-on because annotation import must link the exact version VEP wrote (#1222) while search
and classification import want the nearest. Keeping `HGVSConverterInfo` on every result is what lets
`ImportedAlleleInfo` record the cdot version a classification was resolved with (#1321) and re-resolve when it moves.

## History

GeneAnnotationRelease and per-release symbol matching arrived in 2021 (#459, #494) when new gene annotation releases
started diverging from what old analyses expected. cdot replaced GTF imports in 2022 (#566), the data migrated in
place by `genes/migrations/0068_one_off_upgrade_cdot_data.py`. `GeneSymbol.symbol` moved from `citext` to a TextField
with the `case_insensitive` collation (`genes/migrations/0076_genesymbol_symbol_genesymbolalias_alias.py`), the
source of the `startswith` trap below. The HGVS layer went through hgvs_shim (#1477) and back in-tree (#1812), PyHGVS
was removed (#1678, 2026), HGVSDisplay/HGVSComponents were split out of the old CHGVS with a corpus test (#1702),
symbolic CNV HGVS stopped reading the reference (#1571) and the m. kind check was added (#1632). PanelApp gained soft
deletion (#405), spoofed headers (#1457) and HGNC-based gene identity (#1667,
`genes/migrations/0084_panel_app_local_cache_gene_symbol_hgnc.py`). GeneInfo was retired (#693), Pfam went lazy via
InterPro (#1554), gene fusions became Variants (#1506), and the release import was automated from the VAV's VEP pins
in mid-2026 (`genes/cdot_data_release.py`, `annotation/gene_annotation_release_matching.py`).

## Traps

`genes/models/models_gene.py:TranscriptVersion.get_transcript_version` defaults to `best_attempt=True` and returns
the nearest higher (or highest) version when the requested one is missing, then raises `MissingTranscript` if the
row has no valid data; `genes/models/models_gene.py:TranscriptVersion.raise_bad_or_missing_transcript` tells
`genes/transcript_errors.py:BadTranscript` (their typo) from `NoTranscript` (our gap) by asking NCBI/Ensembl, and
caches both answers for a week (`variantgrid/test_runner.py` swaps in
`genes/tests/utils/mock_transcript_sequence_retrieval.py:MockTranscriptSequenceFetcher` so tests never make them).
With `HGVS_VALIDATE_REFSEQ_TRANSCRIPT_LENGTH` on, `TranscriptVersion.hgvs_ok` excludes a RefSeq transcript whose exon
lengths do not sum to the fetched sequence length and silently routes it to ClinGen;
`TranscriptVersion.hgvs_data_errors` is where to look when a transcript "works on the other server".

`genes/models/models_gene.py:GeneSymbol.get_deterministic_queryset` exists because `symbol__startswith` raises
`NotSupportedError` on the collated column. `GeneSymbolMatcher` and `ReleaseGeneMatcher` read whole tables into
dicts on first use and `GeneSymbolMatcher.create_gene_list_gene_symbols` re-matches every release on each save, so
a loop creating lists should share one matcher. Inserting
`GeneListGeneSymbol` rows any other way leaves symbols with no release rows and the list matches nothing in analyses
until `genes/management/commands/rematch_unmatched_gene_list_symbols.py:Command` runs;
`genes/management/commands/fix_rematch_release_symbols_to_genes.py:Command` re-runs the alias walk for symbols that
have a release row but no gene.

Genes prefixed `unknown_` are legacy placeholders from pre-GFF imports; `genes/management/commands/fix_fake_genes.py:Command`
re-points their transcripts where another version names the gene and `Gene.delete_orphaned_fake_genes` removes the
rest. `LRGRefSeqGene.get_transcript_version` requires the `t` part (`LRG_199t1`); a bare `LRG_199` raises.
`genes/models/models_gene_list.py:GeneList.get_q` and
`genes/models/models_gene_annotation_release.py:GeneAnnotationRelease.get_for_latest_annotation_versions_for_builds`
import from annotation inside the method: the genes/annotation cycle is real, and new cross-imports belong at the
call site, not module level. `PanelAppPanel.cache_valid` is time-based, so a test that touches a stale panel makes a
network call; `MANE` is GRCh38-only and `MANE.get_mane_and_aliases_list_from_symbol` raises when the table is empty,
so gene-symbol HGVS search errors on a deployment that has not loaded MANE.
