# genes — agent notes
Owns: Gene/GeneVersion, GeneSymbol/GeneSymbolAlias, HGNC, Transcript/TranscriptVersion (cdot data), GeneAnnotationRelease and per-release symbol→gene matching, HGVS parsing and resolution (biocommons + ClinGen fallback), LRG, MANE, canonical transcripts, gene lists (custom text, sample, PanelApp cache), gene coverage, gene fusions, gnomAD constraint, Pfam/UniProt.
Start with:
- models/models_gene.py — GeneSymbol, Gene, GeneVersion, Transcript, TranscriptVersion, TranscriptVersionSequenceInfo, LRGRefSeqGene, HGNC, MANE. models/__init__.py star-imports every models_*.py, so `from genes.models import X` still works.
- models/models_gene_annotation_release.py — GeneAnnotationRelease, ReleaseGeneVersion/ReleaseTranscriptVersion, ReleaseGeneSymbol/ReleaseGeneSymbolGene
- models/models_gene_list.py — GeneList, GeneListGeneSymbol, GeneListCategory, CustomTextGeneList, SampleGeneList/ActiveSampleGeneList
- hgvs/hgvs_matcher.py — HGVSMatcher (string → VariantCoordinate and back); hgvs/hgvs.py — HGVSComponents (lenient string parsing, no DB)
- gene_matching.py — GeneSymbolMatcher (text → GeneListGeneSymbol), ReleaseGeneMatcher (symbol → genes per release), HGNCMatcher
- management/commands/import_cdot_latest.py wraps management/commands/import_gene_annotation.py (Gene/Transcript versions from cdot JSON); management/commands/import_cdot_gene_annotation_release.py then creates the GeneAnnotationRelease a VEP build uses.
Patterns here:
- Gene and Transcript are stable build-independent ids; GeneVersion and TranscriptVersion carry the build and are unique on (id, version, genome_build) (models/models_gene.py:GeneVersion, models/models_gene.py:TranscriptVersion). Always filter versions by genome_build.
- The symbol lives on GeneVersion, not Gene: use `Gene.get_gene_symbol(genome_build)` and expect `Gene.get_symbols()` to return several over time (models/models_gene.py:Gene).
- Resolve a string to a GeneSymbol with models/models_gene.py:GeneSymbol.cast (cached); resolve aliases and bulk text with gene_matching.py:GeneSymbolMatcher, which also writes the per-release ReleaseGeneSymbolGene rows.
- Symbol→gene lookups are per release: use `GeneAnnotationRelease.genes_for_symbol()` / `GeneList.get_genes(release)` (models/models_gene_annotation_release.py:GeneAnnotationRelease, models/models_gene_list.py:GeneList.get_genes). Take the release from `VariantAnnotationVersion.gene_annotation_release`, never "latest".
- Transcript geometry and tags come from cdot JSON in `TranscriptVersion.data["genome_builds"][build]`; read tags via models/models_gene.py:TranscriptVersion.tags and canonical-ness via models/models_gene.py:TranscriptVersion.CANONICAL_SCORES (MANE Select 2, RefSeq Select 1; "basic" is stripped, not scored).
- Get an HGVSMatcher with hgvs/hgvs_matcher.py:HGVSMatcher.instance (lru-cached per build; construction opens the genome fasta). Pass `clingen_resolution=False` in batch/offline code.
- hgvs/hgvs_matcher.py:HGVSMatcher.get_variant_coordinate_and_details tries local biocommons first, then the ClinGen Allele Registry, across neighbouring transcript versions (hgvs/hgvs_matcher.py:HGVSMatcher.filter_best_transcripts_and_converter_type_by_accession); version-distance ranking is delegated to the external cdot package, so change ranking there.
- Symbol-only HGVS ("BRCA1:c.100A>G") is a search feature, not a matcher feature: snpdb/signals/variant_search.py ranks transcripts with hgvs/hgvs_matcher.py:HGVSMatcher.rank_gene_symbol_transcripts under the `SEARCH_HGVS_GENE_SYMBOL*` settings; the matcher itself raises on a transcript-less c.HGVS.
- Parse or rewrite an HGVS string without touching the DB using hgvs/hgvs.py:HGVSComponents; use HGVSMatcher only when you need coordinates.
- Consortium is the single-letter `AnnotationConsortium` (R/E) on Gene, Transcript, GeneAnnotationRelease and CanonicalTranscriptCollection (models_enums.py:AnnotationConsortium); infer it from an accession with `AnnotationConsortium.get_from_transcript_accession()`.
- Reference data (HGNC, MANE, LRG, Pfam, gnomAD constraint, PanelApp panel lists) loads through CachedWebResource tasks (tasks/cached_web_resource_tasks.py:HGNCWebResourceTask, cached_web_resource/hgnc.py:store_hgnc_from_web), not management commands.
Gotchas:
- "Canonical" means three different things: VEP's per-annotation flag (models/models_gene.py:Gene.get_vep_canonical_transcript), cdot MANE/RefSeq Select tags (TranscriptVersion.canonical_score), and CanonicalTranscriptCollection, a per-enrichment-kit list used for coverage whose default is `settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID` (canonical_transcripts/canonical_transcript_manager.py:CanonicalTranscriptManager). Name which one you mean.
- `GeneSymbol.symbol` is a case-insensitive collation PK: `symbol__startswith` raises NotSupportedError; use models/models_gene.py:GeneSymbol.get_deterministic_queryset and filter on `symbol_deterministic`.
- models/models_gene.py:TranscriptVersion.get_transcript_version with `best_attempt=True` (the default) silently returns the nearest higher (or highest) version when the requested one is missing; pass `best_attempt=False` when the exact version matters.
- A TranscriptVersion with `data == {}` or missing build data is legacy: check models/models_gene.py:TranscriptVersion.hgvs_ok before HGVS work — the matcher skips those and falls back to ClinGen via a hgvs/hgvs_matcher.py:FakeTranscriptVersion.
- LRG c.HGVS is rewritten to the mapped RefSeq accession via models/models_gene.py:LRGRefSeqGene before resolution; models/models_gene.py:TranscriptVersion.get_for_lrg only knows LRG_199t1 and is not the general path.
- Genes whose id starts with `Gene.FAKE_GENE_ID_PREFIX` ("unknown_") are legacy placeholders (models/models_gene.py:Gene); hide them from users and expect them to lack versions.
- A GeneList's symbols only reach analyses once matched into every release (gene_matching.py:ReleaseGeneMatcher); after inserting GeneListGeneSymbol rows outside GeneSymbolMatcher, run `match_unmatched_in_hgnc_and_gene_lists()` or `manage.py rematch_unmatched_gene_list_symbols`.
- models/models_gene_list.py:GeneList.get_q imports from annotation inside the method — the genes/annotation import cycle is real; keep new cross-imports out of module level.
- Creating a second SampleGeneList for a sample deletes the ActiveSampleGeneList instead of switching it (models/models_gene_list.py:sample_gene_list_created); set the active one explicitly.
- `PanelAppPanel.cache_valid` expires after `settings.PANEL_APP_CACHE_DAYS` (models/models_panel_app.py:PanelAppPanel.cache_valid); panel_app.py:get_panel_app_local_cache re-fetches from the live API when stale, so tests must not depend on it.
- GeneCoverageCollection is a partitioned model (models/models_gene_coverage.py:GeneCoverageCollection); delete via the model so partitions are dropped.
- gene_matching.py:GeneSymbolMatcher and gene_matching.py:ReleaseGeneMatcher cache whole-table dicts on first use; build one per import, not per symbol.
Tests:
- annotation/tests/test_data_fake_genes.py:create_fake_transcript_version builds Gene/GeneVersion/Transcript/TranscriptVersion (RUNX1, ENST00000300305.7) for a build; `create_gata2_transcript_version` / `create_pten_transcript_version` add RefSeq examples.
- Pair those with annotation/fake_annotation.py:get_fake_annotation_version, which creates the GeneAnnotationRelease and VariantAnnotationVersion that release-scoped code needs.
- Transcript sequence fetches are mocked for the whole suite by variantgrid/test_runner.py setting `TranscriptSequenceFetcher.override_class` to tests/utils/mock_transcript_sequence_retrieval.py:MockTranscriptSequenceFetcher — add new accessions to the fasta files in tests/test_data rather than hitting NCBI/Ensembl.
- tests/test_urls.py:Test is the URLTestCase for every genes page (owner vs non-owner permission checks included).
- tests/test_hgvs_corpus.py:HGVSCorpusTests runs tests/test_data/hgvs_corpus.tsv (400+ strings, malformed on purpose) through HGVSComponents; add new edge cases there.
- tests/test_hgvs.py needs the fake annotation version and a genome fasta for coordinate resolution, and is the slow one.
Deep reference: claude/research/genes.md · claude/maps/models.md#genes
