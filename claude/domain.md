# VariantGrid domain glossary

The L0 vocabulary from [agent_system.md](plans/agent_system.md) §3: one paragraph per noun - what it is, the
model, the invariant that is easy to get wrong. Every other doc and every `vg` output uses these names.
Citations are `path:Symbol`; `scripts/vg outline <path>` shows the rest of the module.

## Genome and variants (snpdb)

**GenomeBuild / Contig** - `snpdb/models/models_genome.py:GenomeBuild`, `Contig`. GRCh37, GRCh38 and T2T-CHM13v2.0;
`GenomeBuild.builds_with_annotation()` is the set this deployment annotates (settings `ANNOTATION[build]["enabled"]`).
A Contig can belong to more than one build (MT, unplaced scaffolds), which is why a Variant has no build FK of its own.

**Locus / Sequence** - `snpdb/models/models_variant.py:Locus`, `Sequence`. A Locus is (contig, position, ref); one per
VCF line, shared by every alt at that position. Sequence rows are unique on `seq_sha256_hash`, filled in by `Sequence.save`,
so always create them through the model.

**Variant** - `snpdb/models/models_variant.py:Variant`. Build-specific: (locus, alt, svlen). Its builds come from
`locus.contig` via GenomeBuildContig (`Variant.genome_builds`), so restrict querysets to a build with
`Variant.get_contigs_q`, never by joining GenomeBuildContig. Never compare Variants across builds - go through the Allele.

**Reference variant** - a Variant whose alt is `=` (`Variant.REFERENCE_ALT`). They exist on purpose (a sample's ref call at a
locus); exclude them with `Variant.get_no_reference_q` when you mean real calls.

**Symbolic variant / SVLEN** - alts at or beyond `settings.VARIANT_SYMBOLIC_ALT_SIZE` are stored as `<DEL>` / `<DUP>` /
`<INV>` with `svlen`; `snpdb/models/models_variant.py:VariantCoordinate.as_internal_canonical_form` does the conversion, and
`Variant.qs_from_variant_coordinate` applies it for you. Gene-level events (fusions) are Variants on a fake contig - guard
coordinate code with `Variant.get_gene_level_q`.

**Variant kind** - what the grids' kind badge says a row is, read off the alt: `FUSION` / `AMP` / `LOSS` for a gene-level
alt (`library/genomics/vcf_enums.py:GeneLevelSymbolicAlt`), `DEL` / `DUP` / `INV` / `CNV` / `INS` with the size from `svlen`
for a symbolic one, and nothing at all for a small variant. Not VEP's `variant_class`
(`library/genomics/vcf_enums.py:VariantClass`), which the Effect node filters on: VEP calls a 1 Mb `<DEL>` and a 1 bp
deletion the same class, so the badge answers "is this a small variant?" where `variant_class` answers "what sort of
change?". Drawn client side by `_variantKind` in `variantgrid/static_files/default_static/js/variantgrid_formats.js`.

**VariantCoordinate** - `snpdb/models/models_variant.py:VariantCoordinate`, a pydantic value object (chrom, position, ref,
alt, svlen), the currency between HGVS, VCF and Variant. Canonicalise before lookup or insert.

**Allele** - `snpdb/models/models_variant.py:Allele`. The build-independent hub: one Allele, one Variant per build,
linked by `VariantAllele` (`Allele.variant_for_build`). Classifications, ClinGen ids and liftover hang off the Allele.
`Allele.grch37` / `grch38` are cached properties - refetch after a liftover or merge.

**VariantAllele** - `snpdb/models/models_variant.py:VariantAllele`. The Variant→Allele link per build, recording how it was
made (`origin`, `allele_linking_tool`). Filter by `genome_build` when joining or shared contigs duplicate rows.

**ClinGenAllele** - `snpdb/models/models_clingen_allele.py:ClinGenAllele`. The ClinGen Allele Registry `CA…` id and its
record; fetched over the network (`snpdb/clingen_allele.py`), bounded by `Variant.can_have_clingen_allele`. Two Alleles that
both have one cannot be merged.

**Liftover** - `snpdb/models/models_variant.py:LiftoverRun`, `AlleleLiftover`; pipelines in `snpdb/liftover.py`. Per
Allele, not per Variant: create the missing build's Variant from an Allele via ClinGen, bcftools/chain, or SAME_CONTIG.
Whole-database liftover is batched (`settings.LIFTOVER_BATCH_SIZE`) and fanned out as celery tasks - see
[operations.md#scale](guides/operations.md#scale).

## Samples and cohorts (snpdb, pedigree)

**VCF / Sample** - `snpdb/models/models_vcf.py:VCF`, `Sample`. A VCF is one imported file with an `import_status`;
a Sample is one genotype column of it. Permissions cascade from the VCF; deleting from the UI is a soft delete
(`ImportStatus.MARKED_FOR_DELETION`) finished by a celery task. Uploads arrive through the upload app (below).

**Cohort / CohortSample** - `snpdb/models/models_cohort.py:Cohort`, `CohortSample`. An ordered set of Samples. Every VCF has
an automatic cohort of all its samples; custom cohorts pick across VCFs. Membership changes go through
`Cohort.set_samples` (one version bump, one genotype rebuild), never a bulk update.

**CohortGenotypeCollection / CohortGenotype** - `snpdb/models/models_cohort.py:CohortGenotypeCollection`, `CohortGenotype`.
Genotypes are packed one row per variant per cohort, in arrays indexed by `CohortSample.cohort_genotype_packed_field_index`;
each collection is its own partition table (`library/django_utils/django_partition.py:RelatedModelsPartitionModel`).
Query with `CohortGenotypeCollection.get_annotation_kwargs` / `get_zygosity_q`.

**Trio / Duo / Quad** - `snpdb/models/models_cohort.py:Trio`, `Duo`, `Quad`. Named family structures over a Cohort
(proband, parents, sibling) with affected flags; analysis inheritance nodes read them. **Pedigree** is the general form:
`pedigree/models.py:Pedigree` over a PED file.

**VariantZygosityCountCollection** - `snpdb/models/models_zygosity_counts.py:VariantZygosityCountCollection`. Per-deployment
counts of hom/het/ref per variant, another partition table, used for the "seen in N samples" columns.

## Genes and transcripts (genes)

**GeneSymbol / Gene / GeneVersion** - `genes/models/models_gene.py:GeneSymbol`, `Gene`, `GeneVersion`. A GeneSymbol is the
HGNC-style name (with `GeneSymbolAlias`); a Gene is the stable Ensembl/RefSeq id; a GeneVersion is that gene in one
annotation release and build. Symbol→gene matching is per release (`genes/models/models_gene_annotation_release.py`).

**Transcript / TranscriptVersion** - `genes/models/models_gene.py:Transcript`, `TranscriptVersion`. Transcript is the
versionless accession; TranscriptVersion (`NM_000059.4`, per build) carries the cdot exon data used for HGVS. HGVS
resolution lives in `genes/hgvs/`, biocommons first with ClinGen fallback.

**GeneAnnotationRelease** - `genes/models/models_gene_annotation_release.py:GeneAnnotationRelease`. The GTF/cdot release a
VEP version was built against; ties `VariantAnnotationVersion` to the gene/transcript versions it reports.

**Canonical transcript** - `genes/models/models_gene_coverage.py:CanonicalTranscriptCollection`: a named per-gene choice of
representative transcript (MANE, a lab's list); `settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID` picks the default.

**GeneList** - `genes/models/models_gene_list.py:GeneList` (custom text, PanelApp cache, per-sample `SampleGeneList`); gene
list nodes and coverage use them.

## Annotation (annotation)

**AnnotationVersion** - `annotation/models/models.py:AnnotationVersion`. The bundle per build of sub-versions
(VariantAnnotationVersion, GeneAnnotationVersion, ClinVarVersion, HPA…) that an analysis is pinned to.

**VariantAnnotationVersion (VAV)** - `annotation/models/models.py:VariantAnnotationVersion`. One VEP configuration per build:
`vep` (code version), `columns_version` (which columns the pipeline writes), consortium, gene annotation release, and a
lifecycle `status` NEW → ACTIVE → HISTORICAL. "The VAV" means the ACTIVE one (`VariantAnnotationVersion.latest`). Creating a
new VAV re-annotates every variant - a shared, hours-long operation.

**AnnotationRun** - `annotation/models/models.py:AnnotationRun`. One batch (an `AnnotationRangeLock` of variant pks) through
dump → VEP → upload, with `AnnotationStatus` (`annotation/models/models_enums.py`); `vg status` counts the ones in flight.

**VariantAnnotation / VariantTranscriptAnnotation** - `annotation/models/models.py:VariantAnnotation`,
`VariantTranscriptAnnotation`. Partitioned per VAV (`SubVersionPartition`): VariantAnnotation is the representative-transcript
row the grids read; VariantTranscriptAnnotation has every transcript. `annotation_variantannotation` is the biggest table
in every deployment.

**ClinVar / ClinVarVersion** - `annotation/models/models.py:ClinVar`, `ClinVarVersion`. The imported ClinVar summary per
version; ClinVarExport (below) is the other direction.

## Classification (classification)

**Classification / ClassificationModification** - `classification/models/classification.py:Classification`,
`ClassificationModification`. A Classification is a lab's record for an allele + condition; every edit is a
ClassificationModification holding the evidence JSON at that point. Outside the owning lab only *published*
modifications are visible (`ClassificationModification.latest_for_user`), and view permission lives on the modification.

**EvidenceKey** - `classification/models/evidence_key.py:EvidenceKey`, `EvidenceKeyMap`. The schema of the evidence JSON:
type, options, share level, immutability. Name keys via `classification/enums/classification_enums.py:SpecialEKeys`.
Change evidence only through `Classification.patch_value`.

**ShareLevel / published** - `classification/enums/classification_enums.py:ShareLevel`: lab → institution → logged-in users
→ public. Publishing (`ClassificationModification.publish`) grants the level's Guardian group read access; only
`ShareLevel.is_discordant_level` levels count toward discordance.

**ImportedAlleleInfo / ResolvedVariantInfo** - `classification/models/classification_variant_info_models.py`. What the lab
sent (HGVS + transcript + build) and what it resolved to per build; the only link from a classification to an Allele. Unique
on the md5 of the imported text, so re-imports share one resolution.

**ClinicalContext / DiscordanceReport** - `classification/models/clinical_context_models.py:ClinicalContext`,
`classification/models/discordance_models.py:DiscordanceReport`. A ClinicalContext groups the published classifications for
one allele + allele-origin bucket + condition name; when their significance buckets disagree it is discordant and a
DiscordanceReport tracks the resolution. Buckets come from the EvidenceKey option metadata, not the enum.

**ClassificationGrouping** - `classification/models/classification_grouping.py:ClassificationGrouping`. The per-lab,
per-allele roll-up the listing grids read instead of scanning modifications.

**ConditionText** - `classification/models/condition_text_matching.py:ConditionText`. Free-text conditions matched to
`ontology/models/models_ontology.py:OntologyTerm` (MONDO/OMIM/HPO), per lab, so the same text resolves once.

**ClinVarExport** - `classification/models/clinvar_export_models.py:ClinVarExport`. The outbound ClinVar submission record
for a grouping, with its batch and status.

## Analysis (analysis)

**Analysis / AnalysisNode** - `analysis/models/models_analysis.py:Analysis`, `analysis/models/nodes/analysis_node.py:AnalysisNode`.
An Analysis is a DAG of nodes over one build and AnnotationVersion; each node is a *filter* contributing a Django Q, composed
down the graph into a queryset - not a stored result set. Source nodes (sample, cohort, trio, pedigree, all variants) start
the graph; filter nodes narrow it.

**NodeVersion / NodeCache / NodeTask** - `analysis/models/nodes/analysis_node.py:NodeVersion`. Each node edit bumps the
version; counts live on NodeVersion (NodeCount was folded in, #1820), caches hold materialised variant sets for slow nodes,
and NodeTask is the lease the celery scheduler works from.

**AnalysisTemplate** - `analysis/models/models_analysis.py`: a saved analysis with `AnalysisVariable`s, run per sample or
cohort (`AnalysisTemplateRun`); auto-analyses on import use them.

**VariantTag** - `analysis/models/models_variant_tag.py:VariantTag`. A user tag on a variant inside an analysis (Artifact,
Candidate…). Repeats heavily - the same artefact is re-tagged in every analysis it appears in - so aggregate in SQL
([operations.md#scale](guides/operations.md#scale)).

**Classify-queue tag** - a live `snpdb/models/models.py:Tag` with `requires_classification` set: tagging a variant with one
is asking for it to be classified, and classifying resolves the tagging rather than deleting it. A lab flags its own on the
tag settings page; `settings.TAG_REQUIRES_CLASSIFICATION` only names the one a fresh install is seeded with. A queue tag is
*for a bucket* when its `allele_origin_bucket` is that bucket or "Both", the same rule
`snpdb/models/models_enums.py:AlleleOriginFilterDefault` uses. Query with `snpdb/models/models.py:Tag.classify_queue_qs` /
`classify_queue_qs_for_bucket` rather than naming a tag.

## People, labs and permissions (snpdb, library)

**Organization / Lab** - `snpdb/models/models.py:Organization`, `Lab`. A Lab belongs to an Organization; each has a
`group_name` (`org/lab`) that names its Django Group, and lab membership *is* group membership. Classifications are owned by a
Lab; sharing levels map to these groups.

**Guardian permissions** - object-level read/write via `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin`
(`can_view`, `check_can_write`, `filter_for_user`). Standard groups `all_users` and `public`
(`library/guardian_utils.py`); `admin_bot` is the system user.

**UserSettings** - `snpdb/models/models_user_settings.py:UserSettings`. Layered preferences (Global → Organization → Lab →
User, later wins) read through `UserSettings.get_for_user`; holds the default build, columns and initial permission groups.

**VariantGridColumn / CustomColumnsCollection** - `snpdb/models/models_columns.py`. The catalogue of grid columns and a
user's chosen set; every variant grid builds its columns from them (`snpdb/grid_columns/custom_columns.py`).

**Patient** - `patients/models.py:Patient` with phenotype text matched to OntologyTerms; a Sample may link to one.

## Pipelines and operations (upload, manual, flags)

**FileUpload / UploadPipeline / UploadStep** - `upload/models/models.py`. A FileUpload is the file; an UploadPipeline is its
processing (steps, status, error) - for a VCF: preprocess → insert unknown variants → insert genotypes → wait for annotation →
finish. `UploadedVCF` and the other `Uploaded*` satellites link the pipeline to what it created.

**ManualMigrationTask** - `manual/models/manual_migration_models.py:ManualMigrationTask`. A deploy-time step a migration
registered (`manual/operations/manual_operations.py:ManualOperation`): a management command or a human action, gated by
`requires`, surfaced by `manage.py manual_outstanding` and the migrator.

**Flag** - `flags/models/models.py:Flag`. A typed, commentable flag on any model (`FlagsMixin`), used for classification
workflow states (suggestions, discordance, withdrawn) as much as for data quality.

**EnrichmentKit / SequencingRun** - `seqauto/models/models_sequencing.py:EnrichmentKit`, `seqauto/models/models_seqauto.py:SequencingRun`.
The sequencing side (SeqAuto): a run on a sequencer with an enrichment kit, whose QC and VCFs flow into VCF/Sample.
