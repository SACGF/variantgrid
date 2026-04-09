# Genes App - VariantGrid Research Document

## Overview

The `genes` Django app manages genetic reference data in VariantGrid. It handles genes, transcripts, gene lists, and related data including Panel App integration, canonical transcript management, and gene coverage data. This app is foundational to variant annotation and analysis workflows.

---

## Key Models

### HGNC

Stores HUGO Gene Nomenclature Committee data.

- **Primary Key:** Integer HGNC ID
- **Fields:**
  - `gene_symbol` (FK to GeneSymbol)
  - `alias_symbols` (text)
  - `approved_name`
  - `status`: choices are `APPROVED`, `SYMBOL_WITHDRAWN`, `ENTRY_WITHDRAWN`
  - `location` (chromosomal location)
  - External IDs: `ensembl_gene_id`, `refseq_ids`, `uniprot_ids`, `omim_ids`, `ccds_ids`, `mgd_ids`, `rgd_ids`
- **Property:** `hgnc_id` returns `"HGNC:{pk}"` format string
- **Caching:** Uses `ObjectManagerCachingRequest` for efficient repeated lookups

---

### GeneSymbol

Canonical gene symbol records.

- **Primary Key:** The symbol string itself (case-insensitive collation)
- **Key Methods:**
  - `cast(value)`: Converts a string to a `GeneSymbol` with caching; used widely for efficient symbol resolution
  - `get_genes()`: Returns associated Gene records
  - `latest_gene_version(genome_build)`: Retrieves the most recent GeneVersion for a given build
  - `get_upper_case_lookup()`: Utility for case-insensitive lookups

---

### GeneSymbolAlias

Links alternative gene names to canonical GeneSymbol records.

- **Fields:**
  - `alias` (case-insensitive field)
  - `gene_symbol` (FK to GeneSymbol)
  - `source`: choices are `NCBI`, `HGNC`, `UCSC`, `Manual`
- **Unique constraint:** `(alias, gene_symbol)`

Used for resolving historical, alternative, or consortium-specific names to their canonical form.

---

### Gene

Represents a gene independent of any genome build or version.

- **Primary Key:** Identifier string (RefSeq GeneID or Ensembl ENSG accession)
- **Fields:**
  - `annotation_consortium`: `REFSEQ` or `ENSEMBL`
- **Constant:** `FAKE_GENE_ID_PREFIX = "unknown_"` for genes that cannot be resolved
- **Key Methods:**
  - `latest_gene_version(build)`: Gets the most current GeneVersion for a build
  - `get_gene_symbol(build)`: Resolves the GeneSymbol for a given genome build
  - `get_vep_canonical_transcript()`: Returns the canonical transcript used by VEP

---

### GeneVersion

Build-specific and consortium-specific gene data.

- **Fields:**
  - `gene` (FK to Gene)
  - `version` (integer; 0 = unversioned)
  - `gene_symbol` (FK to GeneSymbol)
  - `hgnc` (FK to HGNC, nullable)
  - `description`
  - `biotype`
  - `genome_build`
  - `import_source`
- **Unique constraint:** `(gene, version, genome_build)`
- **Properties:** `accession`, `coordinate`, `chrom`, `start`, `end`, `strand`

---

### Transcript

Represents a transcript independent of version or genome build.

- **Primary Key:** Identifier string (ENST accession or NM_ accession)
- **Fields:**
  - `annotation_consortium`: `REFSEQ` or `ENSEMBL`
- **Key Methods:**
  - `latest_version(build)`: Returns the most recent TranscriptVersion for a build
  - `known_transcript_ids(build, consortium)`: Returns set of known transcript IDs for a build/consortium combination

---

### TranscriptVersion

The most important model for variant annotation. Stores versioned, build-specific transcript data in cdot JSON format.

- **Unique constraint:** `(transcript, version, genome_build)`
- **Fields:**
  - `data` (JSONField): cdot-format data. Structure: `genome_builds[build_name]` = `{contig, strand, exons: [(start, end, exon_id, cdna_start, cdna_end, gap), ...], tag, ...}`. Optional `start_codon`, `stop_codon`.
  - `biotype`
  - `contig` (FK)
- **Canonical Scores (`CANONICAL_SCORES`):**
  - `MANE_Select` = 2
  - `RefSeq_Select` = 1
  - `basic` = 0
- **Properties:**
  - `accession`: Versioned accession string
  - `is_coding`: Whether transcript encodes a protein
  - `tags`: Tags from cdot data (e.g., MANE_Select, basic)
  - `canonical_score`: Integer from `CANONICAL_SCORES` based on tags
  - `fivep_utr`, `cds`, `threep_utr`: Coordinate regions
  - `length`: Total transcript length
  - `coding_length`: CDS length
  - `alignment_gap`: Whether alignment gaps are present
  - `has_valid_data`: Data integrity check
  - `hgvs_ok`: Whether HGVS nomenclature can be computed
- **Key Methods:**
  - `get_transcript_version(build, name)`: Retrieve a specific version
  - `pyhgvs_data`: Returns data formatted for pyhgvs HGVS calculation

#### CDot Integration

`TranscriptVersion.data` stores JSON in cdot format. The exon array encodes all coordinate information: `(start, end, exon_id, cdna_start, cdna_end, gap)`. `start_codon` and `stop_codon` are optional fields.

---

### TranscriptVersionSequenceInfo

Stores transcript nucleotide sequences.

- **Key Methods:**
  - `get(accession, retrieve=True)`: Fetches sequence from DB or API if not present
  - `_get_and_store_from_refseq_api()`: Fetches from NCBI Entrez API using `Entrez.efetch`
  - `_insert_from_genbank_handle()`: Bulk insert from GenBank file handle (used for bulk imports)

Sequence retrieval falls back to external APIs (Ensembl REST for Ensembl transcripts, NCBI Entrez for RefSeq) when not cached locally.

---

### GeneAnnotationRelease

Links GTF/GFF3 imports to specific versioned releases of gene annotations.

- **Fields:**
  - `version` (e.g., `"109.20190607"`)
  - `annotation_consortium`: `REFSEQ` or `ENSEMBL`
  - `genome_build`
  - `gene_annotation_import` (FK to GeneAnnotationImport)
- **Unique constraint:** `(version, consortium, build)`
- **Key Methods:**
  - `get_genes()`: Returns all Gene records in this release
  - `transcript_versions_for_gene(gene)`: TranscriptVersions for a gene in this release
  - `transcript_versions_for_symbol(symbol)`: TranscriptVersions for a symbol in this release
  - `transcript_versions_for_transcript(transcript)`: TranscriptVersions for a transcript in this release

---

### ReleaseGeneVersion / ReleaseTranscriptVersion

Join models that track which GeneVersions and TranscriptVersions were included in each GeneAnnotationRelease. Used for auditing and reproducibility.

---

### ReleaseGeneSymbol / ReleaseGeneSymbolGene

Cache tables storing symbol-to-gene mappings per release. Enables fast lookup of which genes a symbol referred to at a specific release point in time.

---

### GeneAnnotationImport

Record of a single GTF or GFF3 file import.

- **Fields:**
  - `annotation_consortium`
  - `genome_build`
  - `url` (source URL of the import file)
  - `timestamp` (when import occurred)

---

## Gene Lists

### GeneList

User-created gene or transcript filters used in analysis nodes and sample filtering.

- **Fields:**
  - `category` (FK to GeneListCategory)
  - `name`
  - `user` (FK to User)
  - `import_status`
  - `error_message`
  - `locked` (bool)
  - `url` (optional source URL)
- **Permissions:** Uses Django Guardian for object-level permissions
- **Key Methods:**
  - `get_q(annotation_version)`: Returns a Django Q object for filtering variants by gene list
  - `get_genes(release)`: Returns genes for a given GeneAnnotationRelease
  - `clone()`: Creates a copy of the gene list
  - `can_view(user)`: Permission check
  - `can_write(user)`: Permission check
  - `get_warnings(release)`: Returns warnings about unresolved symbols
  - `filter_for_user(user)`: Queryset filtered to visible lists for a user

---

### GeneListCategory

Organizes gene lists by purpose. Known categories:
- `NodeCustomText`
- `SampleGeneList`
- `QCCoverageCustomText`
- `PathologyTest`
- `PanelAppCache`
- `GeneInfo`

---

### GeneListGeneSymbol

Links a GeneList to individual GeneSymbol records.

- **Fields:**
  - `gene_list` (FK)
  - `original_name` (the name as entered by the user)
  - `gene_symbol` (FK, nullable - null if symbol could not be resolved)
  - `gene_symbol_alias` (FK to GeneSymbolAlias, nullable)
- **Unique constraint:** `(gene_list, original_name)`

---

### CustomTextGeneList

Converts raw text input into a GeneList.

- **Fields:**
  - `sha256_hash`: Deduplicates identical text inputs
  - `name`
  - `text`: Raw gene name text
  - `gene_list` (FK to GeneList)

---

### SampleGeneList / ActiveSampleGeneList

Associates a GeneList with a Sample. A signal automatically creates an `ActiveSampleGeneList` when a `SampleGeneList` is created, tracking which gene list is currently active for a sample.

---

## Panel App Integration

### PanelAppServer

Represents a Panel App instance (e.g., Genomics England, Australasian Genomics).

---

### PanelAppPanel

Represents an individual panel within a Panel App server.

- **Fields:**
  - `server` (FK to PanelAppServer)
  - `panel_id` (integer)
  - `disease_group`, `disease_sub_group`
  - `name`
  - `status`
  - `current_version`
- **Unique constraint:** `(server, panel_id)`
- **Property:** `cache_valid`: Whether the local cache is still current

---

### PanelAppPanelLocalCache / PanelAppPanelLocalCacheGeneSymbol

Local cache of panel gene data fetched from the Panel App API.

- `get_gene_list(confidence_level)`: Constructs a filtered GeneList based on confidence:
  - `1` = Low confidence
  - `2` = Intermediate confidence
  - `3` = High confidence (most trusted)

#### Panel App Flow

1. Background Celery task fetches panel list from Panel App API
2. `store_panel_app_panels_from_web()` makes paginated API calls
3. `PanelAppPanel` records are created or updated
4. `get_panel_app_local_cache()` fetches individual panel gene data
5. Genes stored in `PanelAppPanelLocalCacheGeneSymbol`
6. `get_gene_list(confidence)` creates filtered GeneList from cache

---

## Canonical Transcripts

### CanonicalTranscriptCollection

Groups canonical transcript designations, typically associated with an enrichment kit or analysis context.

- **Fields:**
  - `description`
  - `filename`
  - `genome_build`
  - `annotation_consortium`
- **Class Method:** `get_default()`: Loads default collection from Django settings

#### Canonical Transcript Priority

Higher score wins when selecting the canonical transcript for a gene:

| Designation       | Priority Score |
|-------------------|---------------|
| MANE Select       | 2             |
| MANE Plus Clinical | 2            |
| RefSeq Select     | 1             |
| basic             | 0             |

---

## Coverage and Annotation

### GeneCoverageCollection

A `RelatedModelsPartitionModel` for storing gene-level coverage data. Partitioned for scalability across many samples and genome builds.

---

### GnomADGeneConstraint

Stores gnomAD population genetics constraint metrics per gene.

- **Fields:**
  - `gene` (FK)
  - `gnomad_oe_lof`: Observed/Expected loss-of-function ratio
  - `gnomad_oe_lof_lower`: Lower confidence interval
  - `gnomad_oe_lof_upper`: Upper confidence interval

---

### MANE

Stores MANE Select and MANE Plus Clinical transcript designations. Used for canonical transcript scoring.

---

### UniProt

Stores protein function and pathway data linked to genes/transcripts.

---

### Pfam / PfamSequence / PfamDomains

Protein domain annotations from the Pfam database. `PfamDomains` links specific domain regions to transcript sequences stored in `PfamSequence`.

---

### LRGRefSeqGene

Maps LRG (Locus Reference Genomic) identifiers to RefSeq transcript accessions. Important for clinical reporting contexts where LRG identifiers are required.

---

## Documentation and Metadata

### GeneInfo

Special tags attached to specific genes with:
- `icon_css_class`: CSS class for display icon
- `description`: Free-text description
- Linked `GeneList`: Genes associated with this info tag

Used for highlighting genes of special clinical interest.

---

### GeneSymbolWiki / GeneListWiki

Wiki-style documentation that can be attached to GeneSymbol or GeneList records. Supports versioned, user-editable content.

---

## Management Commands

| Command | Purpose |
|---|---|
| `import_cdot_latest` | Imports latest cdot transcript data |
| `import_gene_annotation` | Imports GTF/GFF3 gene annotation files |
| `import_canonical_transcript` | Imports canonical transcript designations |
| `import_hgnc` | Imports HGNC gene nomenclature data |
| `rematch_unmatched_gene_list_symbols` | Retries matching gene list symbols that previously failed |
| `fix_fake_genes` | Resolves genes with FAKE_GENE_ID_PREFIX |
| `fix_hgnc` | Fixes inconsistencies in HGNC data |

---

## URL Patterns

### View URLs

| URL Pattern | Purpose |
|---|---|
| `view_gene/<id>` | Gene detail page |
| `view_gene_symbol/<symbol>` | Gene symbol detail page |
| `view_transcript/<id>` | Transcript detail page |
| `gene_lists` | Gene list index |
| `view_gene_list/<id>` | Gene list detail page |
| `view_canonical_transcript_collection/<pk>` | Canonical transcript collection detail |

### Autocomplete Endpoints

- Gene
- GeneSymbol
- Transcript
- PanelAppPanel
- GeneList

### REST API Endpoints

- Gene list CRUD
- Panel App data
- `text_to_gene_list`: Converts free text to a GeneList
- Gene annotation release data

---

## Data Flow Summary

1. **Import:** GTF/GFF3 files imported via `import_gene_annotation` → creates `GeneAnnotationImport`, `GeneAnnotationRelease`, `Gene`, `GeneVersion`, `Transcript`, `TranscriptVersion` records
2. **cdot:** `TranscriptVersion.data` stores coordinate and exon data in cdot JSON format for HGVS calculation
3. **Symbol resolution:** `GeneSymbol.cast()` and `GeneSymbolAlias` handle alias resolution across consortia
4. **Gene lists:** Users create lists via text input → `CustomTextGeneList` → `GeneListGeneSymbol` with symbol matching
5. **Panel App:** API fetch → local cache → confidence-filtered `GeneList`
6. **Canonical transcripts:** Scored by MANE/RefSeq designations; used by VEP and variant annotation
