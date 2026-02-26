# Ontology App Research

## Purpose

The `ontology` Django app manages medical ontologies for disease and phenotype annotation within VariantGrid. It stores ontology terms, tracks the relationships between those terms, and enables path-finding through ontology graphs to map genes to diseases. This underpins phenotype-driven analysis workflows throughout the platform.

## Ontologies Supported

### Locally Stored (fully imported and queryable)
- **MONDO** - Monarch Disease Ontology; primary disease ontology
- **OMIM** - Online Mendelian Inheritance in Man; gene-disease associations
- **HPO** - Human Phenotype Ontology; phenotype terms
- **HGNC** - HUGO Gene Nomenclature Committee; gene symbols

### External References Only (terms stored as stubs, no local relationship data)
- **DOID** - Disease Ontology
- **ORPHANET** - Rare disease ontology
- **MEDGEN** - NCBI MedGen
- **MeSH** - Medical Subject Headings

---

## Key Models

### OntologyTerm

The central model. Primary key is the full ontology identifier string (e.g., `"MONDO:0000043"`).

**Fields:**
- `ontology_service` - Enum: MONDO / OMIM / HPO / HGNC / DOID / ORPHA / MedGen / MeSH
- `index` - Numeric portion of the identifier
- `name` - Human-readable term name
- `definition` - Term definition text
- `extra` - JSONField for additional metadata
- `aliases` - ArrayField of alternative term names
- `status` - Enum: CONDITION / DEPRECATED / NON_CONDITION / STUB
- `from_import` - FK to OntologyImport

**Methods:**
- `get_from_slug()` - Looks up a term by its slug form
- `get_gene_symbol()` - Returns the gene symbol associated with this term (relevant for HGNC terms)
- `get_or_stub()` - Returns the term if it exists, otherwise creates a stub record
- `get_or_stub_cached()` - Cached version of `get_or_stub()`

**Properties:**
- `is_stub` - True if status is STUB (referenced but not yet imported)
- `is_obsolete` - True if status is DEPRECATED
- `is_leaf` - True if the term has no children in the ontology graph
- `warning_text` - Human-readable warning if the term is deprecated or a stub
- `url_safe_id` - ID safe for use in URLs (colon replaced)
- `external_url` - Link to the term on its source database website
- `best_url` - Prefers local URL if available, falls back to external_url

---

### OntologyImport

Tracks each import event, providing version history and change detection.

**Fields:**
- `import_source` - Enum: PANEL_APP_AU / MONDO / OMIM / HPO / HGNC / GENCC
- `filename` - Source file imported
- `version` - Integer, auto-incremented per (import_source, filename) pair
- `context` - Free-text context description
- `hash` - MD5 hash of the source file for change detection
- `processor_version` - Version of the import processor code used
- `processed_date` - Timestamp of processing
- `completed` - Boolean indicating successful completion

**Unique constraint:** `(import_source, filename, version)`

**Behavior on save:** Creates a PostgreSQL partition for `OntologyTermRelation` keyed on this import's ID, enabling efficient partitioned queries and bulk deletion of old import data.

---

### OntologyTermRelation

Represents a directed relationship between two ontology terms. The table is partitioned by `from_import_id` for performance.

**Fields:**
- `source_term` - FK to OntologyTerm (the "from" end)
- `dest_term` - FK to OntologyTerm (the "to" end)
- `relation` - Enum of relationship types (see below)
- `extra` - JSONField containing:
  - `strongest_classification` - Best GeneDiseaseClassification seen
  - `sources` - List of dicts with `moi`, `classification`, `submitter`
- `from_import` - FK to OntologyImport (determines partition)

**Unique constraint:** `(from_import, source_term, dest_term, relation)`

**Relation types:**
- `IS_A` - Subclass/subtype relationship
- `EXACT` - Exact mapping between terms in different ontologies
- `EXACT_SYNONYM` - Exact synonym
- `RELATED` - Related term
- `RELATED_SYNONYM` - Related synonym
- `CLOSE` - Close mapping
- `BROAD` - Broader term
- `NARROW` - Narrower term
- `ALTERNATIVE` - Alternative term
- `XREF` - Cross-reference
- `REPLACED` - Term replaced by another
- `FREQUENCY` - Frequency association (HPO phenotype frequency)
- `ASSOCIATED` - General association
- `ENTREZ_ASSOCIATION` - Gene association via Entrez
- `PANEL_APP_AU` - Association sourced from PanelApp Australia

**Methods:**
- `other_end()` - Returns the term at the opposite end of the relationship from a given starting term
- `as_mondo()` - Returns the MONDO term from either end of a cross-ontology relation
- `as_omim()` - Returns the OMIM term from either end of a cross-ontology relation
- `get_gene_disease_moi_classifications()` - Extracts mode-of-inheritance and classification data from `extra.sources`

---

### OntologyVersion

Groups a set of imports together into a coherent versioned snapshot of the ontology.

**Fields:**
- `gencc_import` - FK to OntologyImport for GenCC data
- `mondo_import` - FK to OntologyImport for MONDO data
- `hp_owl_import` - FK to OntologyImport for HPO OWL data
- `hp_phenotype_to_genes_import` - FK to OntologyImport for HPO phenotype-to-gene mappings
- `omim_import` - FK to OntologyImport for OMIM data (nullable)

**Methods:**
- `latest()` - Class method returning the most recent OntologyVersion
- `get_ontology_term_relations()` - Returns the QuerySet of OntologyTermRelation for all constituent imports
- `cached_gene_symbols_for_terms()` - Returns a cached mapping of term IDs to gene symbols
- `terms_for_gene_symbol()` - Returns all ontology terms associated with a given gene symbol
- `gene_disease_relations()` - Returns gene-to-disease relation records for this version
- `moi_and_submitters()` - Returns mode-of-inheritance and submitter data for gene-disease pairs

---

### OntologySnake

An immutable path-finding structure that traverses the ontology graph. Represents a traversal path (a "snake") from a source term through a sequence of relationships to a leaf term.

**Fields:**
- `source_term` - The starting OntologyTerm
- `leaf_term` - The current end OntologyTerm of the path
- `paths` - Ordered list of OntologyTermRelation steps taken

**Methods:**
- `snake_step()` - Extends the path by one relationship step, returning a new OntologySnake
- `all_descendants_of()` - Class method returning all descendant terms reachable from a given term
- `check_if_ancestor()` - Tests whether a given term is an ancestor of the snake's leaf
- `snake_from()` - Class method creating a snake starting from a given term
- `terms_for_gene_symbol()` - Traverses graph to find ontology terms associated with a gene
- `mondo_terms_for_gene_symbol()` - As above, filtered to MONDO terms only
- `get_children()` - Returns direct child terms of the leaf term
- `get_parents()` - Returns direct parent terms of the leaf term

**Properties:**
- `is_strong_enough` - True if the leaf relationship meets the minimum classification strength (STRONG or DEFINITIVE)
- `leaf_relationship` - The final OntologyTermRelation in the path

---

### GeneDiseaseClassification

Enum ranking the strength of gene-disease associations:

| Value | Meaning |
|-------|---------|
| REFUTED | Claim has been refuted |
| NO_KNOWN | No known disease relationship |
| ANIMAL | Evidence only in animal models |
| DISPUTED | Evidence is disputed |
| LIMITED | Limited evidence |
| SUPPORTIVE | Supportive evidence |
| MODERATE | Moderate evidence |
| STRONG | Strong evidence — meets `is_strong_enough` |
| DEFINITIVE | Definitive evidence — meets `is_strong_enough` |

STRONG and DEFINITIVE are considered strong enough for clinical use in path-finding.

---

### OntologyRelationshipQualityFilter

Preset filter levels controlling which relationships are traversed during graph path-finding:

- `NO_QUALITY` - No filtering; traverse all relationships
- `MINIMUM_QUALITY` - Only the weakest filter applied
- `MEDIUM_QUALITY` - Intermediate filtering
- `STANDARD_QUALITY` - Full clinical-quality filter (default for most workflows)

---

## Import Process

### Management Command

The `ontology_import` management command orchestrates the full import. It calls loader functions in sequence:

- `load_mondo()` - Parses MONDO OWL/OBO file
- `load_hpo()` - Parses HPO OWL and phenotype-to-gene mapping files
- `load_omim()` - Parses OMIM gene-to-phenotype data
- `load_gencc()` - Parses GenCC gene-disease classification data

### OntologyBuilder

A helper class used by each loader. Key behaviors:

- Computes MD5 hash of the source file and compares with the stored hash on the existing OntologyImport record
- Skips re-import if the hash and processor version match (no changes detected)
- Uses `bulk_create()` and `bulk_update()` for efficient database writes
- Tracks operation counts via `OperationCounter` (inserts / updates / deletes)

### Processor Versions

Each loader has a hardcoded processor version. When the processor logic changes, bumping this version forces re-import even if the source file has not changed:

- MONDO processor version: **17**
- GenCC processor version: **11**

---

## URL Patterns

| URL | Purpose |
|-----|---------|
| `/term/<slug:term>` | View a single ontology term detail page |
| `/autocomplete/HPO` | Autocomplete endpoint for HPO terms |
| `/autocomplete/OMIM` | Autocomplete endpoint for OMIM terms |
| `/autocomplete/HGNC` | Autocomplete endpoint for HGNC terms |
| `/autocomplete/MONDO` | Autocomplete endpoint for MONDO terms |
| `/api/mondo/search` | Search MONDO terms via REST API |
| `/api/ontology_term/<term>/gene_list` | Retrieve gene list for an ontology term |
| `/api/disease_relationship/<gene_symbol>` | Retrieve disease relationships for a gene symbol |
