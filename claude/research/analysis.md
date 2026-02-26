# VariantGrid: `analysis` Django App

## Overview

The `analysis` app provides an interactive, graph-based variant filtering pipeline. Users construct DAG (directed acyclic graph) workflows where nodes filter and transform variants, enabling complex multi-step analysis of genomic data. Supports single-sample, trio, cohort, and pedigree analysis modes. Results are displayed in JQGrid datatables with dynamic column configurations. The system is built on `django-dag` for graph structure and Celery for asynchronous execution of node computations.

---

## 1. PURPOSE

The `analysis` app is the core interactive analysis engine of VariantGrid. It allows users to:

- Build visual graph workflows where each node represents a filtering or transformation step applied to genomic variants.
- Chain nodes together to form complex multi-step pipelines with branching and merging logic.
- Analyse single samples, trios (mother/father/proband), cohorts (groups of samples), and pedigrees.
- Save and reuse analysis workflows as **templates** that can be parameterised with new inputs.
- Display filtered variant results in paginated, sortable, filterable JQGrid datatables with dynamically added columns based on node type.
- Cache intermediate results at multiple levels (Redis, database VariantCollections) for performance.
- Automatically trigger analyses on new VCF imports via `AutoLaunchAnalysisTemplate`.

Architecturally, the app is built around a Q-object pipeline system where each node produces Django Q objects representing its filtering logic. These Q objects are composed and applied to database querysets in a staged manner to produce filtered variant result sets.

---

## 2. CORE MODELS

### `Analysis`
Container for an analysis workflow. Defined in `analysis/models/models.py`.

**Key Fields:**
- `name` (str): Human-readable name for the analysis.
- `description` (str): Optional extended description.
- `genome_build` (FK to `GenomeBuild`): The genome build (GRCh37/GRCh38) the analysis operates on.
- `user` (FK to `User`): Owner of the analysis.
- `custom_columns_collection` (FK): Column configuration for the result grid.
- `annotation_version` (FK): Pinned annotation version used for this analysis.
- `version` (int): Incremented to invalidate all node caches when analysis-level settings change.
- `lock_input_sources` (bool): Prevents modification of source node inputs.
- `visible` (bool): Controls whether the analysis appears in listings.
- `template_type` (enum): One of `None` (normal analysis), `TEMPLATE` (reusable template), or `SNAPSHOT` (immutable snapshot of a template version).

**Permissions:** Uses Guardian object-level permissions. Audit-logged via `django-auditlog`.

**Methods:**
- Methods for cloning the analysis (used by the template system).
- Methods for retrieving the node graph structure.

---

### `AnalysisLock`
Tracks the history of lock and unlock actions on an analysis.

**Fields:**
- `analysis` (FK to `Analysis`)
- `locked` (bool): Whether the action was a lock (`True`) or unlock (`False`).
- `date` (datetime): Timestamp of the action.

---

### `AnalysisNodeCountConfiguration`
OneToOne with `Analysis`. Defines which count types are displayed in the node count badges on the graph UI.

**Count types include:** TOTAL, CLINVAR, OMIM, and other domain-specific category counts.

---

### `AnalysisVariable`
Extracts specific node fields as templatable variables, enabling analysis templates to expose fields that callers must populate when instantiating the template.

**Fields:**
- `node` (FK to `AnalysisNode`): The node containing the field.
- `field` (str): The field name on the node (e.g. `sample`).
- `class_name` (str): The fully qualified class name of the field's type (e.g. `snpdb.Sample`), used for type-safe population of values.

**Example:** A `SampleNode` with its `sample` field marked as an `AnalysisVariable` means any template run must provide a `Sample` object for that slot.

---

### `AnalysisTemplate`
A locked-down, reusable analysis workflow.

**Fields:**
- `name` (unique str): Identifier for the template.
- `analysis` (OneToOne FK to `Analysis`): The live template analysis (with `template_type=TEMPLATE`).
- `user` (FK): Owner/author.
- `deleted` (bool): Soft-delete flag.
- `active` (bool): Whether the template is currently available for use.

**Key Methods:**
- `new_version()`: Creates an immutable snapshot of the current template state and registers it as a new `AnalysisTemplateVersion`. The snapshot has `template_type=SNAPSHOT`.
- `requires_sample_somatic` (property): Whether any source node requires a somatic sample.
- `requires_sample_gene_list` (property): Whether any source node requires a gene list.

---

### `AnalysisTemplateVersion`
Represents a specific versioned snapshot of an `AnalysisTemplate`.

**Fields:**
- `version` (int): Monotonically increasing version number within the template.
- `analysis_snapshot` (OneToOne FK to `Analysis`): The immutable snapshot analysis.
- `active` (bool): Only one version per template can be active at a time.
- `analysis_name_template` (str): A Python string template used to name analyses created from this version. Example: `"%(template)s for %(input)s"`.
- `appears_in_autocomplete` (bool): Controls whether this template version shows up in autocomplete suggestions.
- `appears_in_links` (bool): Controls whether this template version appears in quick-launch link lists.

**Note:** Uses `PROTECT` foreign key constraint to prevent deletion of snapshots that are referenced by runs.

---

### `AnalysisTemplateRun`
Created when a user instantiates a template to run it against specific inputs.

**Fields:**
- `template_version` (FK to `AnalysisTemplateVersion`): The template version being run.
- `analysis` (FK to `Analysis`): The newly cloned analysis created for this run.

**Factory:** `AnalysisTemplateRun.create()` clones the snapshot analysis into a new editable analysis, then calls `populate_arguments()` to validate and set all `AnalysisVariable` values.

---

### `AnalysisTemplateRunArgument`
Stores the populated values for each `AnalysisVariable` in a template run.

**Fields:**
- `variable` (FK to `AnalysisVariable`): The variable being populated.
- `object_pk` (str): Primary key of the object being assigned (serialised as string for generic FK support).
- `value` (str): String representation of the assigned value.
- `error` (str): Error message if population failed validation.

---

### `AnalysisNode` (Base Class)
The central abstract base class for all node types. Uses `django-dag` for graph structure. All concrete node types extend this class.

**Positional/Display Fields:**
- `analysis` (FK to `Analysis`)
- `name` (str): Display name of the node.
- `x`, `y` (int): Canvas position coordinates for the graph editor UI.
- `visible` (bool): Whether the node is shown in the graph.
- `output_node` (bool): Marks nodes whose results are displayed in the main result grid.
- `auto_node_name` (bool): Whether the name is auto-generated from node configuration.
- `appearance_version` (int): Incremented when display-only properties change (avoids full cache invalidation).
- `shadow_color` (str): Visual styling hint for the node on the canvas.

**Execution State Fields:**
- `version` (int): Incremented when node configuration changes (cache invalidation key).
- `ready` (bool): Whether the node has been computed and results are available.
- `valid` (bool): Whether the node configuration is currently valid.
- `count` (int, nullable): The number of variants passing through this node.
- `errors` (JSONField): Validation/execution errors.
- `load_seconds` (float): Last recorded execution time.
- `cloned_from` (FK to `NodeVersion`): Tracks provenance when nodes are cloned from templates.
- `status` (`NodeStatus` enum): Current execution state. Values:
  - `DIRTY`: Configuration changed, needs recomputation.
  - `QUEUED`: Celery task submitted, waiting to execute.
  - `LOADING_CACHE`: Pre-caching variant collection.
  - `LOADING`: Actively computing.
  - `READY`: Computation complete, results available.
  - `ERROR_CONFIGURATION`: Node has invalid configuration.
  - `ERROR_WITH_PARENT`: A parent node has an error blocking this node.
  - `ERROR_TECHNICAL`: Unexpected technical error during computation.
  - `CANCELLED`: Task was cancelled (e.g. superseded by newer version).

**Graph Configuration Fields:**
- `min_inputs` (int): Minimum number of parent nodes required (0 for source nodes, 1 for filter nodes).
- `max_inputs` (int): Maximum number of parent nodes allowed (`PARENT_CAP_NOT_SET` for unlimited).
- `uses_parent_queryset` (bool): Whether this node builds on the parent's queryset.
- `disabled` (bool): Temporarily disables the node without removing it from the graph.
- `queryset_requires_distinct` (bool): Whether the resulting queryset must apply `.distinct()`.

**Key Methods:**
- `_get_node_q()`: **Primary override point for subclasses.** Returns the Q object representing this node's filtering logic (without parent filters). Must be implemented by all concrete node types.
- `_get_node_arg_q_dict()`: Returns an `arg_q_dict` (see Q Object System section) for this node's own filters.
- `get_arg_q_dict(disable_cache=False)`: Composes parent and own arg_q_dicts. Caches the result in Redis using the node's version as part of the cache key. Returns `dict[Optional[str], dict[str, Q]]`.
- `get_queryset()`: Applies annotation kwargs and Q filters from the composed arg_q_dict in stages, grouped to prevent double-joins. Returns a Django queryset of variants.
- `_get_model_queryset()`: Returns the base queryset (typically `Variant.objects.filter(...)`) before node-specific filters are applied.
- `_get_node_contigs()`: Returns the set of contigs (chromosomes) relevant to this node, used for optimisation.
- `get_parent_subclasses()`: Returns the concrete subclass instances of all parent nodes.
- `bump_version()`: Increments `version`, which deletes the current `NodeVersion` and cascades to delete `NodeCache` and `NodeCount` entries, forcing full recomputation.
- `load()`: Triggers the node's Celery task chain.
- `node_counts()`: Returns variant count breakdowns by configured label types.

---

### `AnalysisEdge`
The concrete edge model for `django-dag`. Represents a directed connection between two nodes.

**Fields:**
- `parent` (FK to `AnalysisNode`): The upstream/source node.
- `child` (FK to `AnalysisNode`): The downstream/consuming node.

Audit-logged via `django-auditlog`.

---

### `NodeTask`
Tracks active Celery tasks for node computation.

**Fields:**
- `node` (FK to `AnalysisNode`)
- `version` (int): The node version this task is computing.
- `analysis_update_uuid` (UUID): Groups tasks belonging to the same analysis update run.
- `celery_task` (str): Celery task ID for revocation/tracking.
- `db_pid` (int): Database process ID for `pg_cancel_backend()` calls.

**Unique constraint:** `(node, version)` - only one task per node version.

---

### `NodeVersion`
A version snapshot object that ties together a node and its version number. Serves as the primary key for cache entries.

**Lifecycle:** Deleted when `node.bump_version()` is called, which cascades to delete all `NodeCache` and `NodeCount` entries for that version, ensuring cache invalidation.

---

### `NodeCache`
OneToOne with `NodeVersion`. Links to a `VariantCollection` that has been pre-computed and stored in the database for expensive nodes.

**Fields:**
- `node_version` (OneToOne FK to `NodeVersion`)
- `variant_collection` (FK to `VariantCollection`): The pre-computed set of variant IDs.
- Processing status tracking fields.

Only used for nodes with `use_cache=True` (currently `VennNode` and `IntersectionNode` when using BED file collections).

---

### `NodeCount`
Stores variant count values grouped by label for a specific node version.

**Fields:**
- `node_version` (FK to `NodeVersion`)
- `label` (str): Count category (e.g. `"total"`, `"clinvar"`, `"omim"`).
- `count` (int): The count value.

**Unique constraint:** `(node_version, label)`.

---

### `NodeColumnSummaryCacheCollection`
Caches value distribution summaries per column for a node. Used to power the column summary panel in the UI (showing top values for a column across variants in a node).

---

### `NodeVCFFilter`
Associates VCF FILTER tags with a node for VCF FILTER-based filtering.

**Fields:**
- `node` (FK to `AnalysisNode`)
- `vcf_filter` (FK to `VCFFilter`, nullable): A specific FILTER tag. `None` means the PASS filter (only variants with PASS or empty FILTER field).

---

### `NodeAlleleFrequencyFilter`
OneToOne with `AnalysisNode`. Defines allele frequency filtering logic that can be applied to a node's variants.

**Fields:**
- `group_operation` (enum): `ANY` (OR) or `ALL` (AND) logic when multiple AF ranges are configured.

**Key Methods:**
- `get_q()`: Returns the composite Q object for all configured AF ranges.
- `get_sample_arg_q_dict()`: Returns the arg_q_dict for sample-level AF filters, using the cohort genotype annotation alias to avoid double-joins.

---

### `NodeAlleleFrequencyRange`
An individual allele frequency range constraint within a `NodeAlleleFrequencyFilter`.

**Fields:**
- `node_allele_frequency_filter` (FK to `NodeAlleleFrequencyFilter`)
- `min_af` (float, nullable): Minimum allele frequency (inclusive).
- `max_af` (float, nullable): Maximum allele frequency (exclusive).
- Population source (which AF column to apply the range to).

---

## 3. Q OBJECT SYSTEM (Core Filtering Mechanism)

### `arg_q_dict` Structure

The Q object system is the heart of the filtering pipeline. Each node produces an `arg_q_dict` of the form:

```python
dict[Optional[str], dict[str, Q]]
```

The outer key is an **annotation alias** (or `None`):
- `None`: Filters that apply directly to the base `Variant` queryset with no extra annotation. Always applied.
- `"cohort_genotype_alias"` (e.g. `"cga_12345"`): Filters that require the `CohortGenotypeCollection` annotation. Grouped to a single annotation call to avoid duplicate expensive JOIN operations.
- `"variant_transcript_annotation"`: Filters that require joining to the transcript annotation table.

The inner dict maps **argument names** (strings, typically the field path) to **Q objects** representing the filter condition.

**Example:**

```python
{
    None: {
        "variant__locus__contig__genomebuild": Q(variant__locus__contig__genomebuild=build),
    },
    "cohort_genotype_alias": {
        "zygosity_het": Q(cohort_genotype_alias__samples_zygosity__contains=HET_MASK),
    },
}
```

### Composition

`get_arg_q_dict()` composes the parent node's `arg_q_dict` with this node's own `_get_node_arg_q_dict()` result. For nodes with `uses_parent_queryset=True`, Q objects are merged by annotation key.

### Caching

`get_arg_q_dict(disable_cache=False)` caches the result in Redis:
- **Cache key format:** `{node_version_id}:q_cache={disable_cache}`
- The version-based key means cache is automatically invalidated when `bump_version()` creates a new `NodeVersion`.

### Queryset Application

`get_queryset()` takes the composed `arg_q_dict` and:
1. Groups Q objects by annotation alias.
2. Calls `queryset.annotate(**annotation_kwargs)` for each alias group.
3. Applies the Q filters for that group via `queryset.filter(Q(...))`.
4. Handles `queryset_requires_distinct` to add `.distinct()` when needed.
5. Returns the final queryset of `Variant` objects (or related model rows).

This staged approach ensures each expensive JOIN annotation is only added once, regardless of how many filter conditions use that annotation.

---

## 4. ALL NODE TYPES

### Source Nodes (`min_inputs=0`)

These nodes produce variant sets from data sources without requiring parent input.

#### `SampleNode`
Loads variants from a single `Sample` object. The most fundamental source node.

**Key Fields:**
- `sample` (FK to `Sample`): The sample to load variants from.
- `sample_gene_list` (FK): Optional gene list to restrict variants.
- `restrict_to_qc_gene_list` (bool): Whether to restrict to QC-defined gene list.
- `min_ad` (int): Minimum allelic depth.
- `min_dp` (int): Minimum read depth.
- `min_gq` (int): Minimum genotype quality.
- `max_pl` (int): Maximum phred-scaled likelihood.
- `zygosity_ref` / `zygosity_het` / `zygosity_hom` / `zygosity_unk` (bool): Which zygosity classes to include.

**Filtering:** Builds the `arg_q_dict` with:
- Zygosity filter Q objects under the `cohort_genotype_alias` key (avoids re-joining the `CohortGenotype` table).
- Genotype quality (`gq`, `dp`, `ad`) filter Q objects.
- `NodeAlleleFrequencyFilter` Q objects if configured.

#### `TrioNode`
Filters variants by Mendelian inheritance patterns within a trio (proband + mother + father).

**Key Fields:**
- `trio` (FK to `Trio`): The trio being analysed.
- `inheritance` (`TrioInheritance` enum):
  - `RECESSIVE`: Both parents het, proband hom-alt.
  - `COMPOUND_HET`: Two different het variants in the same gene, one from each parent.
  - `DOMINANT`: At least one parent affected, proband het or hom-alt.
  - `DENOVO`: Neither parent has the variant, proband does.
  - `XLINKED_RECESSIVE`: Recessive on the X chromosome.
- `require_zygosity` (bool): Strictly enforce zygosity expectations.

**Compound Het Complexity:** The `COMPOUND_HET` mode is the most computationally expensive. It identifies genes where the proband has two or more het variants and each parent contributes a different one (i.e. one het from mum, one het from dad). This requires a self-join or subquery on the gene annotation to find variants that co-occur in the same gene.

#### `CohortNode`
Filters variants from a `Cohort` (a named group of samples).

**Key Fields:**
- `cohort` (FK to `Cohort`)
- `zygosity` (`SimpleZygosity`): HET, HOM_ALT, or both.
- `zygosity_op` (enum): `ALL` (all samples must match) or `ANY` (at least one sample must match).
- `accordion_panel` (enum): `COUNT` (filter by het/hom counts), `SIMPLE_ZYGOSITY`, or `PER_SAMPLE_ZYGOSITY` (individual sample zygosity filters).
- Per-sample filter fields for the `PER_SAMPLE_ZYGOSITY` mode.

**Dynamic Columns:** Adds `ref_count`, `het_count`, `hom_count` annotation columns to the result grid, showing zygosity distribution across cohort samples.

#### `PedigreeNode`
Similar to `TrioNode` but operates on a `Pedigree` object that can contain any number of related individuals. Supports the same inheritance modes but with more flexible family structures.

#### `ClassificationsNode`
Loads variants that have associated `Classification` objects (clinical variant classifications). Queries `Classification` objects filtered by criteria such as clinical significance and links back to variants.

#### `AllVariantsNode`
Returns all variants in the genome build without any source-level filtering. Typically used as a starting point when the goal is to filter down from the entire variant catalogue.

---

### Filter Nodes (`min_inputs=1`)

These nodes take the variant set from their parent node(s) and apply additional filtering.

#### `FilterNode`
Generic column-based filtering using JQGrid's filter JSON format.

**Mechanism:** The JQGrid UI sends filter rules as JSON (column name, operator, value). `FilterNode` deserialises this JSON and converts it to Django Q objects, supporting operators like `eq`, `ne`, `lt`, `gt`, `bw` (begins with), `cn` (contains), `in`, etc.

This is the most flexible filter node, allowing users to filter on any column available in the analysis grid.

#### `ZygosityNode`
Filters variants by zygosity for a specific sample.

**Key Fields:**
- `sample` (FK to `Sample`)
- `zygosity` (enum): `HET`, `HOM_ALT`, `MULTIPLE_HIT`, etc.
- `exclude` (bool): Inverts the filter (exclude variants with this zygosity).

**Multiple Hit Mode:** `MULTIPLE_HIT` queries `gene_counts` to find genes where the sample has two or more hits (relevant for compound het analysis without the full trio structure).

#### `TagNode`
Filters variants by variant tags applied by users.

**Key Fields:**
- `parent_input` (bool): If True, only considers tags applied within the current analysis context; if False, considers tags globally.
- `exclude` (bool): Inverts the filter (excludes tagged variants).
- `mode` (enum): `THIS_ANALYSIS` (tags in this analysis only) or `ALL_TAGS` (tags in any analysis).

**Implementation:** Queries `VariantTag` objects and builds Q objects on variant PKs.

#### `PhenotypeNode`
Filters variants to those in genes associated with phenotypes or diseases.

**Key Fields:**
- `text_phenotype` (str): Free-text phenotype search.
- `patient` (FK to `Patient`): If set, uses the patient's recorded phenotypes.
- `accordion_panel` (enum): `ONTOLOGY` (HPO/OMIM/MONDO terms) or `TEXT` (free text search).
- `phenotypenodeontologyterm_set`: Related set of specific ontology terms (HPO, OMIM, MONDO) with their match mode.

**Mechanism:**
1. Resolves ontology terms (HPO/OMIM/MONDO) to associated gene lists using the ontology database.
2. Free-text mode searches gene names, phenotype descriptions.
3. Builds a Q object on `VariantTranscriptAnnotation.gene_symbol` (or gene ID) using the `variant_transcript_annotation` annotation alias key.

#### `GeneListNode`
Filters variants to those falling within genes in one or more gene lists.

**Key Fields:**
- Links to `GeneList` objects (user-created or imported).
- Links to PanelApp panels (via PanelApp integration).

**Mechanism:** Resolves gene lists to sets of gene symbols/IDs, then builds a Q on the transcript annotation gene field.

#### `PopulationNode`
Filters variants by gnomAD (or other population database) allele frequency.

**Key Fields:**
- gnomAD AF thresholds (overall and per-population sub-group).
- Which gnomAD dataset (exomes, genomes).
- Frequency comparison operator (less than, greater than, etc.).

#### `ConservationNode`
Filters variants by evolutionary conservation scores (e.g. GERP, PhyloP, SiPhy).

#### `TissueNode`
Filters genes by tissue expression data. Allows restricting to variants in genes expressed in a specified tissue (e.g. using GTEx data).

#### `DamageNode`
Filters variants by in-silico pathogenicity prediction scores (e.g. CADD, SIFT, PolyPhen). Applies thresholds to prediction score columns in the annotation.

#### `AlleleFrequencyNode`
Filters variants by allele frequency from multiple population databases (gnomAD, UK Biobank, etc.) with a configurable filter structure. Distinct from `PopulationNode` in that it uses the `NodeAlleleFrequencyFilter` / `NodeAlleleFrequencyRange` mechanism.

#### `BuiltInFilterNode`
Applies predefined, named filter sets maintained by the system. Examples:
- `CLINVAR`: Restrict to variants with ClinVar significance annotations.
- `OMIM`: Restrict to variants in OMIM disease genes.

**Mechanism:** Maps built-in filter names to pre-defined Q object factory functions.

#### `MOINode`
Filters variants by Mode of Inheritance (MOI) using ontology terms (HP:000XXXX terms for inheritance modes).

**Mechanism:** Resolves MOI ontology terms to gene lists (via the ontology database which links MOI terms to OMIM gene-disease records) and filters accordingly.

#### `SelectedInParentNode`
Filters variants to only those that were manually selected (marked/highlighted) by the user in the parent node's result grid.

**Mechanism:** Queries a per-node selection state stored in the database and builds a Q on variant PKs.

---

### Multi-Input Nodes

These nodes combine or intersect results from multiple parent nodes.

#### `MergeNode`
Combines variants from multiple parent nodes using set union (OR logic).

**Configuration:**
- `min_inputs=1`, `max_inputs=PARENT_CAP_NOT_SET` (unlimited parents).

**Optimisation:** If the total combined count is at or below `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX`, the variant IDs are stored explicitly and the Q object becomes a simple `variant_id__in=[...]` lookup. For larger result sets, the individual parent querysets are OR-combined.

#### `VennNode`
Performs set operations on exactly two parent nodes.

**Configuration:**
- Requires exactly 2 parents.
- `use_cache=True`: Pre-computes and caches the three regions (A-only, B-only, intersection).

**Operations:**
- `INTERSECTION`: Variants in both A and B.
- `UNION`: Variants in either A or B.
- `A_NOT_B`: Variants in A but not B.
- `B_NOT_A`: Variants in B but not A.

**Caching:** Pre-computes a `VennNodeCache` (or `NodeCache` with `VariantCollection`) containing the set memberships, so that changing the operation (e.g. from INTERSECTION to A_NOT_B) only requires switching which pre-computed set to use, not re-running the parent queries.

#### `IntersectionNode`
Filters variants to those falling within specified genomic intervals (regions).

**Key Fields:**
- `genomic_intervals_collection` (FK to `GenomicIntervalsCollection`): A BED-format interval set.
- `hgvs_string` (str): HGVS coordinates for manual interval specification.
- `accordion_panel` (enum): `SELECTED`, `CUSTOM` (user-drawn intervals), `HGVS`, or `BACKEND_KIT` (sequencing kit capture regions).

**Caching:** `use_cache=True`. The `write_cache()` method performs a genomic interval intersection using a `bedtools intersect` subprocess call, piping variant coordinates through `intersectBed` to find overlapping variants. The result is stored as a `VariantCollection` in the database.

---

## 5. CACHING SYSTEM

The analysis app uses three distinct caching levels, each targeting a different granularity:

### Level 1: Redis Q-Object Cache

**What is cached:** The serialised `arg_q_dict` (composed Q object dictionary) for each node.

**Cache key format:** `{node_version_id}:q_cache={disable_cache}`

**Purpose:** Avoids recomputing the composition of parent and child Q objects on every grid page request. Since the Q composition can involve traversing the entire ancestor chain and merging complex Q objects, caching this result in Redis provides significant speed improvements for large graphs.

**Invalidation:** Automatic. When `bump_version()` creates a new `NodeVersion`, the old version ID is no longer valid and its Redis entries are simply never looked up again (they expire naturally via TTL, or are replaced when the new version is computed).

### Level 2: VariantCollection Cache (Database)

**What is cached:** The actual set of variant IDs (as a `VariantCollection`) for nodes with `use_cache=True`.

**Models involved:** `NodeCache` (OneToOne with `NodeVersion`), `VariantCollection`.

**Nodes using this level:** `VennNode`, `IntersectionNode` (when using BED file collections).

**Purpose:** For nodes requiring expensive operations (bedtools subprocess, complex set intersection of large variant collections), the resulting variant IDs are stored in the database as a `VariantCollection`. Subsequent queryset builds for this node simply query `variant_id IN (SELECT variant_id FROM variant_collection_records WHERE collection_id=...)`.

**Invalidation:** When `bump_version()` deletes the `NodeVersion`, the `NodeCache` is cascade-deleted. The associated `VariantCollection` is also cleaned up.

### Level 3: Count Cache (Database)

**What is cached:** `NodeCount` records storing variant counts by label for each `NodeVersion`.

**Models involved:** `NodeCount` (FK to `NodeVersion`), `NodeColumnSummaryCacheCollection`.

**Purpose:** Avoids re-running `COUNT(*)` queries against potentially large querysets on every UI refresh. The node count badges shown in the graph editor (showing how many variants pass each node) are served from this cache.

**Labels:** Configurable per-analysis via `AnalysisNodeCountConfiguration`. Common labels include total count, ClinVar pathogenic count, OMIM gene count.

**Invalidation:** Cascade delete from `NodeVersion` deletion.

---

## 6. EXECUTION FLOW

### Analysis Update Entry Point

```
update_analysis(analysis_id)
    → delete_analysis_old_node_versions   [cleanup stale tasks/versions]
    → create_and_launch_analysis_tasks    [single-worker serialised entry point]
```

`create_and_launch_analysis_tasks` runs on a single Celery worker (via dedicated queue or locking) to prevent race conditions when multiple concurrent events try to update the same analysis simultaneously (e.g. user edits node A and node B within milliseconds of each other).

### Task Graph Construction

Inside `create_and_launch_analysis_tasks`:
1. Performs a **topological sort** of the analysis node graph using the `toposort` library.
2. For each node that needs recomputation (status is DIRTY), creates a Celery task chain.
3. Nodes at the same topological level can execute in **parallel** (independent branches).
4. Nodes with parent dependencies are chained so they execute after their parents complete.

### Celery Tasks

Defined in `analysis/tasks/node_update_tasks.py`:

**`update_node_task`** (main computation task):
- Implements `AbortableTask` to support cancellation.
- Calls `node.get_queryset()` to build the filtered queryset.
- Executes the queryset to get variant counts.
- Stores results and updates node `status` to `READY` or `ERROR_TECHNICAL`.
- Updates `NodeVersion`, `NodeCount` entries.

**`node_cache_task`** (pre-caching):
- For nodes with `use_cache=True`.
- Calls the node's `write_cache()` method.
- Creates the `NodeCache` and underlying `VariantCollection`.
- Transitions node to `LOADING` status.

**`wait_for_cache_task`** (polling waiter):
- Polls for `NodeCache` processing completion.
- Maximum 60 checks before timing out.
- Re-queues itself if cache is not yet ready (asynchronous polling pattern).

**`wait_for_node`** (parent dependency waiter):
- Called at the start of a child node's chain when the parent node was also being recomputed.
- Waits until the parent node's `status` reaches `READY` before proceeding.
- Handles parent error states by marking child as `ERROR_WITH_PARENT`.

### Node Status Transitions

```
DIRTY → QUEUED           [task submitted]
QUEUED → LOADING_CACHE   [node_cache_task starts, use_cache=True nodes]
LOADING_CACHE → LOADING  [cache written, starting count computation]
QUEUED/LOADING → READY   [computation complete]
* → ERROR_CONFIGURATION  [invalid node configuration detected]
* → ERROR_WITH_PARENT    [parent node in error state]
* → ERROR_TECHNICAL      [unexpected exception during computation]
* → CANCELLED            [superseded by newer version task]
```

---

## 7. TEMPLATE SYSTEM

The template system allows saving analysis workflows as reusable, parameterised templates.

### Template Creation

1. User builds an analysis graph with nodes configured as desired.
2. Fields that should be parameterised (e.g. the `sample` field on a `SampleNode`) are marked as `AnalysisVariable` instances.
3. The `Analysis` is set to `template_type=TEMPLATE`, creating an `AnalysisTemplate` record.

### Versioning a Template

Calling `AnalysisTemplate.new_version()`:
1. Clones the current template analysis into an immutable **snapshot** (`template_type=SNAPSHOT`).
2. Creates an `AnalysisTemplateVersion` record pointing to the snapshot.
3. Sets the new version as `active=True`, deactivating the previous active version.
4. The snapshot cannot be modified; it serves as a reproducible baseline.

### Running a Template

`AnalysisTemplateRun.create(template_version, arguments)`:
1. Clones the snapshot analysis into a new, fully editable `Analysis`.
2. Creates an `AnalysisTemplateRun` record linking the template version to the new analysis.
3. Calls `populate_arguments()` which iterates over `AnalysisVariable` instances and sets the corresponding node fields to the provided values.
4. If any argument is invalid (wrong type, object not found, permission denied), stores the error in `AnalysisTemplateRunArgument.error`.
5. Triggers `update_analysis()` to compute the new analysis.

### Name Templating

The `analysis_name_template` field uses Python string formatting:

```python
"%(template)s for %(input)s" % {"template": template.name, "input": sample.name}
# Result: "Trio Analysis for Patient_001"
```

### Auto-Launch Templates

`AutoLaunchAnalysisTemplate` links an `AnalysisTemplate` to:
- An `enrichment_kit`: Only trigger for samples sequenced with this kit.
- A `sample_regex`: Only trigger for samples whose name matches this regex.

When a new VCF is imported, `auto_run_analyses_for_vcf` (Celery task) checks all `AutoLaunchAnalysisTemplate` records and automatically creates `AnalysisTemplateRun` instances for matching samples.

---

## 8. VIEWS & URLS

### Analysis List and Main View

- `analyses/list/` → `AnalysesGrid`: Paginated JQGrid listing all analyses visible to the user, with Guardian permission filtering.
- `<analysis_id>/` → `view_analysis`: The main analysis editor UI. Renders the graph canvas with nodes, edges, and the result data grid.

### Node Views

- `node/view/<node_id>/<version>/<extra_filters>/` → Node detail view. The `version` parameter allows the UI to detect stale views and prompt for refresh. `extra_filters` provides additional JQGrid filter context.
- `node/<id>/update/` → `NodeUpdate`: JSON endpoint for saving node configuration changes. Returns the new node state including version and status.
- `node/<id>/data/` → Node data endpoint for loading node-specific configuration data into the UI.

### Grid Endpoints

- `node_grid/handler/` → `NodeGridHandler`: The primary data endpoint for the JQGrid result table.
  - Implements **per-user locking** with a 10-minute lock to prevent duplicate expensive queries when the same user triggers multiple rapid page requests.
  - Caches grid responses for up to **1 week** (since results are deterministic for a given node version).
  - Returns paginated, sorted, filtered variant rows in JQGrid JSON format.
- `node_grid/cfg/` → `NodeGridConfig`: Returns the JQGrid column configuration for a specific node. Dynamic columns (e.g. sample zygosity columns for `CohortNode`, per-sample AF columns) are added based on the node type.
- `node_grid/export/` → Export endpoint. Supports CSV and VCF export formats via Celery background tasks.

### Template Views

- `analysis_templates/` → Template management list view.
- `analysis_template/<pk>/save/` → Save/version a template.
- `templates/variable/<node_id>/` → Manage `AnalysisVariable` entries for a node.

### JQGrid Integration

The `VariantGrid` class (in `analysis/views/`) orchestrates the JQGrid integration:

- `get_config()`: Produces the column configuration object consumed by the JQGrid JavaScript plugin. Includes:
  - Standard variant columns (chromosome, position, ref, alt, gene, consequence, etc.).
  - Annotation columns (ClinVar, gnomAD, etc.).
  - Node-type-specific dynamic columns (zygosity per sample, AF columns, count columns).
- `get_data()`: Executes the node's queryset with pagination, sorting, and filter parameters from the JQGrid request. Returns serialised row data.

---

## 9. AUDIT LOG

The analysis app uses `django-auditlog` to record all significant mutations.

**Registered models:**
- `Analysis`: Records analysis creation, updates, name changes, locking.
- `AnalysisEdge`: Records edge additions and removals (graph topology changes).
- `AnalysisNode`: Records node creation, configuration changes, deletions.

**`NodeAuditLogMixin`**: A mixin applied to node-modifying views that enriches audit log entries with additional context fields:
- `analysis_id`: The analysis containing the node.
- `node_id`: The specific node being modified.

**Audit Log URL:** `node_audit_log/<node_id>/<version>/<extra_filters>/` renders a paginated history of audit events for a specific node, useful for debugging and change tracking.

---

## 10. CELERY TASKS

All Celery tasks are defined under `analysis/tasks/`.

### Node Update Tasks (`node_update_tasks.py`)

| Task | Description |
|---|---|
| `update_node_task` | Main node computation. `AbortableTask`. Executes the node's queryset, computes counts, updates status. |
| `node_cache_task` | Pre-caches a `VariantCollection` for nodes with `use_cache=True`. |
| `wait_for_cache_task` | Polling waiter (max 60 iterations) that waits for `NodeCache` to finish processing. |
| `wait_for_node` | Dependency waiter that blocks a child node's chain until its parent node reaches `READY` status. |
| `delete_analysis_old_node_versions` | Cleans up stale `NodeVersion`, `NodeCache`, `NodeCount` entries from superseded task runs. |

### Analysis Update Tasks (`analysis_update_tasks.py`)

| Task | Description |
|---|---|
| `create_and_launch_analysis_tasks` | Serialised single-worker entry point. Performs topological sort and dispatches node task chains. |

### Variant Tag Tasks (`variant_tag_tasks.py`)

| Task | Description |
|---|---|
| `variant_tag_created_task` | Called when a `VariantTag` is created. Invalidates `TagNode` caches in affected analyses. |
| `variant_tag_deleted_in_analysis_task` | Called when a `VariantTag` is deleted. Invalidates `TagNode` caches. |

### Auto-Analysis Tasks (`auto_analysis_tasks.py`)

| Task | Description |
|---|---|
| `auto_run_analyses_for_vcf` | Triggered on VCF import. Checks `AutoLaunchAnalysisTemplate` records and creates runs for matching samples. |
| `auto_run_analyses_for_sample` | Triggered on sample creation. Similar logic for sample-triggered auto-launch. |

### Export Tasks (`analysis_grid_export_tasks.py`)

| Task | Description |
|---|---|
| Grid export tasks | Handles asynchronous CSV and VCF export of node result sets. Streams large result sets to file without blocking the web process. |

---

## 11. KEY FILE LOCATIONS

```
analysis/
├── models/
│   ├── models.py                    # Analysis, AnalysisLock, AnalysisTemplate*, AnalysisVariable
│   ├── analysis_node.py             # AnalysisNode base class, AnalysisEdge, NodeVersion, NodeCache, NodeCount
│   └── nodes/
│       ├── source_nodes/
│       │   ├── sample_node.py       # SampleNode
│       │   ├── trio_node.py         # TrioNode
│       │   ├── cohort_node.py       # CohortNode
│       │   ├── pedigree_node.py     # PedigreeNode
│       │   └── ...
│       └── filter_nodes/
│           ├── filter_node.py       # FilterNode (JQGrid JSON → Q)
│           ├── phenotype_node.py    # PhenotypeNode
│           ├── gene_list_node.py    # GeneListNode
│           ├── venn_node.py         # VennNode
│           ├── intersection_node.py # IntersectionNode
│           └── ...
├── tasks/
│   ├── node_update_tasks.py         # update_node_task, node_cache_task, wait_for_*
│   ├── analysis_update_tasks.py     # create_and_launch_analysis_tasks
│   ├── variant_tag_tasks.py         # Tag invalidation tasks
│   ├── auto_analysis_tasks.py       # AutoLaunchAnalysisTemplate tasks
│   └── analysis_grid_export_tasks.py
├── views/
│   ├── views.py                     # view_analysis, AnalysesGrid, template views
│   ├── node_view.py                 # Node update/view endpoints
│   └── analysis_grid_view.py        # NodeGridHandler, NodeGridConfig, VariantGrid
├── urls/
│   └── urls.py                      # URL routing
└── admin.py                         # Django admin registrations
```

---

## 12. DESIGN PATTERNS AND ARCHITECTURAL NOTES

### DAG-Based Workflow

The use of `django-dag` (directed acyclic graph) for the node graph structure means:
- Parent-child relationships are stored as `AnalysisEdge` records.
- The graph is validated to be acyclic (no cycles allowed).
- Topological sort is straightforward and used directly in task scheduling.
- Each node can have multiple parents (for `MergeNode`, `VennNode`) or multiple children (branching workflows).

### Queryset Composition vs. Materialisation

A key design decision is that most nodes compose Q objects rather than materialising intermediate results. This means:
- For linear chains, the entire chain's filter logic is expressed as a single Django ORM queryset with multiple `.filter()` calls and annotations.
- The database handles the combined filtering in a single SQL query (with multiple JOINs/WHERE clauses).
- Only `VennNode` and `IntersectionNode` materialise intermediate results (as `VariantCollection` records) because their logic cannot be expressed efficiently as a single SQL query.

### Version-Based Cache Invalidation

The `version` integer on `AnalysisNode` + `NodeVersion` provides a simple but effective cache invalidation mechanism:
- No explicit cache clearing calls are needed for Q-object cache entries in Redis.
- No explicit cleanup queries are needed for `NodeCount` records.
- Simply deleting the `NodeVersion` record cascades all dependent cache records away.
- The cache key incorporates the version ID, so stale entries are simply never referenced.

### Single-Worker Serialisation for Task Launch

The `create_and_launch_analysis_tasks` task runs on a dedicated single-worker queue. This prevents the race condition where:
1. User edits node A → triggers update.
2. User edits node B (milliseconds later) → triggers second update.
3. Both updates run concurrently and try to create overlapping task chains.

By serialising the task-creation step, the second update sees the state left by the first and produces a coherent, non-overlapping task graph.

### Guardian Object Permissions

`Analysis` objects use Django Guardian for row-level permissions, allowing fine-grained sharing of analyses between users while preventing unauthorised access to genomic data.

### AbortableTask for Node Computation

Node computation tasks implement `AbortableTask`, which allows them to be cancelled when a node's version is bumped (i.e. the user makes another configuration change while a computation is still running). The task periodically checks whether it has been aborted and exits cleanly if so, rather than completing and writing stale results.
