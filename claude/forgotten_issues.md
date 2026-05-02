# Forgotten Issues — Top 20 by Impact per Effort

**Written:** 2026-03-26
**Method:** Scanned all 766 open issues across SACGF/variantgrid and SACGF/variantgrid_private. Selected old issues (opened 2022 or earlier) with zero or few comments that appear to have slipped off the radar.

---

## 1. variantgrid#170 — Import/Export Analyses via Web UI
**Opened:** Jan 2021 · **Comments:** 0 · **Effort:** Very Low

The management command to export/import analyses already exists (`analysis_import_export`). There is no web UI to trigger it. Scientists can't move an analysis from one server to another (e.g. from SA Path to Shariant) without admin CLI access.

**Why forgotten:** It's one of those "command already exists, just needs a button" tasks that never quite rises to the top.

---

## 2. variantgrid#688 — Redis Cache Bloat (16 GB Unnecessary Usage)
**Opened:** Sep 2022 · **Comments:** 0 · **Effort:** Low

Redis was found to be using ~16 GB unnecessarily on production. Celery task results and old node caches were identified as suspects. Cache TTLs were likely never tuned.

**Why forgotten:** Infrastructure issues that don't cause an outage tend to sit.

---

## 3. variantgrid#658 — Link to AlphaFold Predictions
**Opened:** Jul 2022 · **Comments:** 0 · **Effort:** Trivial

Add a link to `alphafold.ebi.ac.uk/entry/{UniProt_ID}` on variant detail pages for missense variants. AlphaFold is now the de facto first-look tool for structural context. One URL template in the variant page template.

**Why forgotten:** It's so small it never got scheduled.

---

## 4. variantgrid#381 — Analysis Serialisation Bug (Phenotype/Gene List Nodes)
**Opened:** Jun 2021 · **Comments:** 0 · **Effort:** Low-Medium

`PhenotypeNode` attempts to recreate ontology terms on deserialisation rather than reference them by ID, causing failures. `GeneListNode` has issues with PanelApp panels and user-owned gene lists. Analyses with these nodes don't round-trip cleanly through export/import or template creation.

**Why forgotten:** The symptom (broken template/import) surfaces rarely; most scientists build analyses from scratch each time.

---

## 5. variantgrid#389 — ClassificationNode: Add ClinVar and Internal Lab Filters
**Opened:** Jun 2021 · **Comments:** 0 · **Effort:** Low

`ClassificationNode` currently returns classifications from all labs. Proposed enhancements:
- Filter by ClinVar as a "lab" source
- Add autocomplete to select specific internal labs
- Clean up how ClinVar submitter IDs are displayed

Useful for validation analyses: "show me variants where our classification disagrees with ClinVar's."

**Why forgotten:** ClassificationNode is a less-used node type so issues against it accumulate silently.

---

## 6. variantgrid#558 — Hide Auto-Generated "Sample Tab" Analyses from Lists
**Opened:** Jan 2022 · **Comments:** 0 · **Effort:** Very Low

When a scientist opens the Sample tab, VariantGrid auto-creates an analysis behind the scenes. These auto-analyses pollute the "Recent Analyses" list and confuse users who see analyses they don't remember creating. Also: MOI nodes should link back to the sample/patient they're filtering for.

**Why forgotten:** Low-severity UX nuisance that's easy to live with.

---

## 7. variantgrid#207 — Automated Allele Linking Across Genome Builds
**Opened:** Feb 2021 · **Comments:** 0 · **Effort:** Low

When a VCF is imported for one genome build and a matching allele exists in the other build (via ClinGen Allele Registry), the link isn't automatically created. A nightly Celery task could check newly imported variants against ClinGen and create missing `VariantAllele` links, improving cross-build allele resolution without any user action.

**Why forgotten:** The manual path (ClinGen lookup on demand) works, so the proactive linking never got prioritised.

---

## 8. variantgrid#299 / #607 — Cache Analysis Node Counts by SQL Hash
**Opened:** Mar 2021 · **Comments:** 0 · **Effort:** Medium

Currently node counts invalidate when any configuration changes, even if the resulting SQL is identical. Caching by a hash of the actual SQL (parent input + node config → query → hash) would mean switching between node configurations reuses cached counts instead of recomputing. Scientists who toggle filters frequently would see dramatic speedups.

Two separate issues (#299 and #607) describe the same root fix — cache by SQL not by version number.

**Why forgotten:** Requires touching the node caching architecture, which is intimidating enough to keep deferring.

---

## 9. variantgrid#322 — Multi-VCF Source Node
**Opened:** Apr 2021 · **Comments:** 0 · **Effort:** Medium

Currently, to analyse multiple VCFs together, scientists must create a custom cohort first — a slow and cumbersome process. A "Multi-VCF Source Node" would let them select multiple VCFs directly as the source for an analysis. Noted as needed "weekly" by the genomics analysis team.

**Why forgotten:** It's a meaningful new node type, so it feels bigger than it probably is given the existing DAG infrastructure.

---

## 10. variantgrid#140 — Permissions of Objects Used Inside an Analysis
**Opened:** Nov 2020 · **Comments:** 0 · **Effort:** Low-Medium (design + implementation)

If a scientist shares an analysis with a colleague, the colleague can see results that depend on gene lists, samples, or classifications the colleague doesn't have direct permission to view. The behaviour is undefined and inconsistent. Needs a policy decision: does analysis-level access grant read access to constituent objects, or should constituent objects be checked individually?

**Why forgotten:** It's a security/design question without an obvious "right" answer, so it sits.

---

## 11. variantgrid#508 — AlleleFrequency Node: Filter by Cohort MAF
**Opened:** Nov 2021 · **Comments:** 0 · **Effort:** Low

The `AlleleFrequencyNode` currently filters on gnomAD/external population frequencies. Extend it to also support filtering by `GlobalVariantZygosityCount` (internal cohort frequency). This lets scientists filter out variants common in their own database without leaving the analysis pipeline.

**Why forgotten:** `GlobalVariantZygosityCount` existed; nobody circled back to wire it into the frequency filter node.

---

## 12. variantgrid_private#3388 — Lab-Level "Lab Specific Explain" Text
**Opened:** Oct 2022 · **Comments:** 0 · **Effort:** Trivial

Allow each lab to store a short explanation text that appears on the diff/comparison page alongside their classifications. Useful for explaining lab-specific methodology quirks ("we use REVEL ≥0.75 rather than ≥0.5 for PP3") so other labs understand systematic differences without having to email.

**Why forgotten:** Small convenience feature; easy to defer.

---

## 13. variantgrid_private#3370 — Discordance "Since Date" Misses Active→Continued Transitions
**Opened:** Sep 2022 · **Comments:** 0 · **Effort:** Low

When a `DiscordanceReport` transitions from Active to Continued status, this change is not picked up by the "since date" filter. Labs filtering for recent discordance activity miss these transitions. The fix is to include `DiscordanceReport` modification timestamps in the since-date calculation.

**Why forgotten:** Edge-case bug in a moderately complex workflow; not reported by users because they likely don't know what they're missing.

---

## 14. variantgrid_private#3368 — Update Ontology from the Annotation Screen
**Opened:** Sep 2022 · **Comments:** 0 · **Effort:** Low

Currently ontology updates require downloading files manually and running management commands via CLI. The annotation screen could offer a "Check for updates / Re-import" button that wraps these commands. Most ontology source URLs are stable; the admin just needs to trigger the import.

**Why forgotten:** Admins have the CLI workaround; no user pressure to build the UI.

---

## 15. variantgrid_private#2752 — Mark Evidence Keys as Deprecated
**Opened:** Aug 2020 · **Comments:** 0 · **Effort:** Very Low

Add a `deprecated` flag to `EvidenceKey`, distinct from `hidden`. Hidden means "don't show by default"; deprecated means "this key is obsolete — show a warning if it's populated, don't show it in empty forms." Useful when retiring a criterion that some historical classifications still reference.

**Why forgotten:** Evidence key management is only touched when setting up a new deployment; deprecation was never needed urgently.

---

## 16. variantgrid_private#3341 — Better Display of Withdrawn ClinVar Exports
**Opened:** Aug 2022 · **Comments:** 0 · **Effort:** Low

Withdrawn `ClinVarExport` records show "?" in the allele column because the classification has been withdrawn. They should display the cached HGVS from when the export was active, allow admins to review the SCV without a linked classification, and make the withdrawn state visually clear rather than confusing.

**Why forgotten:** Withdrawn records are rare and only visible to admins.

---

## 17. variantgrid#414 — Annotation Page Toggle Shows Spurious Error
**Opened:** Jun 2021 · **Comments:** 0 · **Effort:** Trivial

Toggling the annotation page displays an error message mixed in with normal instructions. Confuses new users into thinking something is wrong. Fix is to separate error states from informational messages in that view.

**Why forgotten:** Classic "not a real bug, just looks bad" issue.

---

## 18. variantgrid#117 — Link Variants to IGV / UCSC Genome Browser
**Opened:** Nov 2020 · **Comments:** 0 · **Effort:** Low

Add links from variant pages and analysis grid rows to open the locus in IGV (via igv.js local launch URL) or UCSC Genome Browser. Useful for manual validation of calls, reviewing nearby variants, and checking repetitive regions. Both IGV and UCSC support URL-based locus navigation.

**Why forgotten:** Requires knowing the IGV launch URL scheme; feels harder than it is.

---

## 19. variantgrid#661 — User-Acknowledged System Messages
**Opened:** Aug 2022 · **Comments:** 0 · **Effort:** Low

When system messages are sent (upgrade notices, feature announcements, scheduled downtime), there's no way for individual users to dismiss/acknowledge them. Messages either persist for everyone or disappear for everyone. A simple per-user acknowledgement record would allow "dismiss this for me" without removing it for others.

**Why forgotten:** Django messages library was being replaced; the replacement never included acknowledgement.

---

## 20. variantgrid#635 — Transcripts Page: Show Gene Annotation Release
**Opened:** Jun 2022 · **Comments:** 0 · **Effort:** Trivial

The transcripts page doesn't show which gene annotation release (Ensembl/RefSeq version, date) the transcript data came from. When investigating why a transcript is or isn't present, or why a canonical transcript changed, this context is essential. The data exists in the `AnnotationVersion` model.

**Why forgotten:** Documentation/display issue with no obvious breakage; easy to ignore.

---

## Quick-Win Summary (could be done in a day or less)

| # | Issue | What it takes |
|---|-------|--------------|
| 3 | AlphaFold link | One URL template line |
| 6 | Hide auto-analyses | Filter one queryset + add patient link to MOI node |
| 12 | Lab explain text | One TextField on Lab model + display on diff page |
| 15 | Deprecated evidence keys | One BooleanField + form filter |
| 17 | Annotation page error | Fix one template/view message |
| 20 | Transcript annotation release | Display one existing field |
