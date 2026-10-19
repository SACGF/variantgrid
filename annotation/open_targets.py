"""
Display parsing for the OpenTargets VEP plugin columns (open_targets_* on VariantAnnotation): the plugin '&'-joins
one value per association into parallel arrays, which open_targets_records zips back into records, and
open_targets_gwas_genes / open_targets_qtl_genes / open_targets_variant shape them for the variant page. Each takes
a VariantAnnotation; VariantAnnotation keeps same-named properties for templates. Tests: annotation/tests/test_open_targets.py
"""
import logging
from collections.abc import Iterable
from typing import TYPE_CHECKING, Optional

from django.urls import reverse

from genes.models import GeneVersion
from ontology.models import OntologyIdNormalized, OntologyTerm

if TYPE_CHECKING:
    from annotation.models.models import VariantAnnotation

OPEN_TARGETS_NA = "NA"  # OpenTargets.pm pads every column for records that have no value for it
OPEN_TARGETS_GWAS = "gwas"
# Model field -> key in a zipped record. The plugin emits every column for every association, so
# these are parallel arrays - @see open_targets_records
OPEN_TARGETS_RECORD_FIELDS = {
    "open_targets_study_type": "study_type",
    "open_targets_study_id": "study_id",
    "open_targets_is_lead": "is_lead",
    "open_targets_variant_id": "variant_id",
    "open_targets_gwas_gene_id": "gwas_gene_id",
    "open_targets_gwas_l2g_scores": "l2g_score",
    "open_targets_qtl_gene_id": "qtl_gene_id",
    "open_targets_qtl_biosample": "qtl_biosample",
}
OPEN_TARGETS_STUDY_TYPE_LABELS = {
    "gwas": "GWAS", "eqtl": "eQTL", "pqtl": "pQTL", "sqtl": "sQTL", "tuqtl": "tuQTL", "sceqtl": "sc-eQTL",
}
OPEN_TARGETS_URL = "https://platform.opentargets.org/"
# Disease ids are ontology CURIEs (EFO, OBA, MONDO, HP, GO, Orphanet, NCIT...) - all have an Open
# Targets disease page, and these two we can also name from our own ontology
OPEN_TARGETS_LOCAL_ONTOLOGY_PREFIXES = ("MONDO_", "HP_")


def _split(value: Optional[str]) -> list[Optional[str]]:
    """ Split one '&'-joined Open Targets column into per-record values. VEP escapes ',' and '|'
        inside a CSQ value to '&' - the same character the plugin joins records with - but a comma
        is always followed by a space, which VEP escapes to '_', so a fragment starting with '_'
        is a continuation of the one before it (e.g. 'CD4-positive, alpha-beta T cell'). """
    if not value:
        return []
    values = []
    for part in value.split("&"):
        if part.startswith("_") and values:
            values[-1] = f"{values[-1]}, {part[1:]}"
        else:
            values.append(part)
    return [None if v in ("", OPEN_TARGETS_NA) else v for v in values]


def _record_diseases(variant_annotation: 'VariantAnnotation', study_types: list[Optional[str]]) -> list[list[str]]:
    """ Diseases are the one Open Targets column holding several values per record, '|'-separated
        before VEP escapes them to '&', so the array can be longer than the record count. Only GWAS
        records carry diseases, so the 'NA' of each QTL record bounds a run of consecutive GWAS
        records: where a run's disease count matches its record count they line up exactly,
        otherwise the run's diseases go to every gene in it (genes sharing a study locus share
        its diseases). """
    diseases: list[list[str]] = [[] for _ in study_types]
    values = _split(variant_annotation.open_targets_gwas_diseases)
    if not values:
        return diseases

    value_index = 0
    record_index = 0
    while record_index < len(study_types) and value_index < len(values):
        if study_types[record_index] != OPEN_TARGETS_GWAS:
            value_index += 1  # the record's 'NA'
            record_index += 1
            continue

        run_start = record_index
        while record_index < len(study_types) and study_types[record_index] == OPEN_TARGETS_GWAS:
            record_index += 1
        run_values = []
        while value_index < len(values) and values[value_index] is not None:
            run_values.append(values[value_index])
            value_index += 1

        run_length = record_index - run_start
        if len(run_values) == run_length:
            for i, disease in enumerate(run_values):
                diseases[run_start + i] = [disease]
        else:
            for i in range(run_start, record_index):
                diseases[i] = list(run_values)
    return diseases


def open_targets_records(variant_annotation: 'VariantAnnotation') -> list[dict]:
    """ The OpenTargets plugin emits one record per overlapping association, every requested column
        for every record, '&'-joined with 'NA' where a record has no value - so the open_targets_*
        fields are parallel arrays. Zip them back into one dict per association.
        open_targets_gwas_l2g_scores is null until backfilled (#1822), in which case records simply
        have no score. """
    study_types = _split(variant_annotation.open_targets_study_type)
    if not study_types:
        return []

    records = [{"diseases": d} for d in _record_diseases(variant_annotation, study_types)]
    for field, key in OPEN_TARGETS_RECORD_FIELDS.items():
        values = _split(getattr(variant_annotation, field))
        if values and len(values) != len(study_types):
            logging.error("%s has %d values for %d Open Targets records",
                          field, len(values), len(study_types))
            values = []
        for i, record in enumerate(records):
            record[key] = values[i] if i < len(values) else None

    for record in records:
        if biosample := record["qtl_biosample"]:
            record["qtl_biosample"] = biosample.replace("_", " ")  # VEP escapes whitespace
        if score := record["l2g_score"]:
            record["l2g_score"] = float(score)
        if (is_lead := record["is_lead"]) is not None:
            record["is_lead"] = is_lead == "true"  # None where the column predates the record
    return records


def _gene_details(variant_annotation: 'VariantAnnotation', gene_ids: Iterable[str]) -> dict[str, dict]:
    """ Open Targets is GRCh38/Ensembl only and Gene.pk is the Ensembl gene id, so look the symbols
        up in one query, falling back to an Ensembl link for genes we don't have locally """
    gene_ids = list(gene_ids)
    symbols = {}
    gv_qs = GeneVersion.objects.filter(gene__in=gene_ids, genome_build=variant_annotation.version.genome_build)
    for gene_id, gene_symbol_id in gv_qs.order_by("version").values_list("gene_id", "gene_symbol_id"):
        symbols[gene_id] = gene_symbol_id  # latest version wins

    details = {}
    for gene_id in gene_ids:
        if gene_id in symbols:
            gene_url = reverse("view_gene", kwargs={"gene_id": gene_id})
        else:
            gene_url = f"https://ensembl.org/Homo_sapiens/Gene/Summary?g={gene_id}"
        details[gene_id] = {
            "gene_id": gene_id,
            "gene_symbol": symbols.get(gene_id) or gene_id,
            "gene_url": gene_url,
        }
    return details


def _diseases(disease_ids: Iterable[str]) -> list[dict]:
    """ Name the MONDO / HP diseases from our own ontology - the rest display as their raw id """
    disease_ids = sorted(set(disease_ids))
    ontology_term_ids = {}
    for disease_id in disease_ids:
        if disease_id.startswith(OPEN_TARGETS_LOCAL_ONTOLOGY_PREFIXES):
            try:
                ontology_term_ids[disease_id] = str(OntologyIdNormalized.normalize(disease_id))
            except ValueError:
                pass

    names = {}
    if ontology_term_ids:
        ontology_names = dict(OntologyTerm.objects.filter(pk__in=ontology_term_ids.values())
                              .values_list("pk", "name"))
        for disease_id, term_id in ontology_term_ids.items():
            if name := ontology_names.get(term_id):
                names[disease_id] = name

    return [{"id": disease_id,
             "name": names.get(disease_id, disease_id),
             "url": f"{OPEN_TARGETS_URL}disease/{disease_id}"} for disease_id in disease_ids]


def open_targets_gwas_genes(variant_annotation: 'VariantAnnotation') -> list[dict]:
    """ GWAS records collapsed to a row per implicated gene - highest L2G score, de-duplicated
        diseases, and how many studies contributed. Sorted by score, unscored last. """
    genes = {}  # keyed by gene id, insertion ordered
    for record in variant_annotation.open_targets_records:
        if record["study_type"] != OPEN_TARGETS_GWAS:
            continue
        if gene_id := record["gwas_gene_id"]:
            gene = genes.setdefault(gene_id, {"scores": [], "diseases": set(), "study_ids": set(),
                                              "lead_study_ids": set()})
            if (score := record["l2g_score"]) is not None:
                gene["scores"].append(score)
            gene["diseases"].update(record["diseases"])
            if study_id := record["study_id"]:
                gene["study_ids"].add(study_id)
                if record["is_lead"]:
                    gene["lead_study_ids"].add(study_id)

    gene_details = _gene_details(variant_annotation, genes)
    gwas_genes = []
    for gene_id, gene in genes.items():
        gwas_genes.append({
            **gene_details[gene_id],
            "l2g_score": max(gene["scores"]) if gene["scores"] else None,
            "diseases": _diseases(gene["diseases"]),
            "study_count": len(gene["study_ids"]),
            "lead_study_count": len(gene["lead_study_ids"]),
        })
    gwas_genes.sort(key=lambda g: (g["l2g_score"] is None, -(g["l2g_score"] or 0)))
    return gwas_genes


def open_targets_qtl_genes(variant_annotation: 'VariantAnnotation') -> list[dict]:
    """ QTL (non-GWAS) records collapsed to a row per (gene, study type), with the biosamples the
        association was measured in """
    genes = {}  # keyed by (gene id, study type)
    for record in variant_annotation.open_targets_records:
        study_type = record["study_type"]
        if not study_type or study_type == OPEN_TARGETS_GWAS:
            continue
        if gene_id := record["qtl_gene_id"]:
            gene = genes.setdefault((gene_id, study_type), {"biosamples": set(), "study_ids": set(),
                                                            "lead_study_ids": set()})
            if biosample := record["qtl_biosample"]:
                gene["biosamples"].add(biosample)
            if study_id := record["study_id"]:
                gene["study_ids"].add(study_id)
                if record["is_lead"]:
                    gene["lead_study_ids"].add(study_id)

    gene_details = _gene_details(variant_annotation, {gene_id for gene_id, _ in genes})
    qtl_genes = []
    for (gene_id, study_type), gene in genes.items():
        qtl_genes.append({
            **gene_details[gene_id],
            "study_type": study_type,
            "study_type_label": OPEN_TARGETS_STUDY_TYPE_LABELS.get(study_type, study_type),
            "biosamples": sorted(gene["biosamples"]),
            "study_count": len(gene["study_ids"]),
            "lead_study_count": len(gene["lead_study_ids"]),
        })
    qtl_genes.sort(key=lambda g: (g["gene_symbol"], g["study_type"]))
    return qtl_genes


def open_targets_variant(variant_annotation: 'VariantAnnotation') -> Optional[dict]:
    """ The variant id is repeated once per record - they all describe this variant """
    for variant_id in _split(variant_annotation.open_targets_variant_id):
        if variant_id:
            return {"id": variant_id, "url": f"{OPEN_TARGETS_URL}variant/{variant_id}"}
    return None
