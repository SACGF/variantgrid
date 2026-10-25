"""
Display parsing for VEP plugin columns stored on VariantAnnotation: ProtVar (protvar_int / protvar_pocket), MaveDB
and Mastermind links, denovo-db records, UTRAnnotator uORF effects and gnomAD-SV overlaps. Plugins '&'-join one
value per record, so these split them back into records (vep_multi_fields_to_list_of_dicts, which the VEP inserter also uses on import).
Functions take a VariantAnnotation or the raw values; VariantAnnotation keeps same-named properties for templates.
Open Targets is annotation/open_targets.py.
"""
import logging
from collections.abc import Iterable
from typing import TYPE_CHECKING, Optional

from more_itertools import all_equal

from library.genomics import parse_gnomad_coord

if TYPE_CHECKING:
    from annotation.models.models import VariantAnnotation


GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS = [
    'gnomad_sv_overlap_af',
    'gnomad_sv_overlap_name',
    'gnomad_sv_overlap_percent',
    'gnomad_sv_overlap_coords',
]


_PROTVAR_CONFIDENCE_CSS = {"low": "secondary", "high": "info", "very high": "success"}


def _protvar_confidence(value, *, low_below: float, high_max: float) -> Optional[dict]:
    """ ProtVar confidence band for a numeric score (pDockQ or pocket score), matching the ProtVar.pm
        bands: below `low_below` is 'low', up to and including `high_max` is 'high', above is
        'very high'. Returns {label, css} for badge display, or None if not numeric. """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if v < low_below:
        label = "low"
    elif v <= high_max:
        label = "high"
    else:
        label = "very high"
    return {"label": label, "css": _PROTVAR_CONFIDENCE_CSS[label]}


def protvar_int_parsed(variant_annotation: 'VariantAnnotation') -> Optional[dict]:
    """ ProtVar interface = 'protein_id&pDockQ' - a single interacting protein and its docking
        confidence (e.g. 'Q6PJI9-2&0.420206'). Split into the interacting protein id + pDockQ,
        with a UniProt entry link for the protein (isoform suffix dropped for the URL) and the
        plugin's pDockQ confidence band. Returns None when there's no interface value. """
    if not variant_annotation.protvar_int:
        return None
    protein_id, _, pdockq = variant_annotation.protvar_int.partition("&")
    accession = protein_id.split("-")[0]  # strip isoform suffix (e.g. Q6PJI9-2 -> Q6PJI9) for the link
    return {
        "protein_id": protein_id,
        "pdockq": pdockq,
        "confidence": _protvar_confidence(pdockq, low_below=0.23, high_max=0.5),  # ProtVar.pm pDockQ bands
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry" if accession else "",
    }


def protvar_pocket_parsed(variant_annotation: 'VariantAnnotation') -> Optional[dict]:
    """ ProtVar pocket = 'id&score&MpLDDT&energy&buriedness&RoG&residues' (e.g.
        'P123&850&95.2&-1.3&0.7&12.4&p10p20p30'). Residues are 'p'-encoded by ProtVar.pm
        (comma -> 'p') to survive the VCF field split; we decode them back to a position list.
        Returns None when there's no pocket value. """
    if not variant_annotation.protvar_pocket:
        return None
    parts = variant_annotation.protvar_pocket.split("&")
    parts += [""] * (7 - len(parts))  # tolerate a short/truncated blob
    pocket_id, score, mplddt, energy, buriedness, rog, residues = parts[:7]
    residue_positions = [r for r in residues.lstrip("p").split("p") if r]
    return {
        "id": pocket_id,
        "score": score,
        "confidence": _protvar_confidence(score, low_below=800, high_max=900),  # ProtVar.pm pocket bands
        "mplddt": mplddt,
        "energy": energy,
        "buriedness": buriedness,
        "rog": rog,
        "residues": ", ".join(residue_positions),
    }


def mavedb_urn_to_urls(mavedb_urn: str) -> dict[str, str]:
    experiment_sets = set()
    EXPERIMENT_URL_PREFIX = "https://www.mavedb.org/#/experiment-sets/"
    if mavedb_urn:
        for urn in mavedb_urn.split("&"):
            es = urn.rsplit("-", 2)[0]
            experiment_sets.add(es)
    return {urn: EXPERIMENT_URL_PREFIX + urn for urn in sorted(experiment_sets)}


def get_mmid3_mastermind_urls(mastermind_mmid3) -> dict:
    mmid3_mastermind_urls = {}
    if mastermind_mmid3:
        for mmid3 in mastermind_mmid3.split("&"):
            mmid3_mastermind_urls[mmid3] = f"https://mastermind.genomenon.com/detail?mutation={mmid3}"
    return mmid3_mastermind_urls


def denovo_db_records(variant_annotation: 'VariantAnnotation') -> list[dict]:
    """ Parse the '&'-separated parallel arrays (one element per contributing
        denovo-db record) and collapse identical (study, phenotype, PubMed)
        entries into a single row with an occurrence 'count'.
        See https://denovo-db.gs.washington.edu/ """
    if not variant_annotation.denovo_db_studies:
        return []

    studies = variant_annotation.denovo_db_studies.split("&")
    pubmed_ids = (variant_annotation.denovo_db_pubmed_ids or "").split("&")
    phenotypes = (variant_annotation.denovo_db_primary_phenotypes or "").split("&")

    records = {}  # keyed by (study, phenotype, pubmed_id), insertion-ordered
    for i, study in enumerate(studies):
        pubmed_id = pubmed_ids[i] if i < len(pubmed_ids) else ""
        phenotype = phenotypes[i] if i < len(phenotypes) else ""
        key = (study, phenotype, pubmed_id)
        record = records.get(key)
        if record is None:
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/" if pubmed_id else ""
            records[key] = {
                "study": study,
                "pubmed_id": pubmed_id,
                "pubmed_url": pubmed_url,
                "primary_phenotype": phenotype,
                "count": 1,
            }
        else:
            record["count"] += 1
    return list(records.values())


def vep_multi_fields_to_list_of_dicts(data: dict, fields: Iterable[str]) -> list[dict]:
    list_values = {}
    for field in fields:
        if value := data.get(field):
            list_values[field] = value.split("&")

    if not list_values:
        return []

    # Ensure they are all same length
    first_field_len = len(next(iter(list_values.values())))
    field_lengths = {k: len(v) for k, v in list_values.items()}

    if not all_equal(field_lengths.values()):
        logging.error("All split VEP multi-values must be equal length")
        for f, l in field_lengths.items():
            logging.error("%s: %d", f, l)
        raise ValueError(
            f"All split multi-values must be equal length: {field_lengths}")

    records: list[dict] = []
    for i in range(first_field_len):
        values_dict = {}
        for field in fields:
            if v := list_values.get(field):
                values_dict[field] = v[i]
        records.append(values_dict)
    return records


# UTRAnnotator writes each record's pairs from a Perl hash, so in a different order every time
_UTR_ANNOTATOR_KEY_ORDER = ["type", "ref_type", "alt_type", "KozakContext", "KozakStrength", "DistanceToCDS",
                            "ref_StartDistanceToCDS", "CapDistanceToStart", "DistanceToStop", "ref_type_length",
                            "alt_type_length", "FrameWithCDS", "AltStop", "AltStopDistanceToCDS", "Evidence"]


def _utr_annotator_key_sort(key: str) -> tuple[int, str]:
    try:
        return _UTR_ANNOTATOR_KEY_ORDER.index(key), key
    except ValueError:
        return len(_UTR_ANNOTATOR_KEY_ORDER), key


def utr_annotator_effects(utr_annotator_annotation: Optional[str]) -> list[dict[str, str]]:
    """ UTRAnnotator's 5UTR_annotation: one record per uORF effect, comma separated (which VEP's VCF
        output turns into '&'), each record 'key=value' pairs joined by ':', keys put in one order so the
        records read alike """
    effects = []
    for record in (utr_annotator_annotation or "").split("&"):
        if pairs := [p.partition("=") for p in record.split(":") if p]:
            effects.append({key: value for key, _, value in sorted(pairs, key=lambda p: _utr_annotator_key_sort(p[0]))})
    return effects


def get_gnomad_url(gnomad_variant: str, dataset: str) -> str:
    return f"http://gnomad.broadinstitute.org/variant/{gnomad_variant}?dataset={dataset}"


def get_gnomad_sv_overlap(data: dict, gnomad_sv_version: str) -> list[dict]:
    """ data contains the GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS (VEP '&' separated multi values)
        gnomad_sv_version is VariantAnnotationVersion.gnomad_sv (parsed from the annotation settings
        data file by vep_config.vep_component_version_kwargs) """
    dataset = None
    remove_prefix = None
    gnomad_sv_version = gnomad_sv_version or ""
    if gnomad_sv_version.startswith("2.1"):
        dataset = "gnomad_sv_r2_1"
        remove_prefix = "gnomAD-SV_v2.1_"
    elif gnomad_sv_version.startswith("4"):
        dataset = "gnomad_sv_r4"
        # gnomAD v4 ships the v3 SV call set, whose records keep their v3 names
        remove_prefix = "gnomAD-SV_v3_"

    sv_overlap_list = []
    for sv_overlap in vep_multi_fields_to_list_of_dicts(
            data, GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS):
        sv_overlap['gnomad_sv_overlap_af'] = float(sv_overlap['gnomad_sv_overlap_af'])

        gnomad_variant = sv_overlap["gnomad_sv_overlap_name"]
        if remove_prefix:
            gnomad_variant = gnomad_variant.replace(remove_prefix, "")

        sv_overlap['gnomad_sv_overlap_url'] = get_gnomad_url(gnomad_variant, dataset)
        if coords := sv_overlap.get("gnomad_sv_overlap_coords"):
            _chrom, start, end = parse_gnomad_coord(coords)
            sv_overlap['gnomad_sv_overlap_length'] = end - start

        sv_overlap_list.append(sv_overlap)
    return sv_overlap_list
