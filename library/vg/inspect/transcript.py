"""
`vg inspect transcript <accession>`: a Transcript (NM_000059) or one TranscriptVersion (NM_000059.4) -
consortium, the gene, and every version per build with its cdot data state (length, exons, gaps,
HGVS-usable), MANE and canonical membership, LRG mapping.
"""
from typing import Any

from genes.models.models_gene import MANE, LRGRefSeqGene, Transcript, TranscriptVersion
from genes.models.models_gene_coverage import CanonicalTranscript
from library.vg.inspect import capped


def load(key: str) -> tuple[Transcript, int | None]:
    transcript_id, version = TranscriptVersion.get_transcript_id_and_version(key.strip())
    try:
        return Transcript.objects.get(pk=transcript_id), version
    except Transcript.DoesNotExist as e:
        raise LookupError(f"No Transcript {transcript_id!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    transcript, requested_version = load(key)
    versions = TranscriptVersion.objects.filter(transcript=transcript).select_related("genome_build", "gene_version__gene_symbol", "gene_version__gene")
    if requested_version is not None:
        versions = versions.filter(version=requested_version)
    versions = list(versions.order_by("version", "genome_build__name"))
    data: dict[str, Any] = {
        "id": transcript.pk,
        "consortium": transcript.get_annotation_consortium_display(),
        "gene": _gene(versions),
        "url": transcript.get_absolute_url(),
        "versions": [_version(tv) for tv in versions],
    }
    if depth >= 2:
        data["mane"] = capped(MANE.objects.filter(refseq_transcript_version__transcript=transcript) | MANE.objects.filter(ensembl_transcript_version__transcript=transcript),
                              lambda m: {"symbol": m.symbol_id, "status": m.get_status_display() if hasattr(m, "get_status_display") else m.status,
                                         "refseq": str(m.refseq_transcript_version), "ensembl": str(m.ensembl_transcript_version)})
        data["canonical_in"] = capped(CanonicalTranscript.objects.filter(transcript=transcript).select_related("collection"),
                                      lambda ct: ct.collection.description or ct.collection.filename)
        data["lrg"] = capped(LRGRefSeqGene.objects.filter(rna__startswith=transcript.pk), lambda l: f"{l.lrg}{l.t or ''}")
    return data


def _gene(versions: list[TranscriptVersion]) -> dict[str, Any] | None:
    if not versions:
        return None
    gene_version = versions[-1].gene_version
    return {"symbol": gene_version.gene_symbol_id, "gene": gene_version.gene_id}


def _version(tv: TranscriptVersion) -> dict[str, Any]:
    info: dict[str, Any] = {"accession": tv.accession, "build": tv.genome_build_id, "gene_symbol": tv.gene_version.gene_symbol_id,
                            "biotype": tv.biotype, "has_data": bool(tv.data), "hgvs_ok": tv.hgvs_ok}
    if tv.data:
        info["length"] = tv.length
        info["exons"] = len(tv.data.get("genome_builds", {}).get(tv.genome_build_id, {}).get("exons", [])) or None
        info["alignment_gap"] = tv.alignment_gap
        info["canonical"] = tv.canonical_tag or None
    return info
