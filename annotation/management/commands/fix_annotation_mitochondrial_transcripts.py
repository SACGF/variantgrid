"""
Backfill RefSeq mitochondrial annotation rows written before #2139.

VEP's RefSeq cache names MT transcripts after the gene ('ND4.1'), which never matched a TranscriptVersion, so every
MT row has a NULL transcript_id and a c.HGVS against that non-accession. This links the 13 coding genes to cdot's
'fake-rna-<gene>' TranscriptVersions (the MT row's symbol is VEP's Feature) and swaps hgvs_c for the variant's
'NC_012920.1:m.' HGVS, as annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py now does on insert.
Only MT rows with no transcript are touched. Needs the fake transcripts imported first (import_cdot_latest).
"""
from django.core.management.base import BaseCommand

from annotation.models.models import (
    VariantAnnotation,
    VariantAnnotationVersion,
    VariantTranscriptAnnotation,
)
from genes.models import TranscriptVersion
from genes.models_enums import AnnotationConsortium
from genes.transcript_parts import CDOT_FAKE_TRANSCRIPT_PREFIX, CDOT_FAKE_TRANSCRIPT_VERSION
from snpdb.models.models_enums import AssemblyMoleculeType

BATCH_SIZE = 2000
HGVS_PLACEHOLDERS = {VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE, VariantAnnotation.SV_HGVS_ERROR_MESSAGE}


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        vav_qs = VariantAnnotationVersion.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ)
        for vav in vav_qs.order_by("pk"):
            if vav.data_archived:
                print(f"Skipping archived {vav}")
                continue
            fix_variant_annotation_version(vav)


def _get_fake_transcript_version_ids_by_symbol(vav: VariantAnnotationVersion) -> dict[str, tuple[str, int]]:
    """ {'ND4': ('fake-rna-ND4', TranscriptVersion.pk)} """
    tv_qs = TranscriptVersion.objects.filter(genome_build=vav.genome_build,
                                             transcript__annotation_consortium=AnnotationConsortium.REFSEQ,
                                             transcript__identifier__startswith=CDOT_FAKE_TRANSCRIPT_PREFIX,
                                             version=CDOT_FAKE_TRANSCRIPT_VERSION)
    return {transcript_id.removeprefix(CDOT_FAKE_TRANSCRIPT_PREFIX): (transcript_id, pk)
            for pk, transcript_id in tv_qs.values_list("pk", "transcript_id")}


def fix_variant_annotation_version(vav: VariantAnnotationVersion):
    mito_contig = vav.genome_build.contigs.filter(molecule_type=AssemblyMoleculeType.MITOCHONDRION).first()
    if mito_contig is None:
        return

    fake_transcripts_by_symbol = _get_fake_transcript_version_ids_by_symbol(vav)
    if not fake_transcripts_by_symbol:
        print(f"{vav}: no cdot fake-rna transcripts - run import_cdot_latest first")

    va_qs = VariantAnnotation.objects.filter(version=vav, variant__locus__contig=mito_contig)
    mitochondrial_hgvs_by_variant = {}
    for variant_id, hgvs_g in va_qs.values_list("variant_id", "hgvs_g"):
        if hgvs_g not in HGVS_PLACEHOLDERS:
            mitochondrial_hgvs_by_variant[variant_id] = hgvs_g

    fields = ["transcript_id", "transcript_version_id", "hgvs_c"]
    for klass in [VariantAnnotation, VariantTranscriptAnnotation]:
        unlinked_qs = klass.objects.filter(version=vav, variant__locus__contig=mito_contig, transcript__isnull=True)
        num_linked = num_hgvs_c = 0
        records = []
        for pk, variant_id, symbol, hgvs_c in unlinked_qs.values_list("pk", "variant_id", "symbol", "hgvs_c"):
            transcript_id, transcript_version_id = fake_transcripts_by_symbol.get(symbol, (None, None))
            if transcript_id:
                num_linked += 1
            if hgvs_c and hgvs_c not in HGVS_PLACEHOLDERS:
                hgvs_c = mitochondrial_hgvs_by_variant.get(variant_id)
                num_hgvs_c += 1
            records.append(klass(pk=pk, transcript_id=transcript_id, transcript_version_id=transcript_version_id,
                                 hgvs_c=hgvs_c))
            if len(records) >= BATCH_SIZE:
                klass.objects.bulk_update(records, fields=fields)
                records = []
        if records:
            klass.objects.bulk_update(records, fields=fields)
        print(f"{vav} {klass.__name__}: linked {num_linked} MT rows to a transcript, set m.HGVS on {num_hgvs_c}")
