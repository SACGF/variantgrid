from django.core.management import BaseCommand
from annotation.models import AnnotationVersion
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.models import Classification, EvidenceKey
from snpdb.models import GenomeBuild
import csv

class Command(BaseCommand):

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        with open('annotation_export.tsv', 'w') as annotation_export_out:
            csv_writer = csv.writer(annotation_export_out, delimiter="\t")

            c: Classification
            genome_build = GenomeBuild.grch38()
            evidence_keys_list = list(EvidenceKey.objects.filter(variantgrid_column__isnull=False).select_related("variantgrid_column"))
            evidence_keys_list_names = list(sorted(x.key for x in evidence_keys_list))

            headers = [
                "cr_id",
                "lab",
                "lab_record_id",
                "c_hgvs_38"
            ] + evidence_keys_list_names

            csv_writer.writerow(headers)

            annotation_version = AnnotationVersion.latest(genome_build)
            for c in Classification.objects.filter(withdrawn=False, allele_info__isnull=False).select_related("lab").iterator():
                if build := c.allele_info.grch38:
                    if c_hgvs := build.c_hgvs_obj:
                        transcript = c_hgvs.transcript
                        if transcript.startswith("N"):
                            if annotation_evidence := get_evidence_fields_for_variant(
                                genome_build=genome_build,
                                variant=build.variant,
                                refseq_transcript_accession=transcript,
                                ensembl_transcript_accession=None,
                                evidence_keys_list=evidence_keys_list,
                                annotation_version=annotation_version
                            ):
                                row = [
                                    c.cr_lab_id,
                                    c.lab.group_name,
                                    c.lab_record_id,
                                    c.chgvs_grch38
                                ]
                                if data := annotation_evidence.data:
                                    for evidence_key in evidence_keys_list_names:
                                        value = None
                                        if cell := data.get(evidence_key):
                                            if isinstance(cell, dict):
                                                value = cell.get('value')
                                        row.append(value)

                                csv_writer.writerow(row)