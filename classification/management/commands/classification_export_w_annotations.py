from typing import Optional

from django.core.management import BaseCommand
from annotation.models import AnnotationVersion
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.models import Classification, EvidenceKey
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild, Variant
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

            transcript_annotation_keys: Optional[list[str]] = [
                "amino_acids", "annotation_run", "cadd_phred", "canonical", "canonical_score", "codons", "consequence",
                "distance", "domains", "ensembl_protein", "ensembl_transcript_accession", "exon",
                "fathmm_pred_most_damaging", "flags", "gene", "gene_id", "gene_symbol", "gerp_pp_rs",
                "gnomad_gene_constraint_method", "gnomad_gene_constraint_oe_lof_summary", "gnomad_gene_constraint_url",
                "grantham", "hgvs_c", "hgvs_p", "id", "impact", "interpro_domain", "intron", "maxentscan_alt",
                "maxentscan_diff", "maxentscan_percent_diff_ref", "maxentscan_ref",
                "mutation_assessor_pred_most_damaging",
                "mutation_taster_pred_most_damaging", "nmd_escaping_variant", "polyphen2_hvar_pred_most_damaging",
                "protein_length", "protein_position", "refseq_transcript_accession", "representative", "revel_score",
                "selected", "sift", "splice_region", "symbol", "tags", "transcript", "transcript_id",
                "transcript_version", "variant", "version"
            ]

            csv_writer.writerow([
                "cr_id",
                "lab",
                "lab_record_id",
                "c_hgvs_38",
                "substituted_transcript_version",
                "status",
                "gnomad_af",
                "gnomad_popmax",
                "gnomad_popmax_af",
                "gnomad_url",
                "splice_ai_acceptor_gain",
                "splice_ai_acceptor_loss",
                "splice_ai_donor_gain",
                "splice_ai_donor_loss"
            ] + transcript_annotation_keys)

            annotation_version = AnnotationVersion.latest(genome_build)
            last_variant: Optional[Variant] = None
            variant_annotations: list
            status = "unresolved"
            count = 0
            for c in Classification.objects.filter(withdrawn=False, allele_info__isnull=False).select_related("lab", "allele_info").order_by('allele_info__allele').iterator():
                count += 1
                if count % 100 == 0:
                    print(f"Processed {count}")
                variant_annotation_cells = [None] * 7
                gnomad_cells = [None] * 4
                transcript_annotations_cells = [None] * len(transcript_annotation_keys)
                substituted_transcript_version = None

                if build := c.allele_info.grch38:
                    if build.variant:
                        status = "matched-but-no-annotations"
                        try:
                            if transcript_version := build.transcript_version:
                                if build.variant != last_variant:
                                    vts = None
                                    last_variant = None

                                    # in case the below line causes an exception make sure we don't cache the vts
                                    vts = VariantTranscriptSelections(
                                        variant=build.variant,
                                        genome_build=genome_build,
                                        annotation_version=annotation_version
                                    )
                                    last_variant = build.variant

                                variant_annotation = vts.variant_annotation
                                if variant_annotation:
                                    status = "variant-annotated"

                                    if variant_annotation.has_spliceai():
                                        variant_annotation_cells = [
                                            variant_annotation.spliceai_pred_ds_ag,
                                            variant_annotation.spliceai_pred_ds_al,
                                            variant_annotation.spliceai_pred_ds_dg,
                                            variant_annotation.spliceai_pred_ds_dl,
                                        ]
                                    gnomad_cells = [
                                        variant_annotation.gnomad_af,
                                        variant_annotation.gnomad_popmax,
                                        variant_annotation.gnomad_popmax_af,
                                        variant_annotation.gnomad_url
                                    ]

                                annotation_data: Optional[dict] = None
                                try:
                                    annotation_data = vts.get_transcript_annotation(transcript_version)
                                    # the above line will throw an exception if the transcript's not annotated
                                    status = "transcript-annotated"
                                except:
                                    all_transcript_versions = transcript_version.transcript.transcriptversion_set.filter(genome_build=genome_build)
                                    lower_versions = list(sorted((v for v in all_transcript_versions if v.version < transcript_version.version), reverse=True))
                                    higher_versions = list(sorted((v for v in all_transcript_versions if v.version > transcript_version.version)))
                                    attempt_versions = higher_versions + lower_versions
                                    for attempt_version in attempt_versions:
                                        try:
                                            annotation_data = vts.get_transcript_annotation(attempt_version)
                                            substituted_transcript_version = attempt_version.version
                                            status = "transcript-annotated (sub)"
                                            break
                                        except:
                                            pass

                                if annotation_data:
                                    transcript_annotations_cells = []
                                    for key in transcript_annotation_keys:
                                        transcript_annotations_cells.append(annotation_data.get(key))

                        except Exception as e:
                            pass

                row = [
                    c.cr_lab_id,
                    c.lab.group_name,
                    c.lab_record_id,
                    c.chgvs_grch38,
                    substituted_transcript_version,
                    status
                ] + gnomad_cells + variant_annotation_cells + transcript_annotations_cells

                csv_writer.writerow(row)