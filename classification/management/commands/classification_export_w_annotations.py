from typing import Optional

from django.core.management import BaseCommand
from annotation.models import AnnotationVersion
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.models import Classification, EvidenceKey
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
                "status",
                "splice_ai_acceptor_gain",
                "splice_ai_acceptor_loss",
                "splice_ai_donor_gain",
                "splice_ai_donor_loss"
            ] + transcript_annotation_keys)

            annotation_version = AnnotationVersion.latest(genome_build)
            last_variant: Optional[Variant] = None
            variant_annotations: list
            status = "unresolved"
            for c in Classification.objects.filter(withdrawn=False, allele_info__isnull=False).select_related("lab", "allele_info").order_by('allele_info__allele').iterator():
                variant_annotation_cells = [None] * 4
                transcript_annotations_cells = [None] * len(transcript_annotation_keys)
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

                                if variant_annotation and variant_annotation.has_spliceai():
                                    variant_annotation_cells = [
                                        variant_annotation.spliceai_pred_ds_ag,
                                        variant_annotation.spliceai_pred_ds_al,
                                        variant_annotation.spliceai_pred_ds_dg,
                                        variant_annotation.spliceai_pred_ds_dl
                                    ]

                                annotation_data = vts.get_transcript_annotation(transcript_version)
                                # the above line will throw an exception if the transcript's not annotated
                                status = "transcript-annotated"

                                transcript_annotations_cells = []
                                for key in transcript_annotation_keys:
                                    transcript_annotations_cells.append(annotation_data.get(key))

                        except Exception:
                            pass

                row = [
                    c.cr_lab_id,
                    c.lab.group_name,
                    c.lab_record_id,
                    c.chgvs_grch38,
                    status
                ] + variant_annotation_cells + transcript_annotations_cells

                csv_writer.writerow(row)