import logging
import os
import shutil
from collections import defaultdict

from django.conf import settings
from django.core.exceptions import (
    MultipleObjectsReturned,
    ObjectDoesNotExist,
)
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.signals import pre_delete
from django.dispatch import receiver
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel

from genes.gene_coverage import load_gene_coverage_df
from genes.models.models_gene import GeneSymbol, Transcript, TranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils.data_archive_mixin import DataArchiveMixin
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.log_utils import log_traceback
from library.utils.file_utils import mk_path
from snpdb.archive import DataArchivedError
from snpdb.models import DataState
from snpdb.models.models_genome import GenomeBuild
from upload.vcf.sql_copy_files import (
    GENE_COVERAGE_HEADER,
    gene_coverage_canonical_transcript_sql_copy_csv,
    gene_coverage_sql_copy_csv,
    write_sql_copy_csv,
)


class CanonicalTranscriptCollection(TimeStampedModel):
    description = models.TextField(blank=True)
    filename = models.TextField(blank=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    file_sha256sum = models.TextField()

    @staticmethod
    def get_default():
        ctc_id = settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID
        ctc = None
        if ctc_id:
            try:
                ctc = CanonicalTranscriptCollection.objects.get(pk=ctc_id)
            except CanonicalTranscriptCollection.DoesNotExist:
                logging.error("setting.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID=%s - record not found", ctc_id)
        return ctc

    def get_absolute_url(self):
        return reverse('view_canonical_transcript_collection', kwargs={"pk": self.pk})

    def __str__(self):
        count = self.canonicaltranscript_set.count()
        name = f"{count} transcripts"
        description = self.description or os.path.basename(self.filename)
        if description:
            name = f"'{description}' ({name})"
        return name


class CanonicalTranscript(models.Model):
    collection = models.ForeignKey(CanonicalTranscriptCollection, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, blank=True, on_delete=SET_NULL)
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, blank=True, on_delete=SET_NULL)
    original_gene_symbol = models.TextField()
    original_transcript = models.TextField()


class GeneCoverageCollection(DataArchiveMixin, RelatedModelsPartitionModel):
    """ Note both GeneCoverage and GeneCoverageCanonicalTranscript point off same collection object """
    RECORDS_BASE_TABLE_NAMES = ["genes_genecoverage", "genes_genecoveragecanonicaltranscript"]
    RECORDS_FK_FIELD_TO_THIS_MODEL = "gene_coverage_collection_id"
    PARTITION_LABEL_TEXT = "collection"

    path = models.TextField()
    data_state = models.CharField(max_length=1, choices=DataState.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    @staticmethod
    def get_gene_coverage_for_sample(sample):
        """ returns coverage (whether via seqauto or standalone uploaded coverage """

        # TODO: Try and get from UploadedGeneCoverage

        gene_coverage = None
        try:
            sequencing_sample = sample.samplefromsequencingsample.sequencing_sample
            bam_file = sequencing_sample.get_single_bam()
            try:
                qc = bam_file.qc_set.get()
                gene_coverage = qc.qcgenecoverage.gene_coverage_collection
            except (ObjectDoesNotExist, MultipleObjectsReturned):
                logging.error("Wasn't exactly 1 qc for bam_file %s", bam_file)
        except ObjectDoesNotExist:
            pass
        except Exception:
            log_traceback()

        return gene_coverage

    def load_from_file(self, enrichment_kit, **kwargs):
        logging.debug("GeneCoverageCollection.load_for_qc()")
        try:
            gene_matcher = kwargs["gene_matcher"]
            canonical_transcript_manager = kwargs.get("canonical_transcript_manager")
            transcript_versions_by_id = kwargs.get("transcript_versions_by_id")
        except KeyError as ke:
            missing_key = ke.args[0]
            logging.error("You need to pass '%s' to kwargs", missing_key)
            raise ke

        missing_genes = 0
        missing_transcripts = 0

        canonical_transcript_collection = canonical_transcript_manager.get_canonical_collection_for_enrichment_kit(enrichment_kit)
        _, original_canonical_transcript_accessions = canonical_transcript_manager.get_canonical_transcripts(canonical_transcript_collection)

        gene_coverage_df = load_gene_coverage_df(self.path)

        gene_coverage_tuples = []
        gene_coverage_canonical_transcript_tuples = []
        warnings = []
        genes_with_matched_canonical = set()
        unmatched_transcripts_by_gene = defaultdict(set)

        for _, row in gene_coverage_df.iterrows():
            original_gene_symbol = row["original_gene_symbol"]
            original_transcript = row["original_transcript"]

            gene_symbol_id = gene_matcher.get_gene_symbol_id(original_gene_symbol)
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(original_transcript)
            if transcript_versions := transcript_versions_by_id.get(transcript_id):
                transcript_version_id = transcript_versions.get(version)
            else:
                transcript_id = None
                transcript_version_id = None

            gc_dict = {
                "gene_coverage_collection_id": self.pk,
                "transcript_id": transcript_id,
                "transcript_version_id": transcript_version_id,
                "gene_symbol_id": gene_symbol_id,
            }
            gc_dict.update(row.to_dict())
            gc_tup = tuple(gc_dict.get(f, '') for f in GENE_COVERAGE_HEADER)

            if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_ALL:
                gene_coverage_tuples.append(gc_tup)

            if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL:
                if original_transcript in original_canonical_transcript_accessions:
                    gene_coverage_canonical_transcript_tuples.append((canonical_transcript_collection.pk,) + gc_tup)
                    genes_with_matched_canonical.add(gene_symbol_id)
                else:
                    unmatched_transcripts_by_gene[gene_symbol_id].add(original_transcript)

        if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL:
            genes_with_unmatched_transcripts = set(unmatched_transcripts_by_gene.keys())
            genes_no_match = genes_with_unmatched_transcripts - genes_with_matched_canonical
            num_unmatched = len(genes_no_match)
            num_matched = len(genes_with_matched_canonical)
            if num_unmatched > num_matched:
                def get_examples(iterable):
                    return ", ".join(list(iterable)[:5])

                first_transcripts_iter = (next(iter(v)) for v in unmatched_transcripts_by_gene.values())
                sample_unmatched_transcripts = get_examples(first_transcripts_iter)
                sample_canonical_transcripts = get_examples(original_canonical_transcript_accessions)

                message = (f"More genes with NO matched transcripts ({num_unmatched}) than matched ({num_matched})."
                           f"{enrichment_kit=} {canonical_transcript_collection=} (may be fall back if not set on kit) "
                           f"Sample unmatched transcripts: {sample_unmatched_transcripts}. "
                           f"Sample transcripts from canonical transcript collection: {sample_canonical_transcripts}. ")
                raise ValueError(message)

        processing_dir = os.path.join(settings.IMPORT_PROCESSING_DIR, "gene_coverage",
                                      f"gene_coverage_collection_{self.pk}")
        mk_path(processing_dir)
        if gene_coverage_tuples:
            csv_filename = os.path.join(processing_dir, f"gene_coverage_{self.pk}.csv")
            write_sql_copy_csv(gene_coverage_tuples, csv_filename)
            partition_table = self.get_partition_table(base_table_name="genes_genecoverage")
            gene_coverage_sql_copy_csv(csv_filename, partition_table)

        if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL:
            if gene_coverage_canonical_transcript_tuples:
                csv_filename = os.path.join(processing_dir, f"gene_coverage_canonical_transcript_{self.pk}.csv")
                write_sql_copy_csv(gene_coverage_canonical_transcript_tuples, csv_filename)
                partition_table = self.get_partition_table(base_table_name="genes_genecoveragecanonicaltranscript")
                gene_coverage_canonical_transcript_sql_copy_csv(csv_filename, partition_table)
            else:
                logging.warning("GeneCoverage had no canonical transcripts")

        if not settings.DEBUG:
            shutil.rmtree(processing_dir)

        logging.info("%d missing genes, %d missing transcripts", missing_genes, missing_transcripts)
        return warnings

    def get_uncovered_gene_symbols(self, gene_symbols, min_coverage):
        if self.data_archived:
            raise DataArchivedError(self)
        # Do as inner query to ensure we restrict to gene coverage partition
        covered_qs = GeneCoverageCanonicalTranscript.objects.filter(gene_coverage_collection=self,
                                                                    min__gte=min_coverage)
        covered_gene_symbols = covered_qs.values_list("gene_symbol", flat=True)
        return gene_symbols.exclude(pk__in=covered_gene_symbols)

    def __str__(self):
        return "GeneCoverageCollection " + os.path.basename(self.path)


@receiver(pre_delete, sender=GeneCoverageCollection)
def gene_coverage_collection_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    try:
        instance.delete_related_objects()
    except Exception:
        pass


class AbstractGeneCoverage(models.Model):
    gene_coverage_collection = models.ForeignKey(GeneCoverageCollection, on_delete=CASCADE)  # rename to "coverage"?
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, blank=True, on_delete=SET_NULL)
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, blank=True, on_delete=SET_NULL)
    original_gene_symbol = models.TextField()  # raw input string
    original_transcript = models.TextField()  # raw input string
    size = models.IntegerField(null=True)
    min = models.IntegerField()
    max = models.IntegerField(null=True)
    mean = models.FloatField()
    std_dev = models.FloatField()
    # As per QCExecSummary - different enrichment kits store different coverages
    # So some will be None
    percent_1x = models.FloatField(null=True)
    percent_2x = models.FloatField(null=True)
    percent_5x = models.FloatField(null=True)
    percent_10x = models.FloatField(null=True)
    percent_15x = models.FloatField(null=True)
    percent_20x = models.FloatField(null=True)
    percent_25x = models.FloatField(null=True)
    percent_30x = models.FloatField(null=True)
    percent_40x = models.FloatField(null=True)
    percent_50x = models.FloatField(null=True)
    percent_60x = models.FloatField(null=True)
    percent_80x = models.FloatField(null=True)
    percent_100x = models.FloatField(null=True)
    percent_150x = models.FloatField(null=True)
    percent_200x = models.FloatField(null=True)
    percent_250x = models.FloatField(null=True)

    class Meta:
        abstract = True

    @classmethod
    def get_for_symbol(cls, genome_build, gene_symbol):
        return cls.objects.filter(gene_symbol=gene_symbol, gene_coverage_collection__genome_build=genome_build)


class GeneCoverage(AbstractGeneCoverage):
    """ We attempt to match transcript from refseq_transcript_id then if that fails,
        match on gene. So it's possible to have a gene match but no transcript """


class GeneCoverageCanonicalTranscript(AbstractGeneCoverage):
    """ A single transcript per gene to be used for some metrics """

    canonical_transcript_collection = models.ForeignKey(CanonicalTranscriptCollection, null=True, on_delete=SET_NULL)

    @staticmethod
    def filter_for_kit_and_gene_symbol(enrichment_kit, genome_build, gene_symbol):
        sequencing_sample = "gene_coverage_collection__qcgenecoverage__qc__bam_file__sequencing_sample"
        kwargs = {sequencing_sample + "__enrichment_kit": enrichment_kit,
                  # Ensure we only get current SampleSheet
                  sequencing_sample + "__sample_sheet__sequencingruncurrentsamplesheet__isnull": False}
        return GeneCoverageCanonicalTranscript.get_for_symbol(genome_build, gene_symbol).filter(**kwargs)
