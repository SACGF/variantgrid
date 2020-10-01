from collections import Counter, defaultdict
import logging

from celery.app.task import Task

from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from annotation.models.models_gene_counts import GeneCountType, GeneValue, CohortGeneCounts,\
    SampleAnnotationVersionVariantSource, GeneValueCountCollection, GeneValueCount
from annotation.models.molecular_consequence_enums import MolecularConsequenceColors
from eventlog.models import create_event
from genes.models import Gene, GeneVersion
from library.django_utils import thread_safe_unique_together_get_or_create
from library.enums.log_level import LogLevel
from library.log_utils import get_traceback
from patients.models_enums import Zygosity
from snpdb.models.models_enums import ProcessingStatus
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from classification.models import VariantClassification
from variantgrid.celery import app


def get_or_create_gene_count_type_and_values(gene_count_type_name: str):
    """ Make sure a GeneCountType is enabled and has related data defined """
    gene_count_type, _ = GeneCountType.objects.update_or_create(name=gene_count_type_name,
                                                                defaults={"enabled": True})
    gene_values = {}
    for label, rgb in MolecularConsequenceColors.ONCOPLOT_COLORS:
        gv, _ = thread_safe_unique_together_get_or_create(GeneValue,
                                                          gene_count_type=gene_count_type,
                                                          label=label,
                                                          defaults={"rgb": rgb})
        gene_values[gv.label] = gv

    other = gene_values[MolecularConsequenceColors.OTHER]
    other.show_counts = False
    other.save()

    not_tested = gene_values[MolecularConsequenceColors.NOT_TESTED]
    not_tested.use_as_empty_value = True
    not_tested.save()

    return gene_count_type, gene_values


class CalculateCohortSampleGeneDamageCountsTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        uploaded_vcf = upload_step.get_uploaded_vcf()
        vcf = uploaded_vcf.vcf
        cohort = vcf.cohort

        for gene_count_type in GeneCountType.objects.filter(enabled=True):
            variant_annotation_version = VariantAnnotationVersion.latest(vcf.genome_build)
            cgc, created = CohortGeneCounts.objects.get_or_create(variant_annotation_version=variant_annotation_version,
                                                                  gene_count_type=gene_count_type,
                                                                  cohort=cohort,
                                                                  cohort_version=cohort.version)

            if created or (cgc.processing_status not in ProcessingStatus.FINISHED_STATES):
                cgc.launch_task()
            else:
                logging.warning("CohortGeneCounts %s already run", cgc)


class CohortSampleGeneDamageCountTask(Task):
    """ Task name in CohortGeneCounts.celery_task_name, called from CohortGeneCounts.launch_task () """

    def run(self, cohort_gene_counts_id, update_gene_id=None):
        cohort_gene_counts = CohortGeneCounts.objects.get(pk=cohort_gene_counts_id)
        cohort_gene_counts.processing_status = ProcessingStatus.PROCESSING
        cohort_gene_counts.save()

        try:
            gene_count_type, gene_values = get_or_create_gene_count_type_and_values(cohort_gene_counts.gene_count_type.pk)
            variant_annotation_version = cohort_gene_counts.variant_annotation_version
            sample_gene_value_counts = self._get_gene_value_counts(cohort_gene_counts, gene_values,
                                                                   variant_annotation_version, update_gene_id)
            self._save_counts(gene_count_type, gene_values, sample_gene_value_counts, variant_annotation_version, update_gene_id)
            cohort_gene_counts.processing_status = ProcessingStatus.SUCCESS
        except:
            details = get_traceback()
            logging.info(details)
            create_event(None, "cohort_sample_gene_damage_counts", details, severity=LogLevel.ERROR)
            cohort_gene_counts.processing_status = ProcessingStatus.ERROR

        cohort_gene_counts.save()

    @classmethod
    def _get_gene_value_counts(cls, cohort_gene_counts, gene_values, variant_annotation_version, update_gene_id):
        cohort = cohort_gene_counts.cohort
        cohort_genotype_collection = cohort.cohort_genotype_collection
        qs = cohort_gene_counts.gene_count_type.get_variant_queryset(variant_annotation_version)
        qs = qs.filter(cohortgenotype__collection=cohort_genotype_collection)
        if update_gene_id:
            gene = Gene.objects.get(pk=update_gene_id)
            qs = qs.filter(**{VariantAnnotation.GENE_COLUMN: gene})
        values_list = qs.values_list(VariantAnnotation.GENE_COLUMN, "variantannotation__consequence",
                                     "cohortgenotype__samples_zygosity")
        mv = gene_values[MolecularConsequenceColors.MISSENSE]
        sg = gene_values[MolecularConsequenceColors.STOPGAIN]
        fv = gene_values[MolecularConsequenceColors.FRAMESHIFT]
        ss = gene_values[MolecularConsequenceColors.SPLICE_SITE]
        ot = gene_values[MolecularConsequenceColors.OTHER]
        sample_ids = cohort.get_sample_ids()

        sample_gene_value_counts = defaultdict(lambda: defaultdict(Counter))
        for (gene_id, consequence, samples_zygosity) in values_list:
            for zyg, sample_id in zip(samples_zygosity, sample_ids):
                if zyg in Zygosity.VARIANT:
                    g_value_counts = sample_gene_value_counts[sample_id][gene_id]

                    has_category = False
                    if 'missense_variant' in consequence:
                        g_value_counts[mv] += 1
                        has_category = True
                    if 'stop_gained' in consequence:
                        g_value_counts[sg] += 1
                        has_category = True
                    if 'frameshift_variant' in consequence:
                        g_value_counts[fv] += 1
                        has_category = True
                    if 'splice_acceptor_variant' in consequence or 'splice_donor_variant' in consequence:
                        g_value_counts[ss] += 1
                        has_category = True

                    if not has_category:
                        g_value_counts[ot] += 1
        return sample_gene_value_counts

    def _save_counts(self, gene_count_type, gene_values, sample_gene_value_counts, variant_annotation_version, update_gene_id):
        variant_type_rank = {gene_values[vt]: i for i, vt in enumerate(MolecularConsequenceColors.VARIANT_TYPES)}
        for sample_id, g_value_counts in sample_gene_value_counts.items():
            sample_source, _ = SampleAnnotationVersionVariantSource.objects.get_or_create(sample_id=sample_id,
                                                                                          variant_annotation_version=variant_annotation_version)
            collection, created = GeneValueCountCollection.objects.get_or_create(source=sample_source,
                                                                                 gene_count_type=gene_count_type)
            if not (created or update_gene_id):
                collection.genevaluecount_set.all().delete()

            gene_value_counts = []
            for gene_id, value_counts in g_value_counts.items():
                if value_counts:
                    value_counts = sorted(value_counts.items(), key=lambda vc: variant_type_rank[vc[0]])
                    value, count = value_counts[0]
                    gene_value_counts.append(GeneValueCount(collection=collection,
                                                            gene_id=gene_id,
                                                            value=value,
                                                            count=count))

            gene_value_counts = GeneValueCount.objects.bulk_create(gene_value_counts, ignore_conflicts=True)
            if update_gene_id:
                # Inserts to existing fail due to DB unique_together constraints - update individually
                for gvc in gene_value_counts:
                    if gvc.pk is None:
                        GeneValueCount.objects.filter(collection=gvc.collection,
                                                      gene=gvc.gene,
                                                      value=gvc.value).update(count=gvc.count)


class CohortSampleClassificationGeneDamageCountTask(CohortSampleGeneDamageCountTask):
    @classmethod
    def _get_gene_value_counts(cls, cohort_gene_counts, gene_values, variant_annotation_version, update_gene_id):
        # There are way less VariantClassifications, and each has a sample, so handle this per-sample
        # As it's classified for this sample, we know the builds are the same (hence can use variant not allele)
        vc_qs = cohort_gene_counts.gene_count_type.get_variant_classification_qs()
        if update_gene_id:
            gene = Gene.objects.get(pk=update_gene_id)
            vc_qs = vc_qs.filter(VariantClassification.get_q_for_gene(gene))

        sample_gene_value_counts = defaultdict(lambda: defaultdict(Counter))

        cohort = cohort_gene_counts.cohort
        for sample in cohort.get_samples():
            for variant_classification in vc_qs.filter(sample=sample):
                ensembl_gene_id = variant_classification.get("ensembl_gene_id")
                gene_id = None
                if ensembl_gene_id:
                    gene_id, _ = GeneVersion.get_gene_id_and_version(ensembl_gene_id)
                else:
                    variant_annotation = variant_classification.get_variant_annotation(variant_annotation_version)
                    if variant_annotation:
                        gene_id = variant_annotation.gene_id

                if gene_id:
                    g_value_counts = sample_gene_value_counts[sample.pk][gene_id]
                    has_category = False
                    for molecular_consequence in variant_classification.get("molecular_consequence", []):
                        category = MolecularConsequenceColors.CONSEQUENCE_LOOKUPS.get(molecular_consequence)
                        if category:
                            gvc = gene_values[category]
                            g_value_counts[gvc] += 1
                            has_category = True

                    if not has_category:
                        ot = gene_values[MolecularConsequenceColors.OTHER]
                        g_value_counts[ot] += 1
        return sample_gene_value_counts


CalculateCohortSampleGeneDamageCountsTask = app.register_task(CalculateCohortSampleGeneDamageCountsTask())
CohortSampleGeneDamageCountTask = app.register_task(CohortSampleGeneDamageCountTask())
CohortSampleClassificationGeneDamageCountTask = app.register_task(CohortSampleClassificationGeneDamageCountTask())
