import logging

from django.db import transaction

from analysis.analysis_templates import auto_launch_analysis_templates_for_sample
from analysis.tasks.variant_tag_tasks import variant_tag_created_task, variant_tag_deleted_in_analysis_task, \
    analysis_tag_nodes_set_dirty
from snpdb.models import ImportStatus


def variant_tag_create(sender, instance, created=False, **kwargs):
    if created:
        if instance.analysis:
            analysis_tag_nodes_set_dirty(instance.analysis, instance.tag, visible=True)
        # want to be as quick as we can so do analysis reload + liftover async
        # Need to launch this at end of transaction so we know VariantTag is in DB for celery job
        celery_task = variant_tag_created_task.si(instance.pk)
        transaction.on_commit(lambda: celery_task.apply_async())


def variant_tag_delete(sender, instance, **kwargs):
    if instance.analysis:
        analysis_tag_nodes_set_dirty(instance.analysis, instance.tag, visible=True)
        # want to be as quick as we can so do analysis reload + liftover async
        celery_task = variant_tag_deleted_in_analysis_task.si(instance.analysis_id, instance.tag_id)
        transaction.on_commit(lambda: celery_task.apply_async())


def handle_vcf_import_success(*args, **kwargs):
    vcf = kwargs["vcf"]

    for sample in vcf.sample_set.all():
        auto_launch_analysis_templates_for_sample(vcf.user, sample,
                                                  analysis_description="Auto Created from vcf_import_success signal",
                                                  skip_already_analysed=True)


def handle_active_sample_gene_list_created(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    # At the moment, VCFs are sent up by API or found via sequencing scan BEFORE QCGeneLists
    # which become SampleGeneList/ActiveSampleGeneList
    # So the 1st time we called auto_launch it would have skipped the templates that requires_sample_gene_list
    # As they would have failed. Now we have them, try again to now run previously skipped

    if created:
        sample = instance.sample
        if sample.import_status == ImportStatus.SUCCESS:
            user = sample.vcf.user
            auto_launch_analysis_templates_for_sample(user, sample,
                                                      analysis_description="Auto Created from sample_gene_list_created signal",
                                                      skip_already_analysed=True)
        else:
            logging.warning("Skipping auto analysis for sample: %s, import_status=%s",
                            sample, sample.get_import_status_display())

