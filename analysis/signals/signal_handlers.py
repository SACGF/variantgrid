from django.db import transaction
from django.db.models import Q

from analysis.analysis_templates import auto_launch_analysis_templates_for_sample
from analysis.models import TagNode, Analysis, Tag
from analysis.tasks.variant_tag_tasks import analysis_tag_created_task, analysis_tag_deleted_task
from library.guardian_utils import admin_bot


def _analysis_tag_nodes_set_dirty(analysis: Analysis, tag: Tag):
    """ Needs to be sync so version is bumped by the time client calls node_versions to see whats dirty """
    tag_filter = Q(tagnodetag__tag__isnull=True) | Q(tagnodetag__tag=tag)
    for node in TagNode.objects.filter(analysis=analysis).filter(tag_filter).distinct():
        node.queryset_dirty = True
        node.save()


def variant_tag_create(sender, instance, created=False, **kwargs):
    if created:
        if instance.analysis:
            _analysis_tag_nodes_set_dirty(instance.analysis, instance.tag)
        # want to be as quick as we can so do analysis reload + liftover async
        # Need to launch this at end of transaction so we know VariantTag is in DB for celery job
        celery_task = analysis_tag_created_task.si(instance.pk)
        transaction.on_commit(lambda: celery_task.apply_async())


def variant_tag_delete(sender, instance, **kwargs):
    if instance.analysis:
        _analysis_tag_nodes_set_dirty(instance.analysis, instance.tag)
        # want to be as quick as we can so do analysis reload + liftover async
        celery_task = analysis_tag_deleted_task.si(instance.analysis_id, instance.tag_id)
        transaction.on_commit(lambda: celery_task.apply_async())


def handle_vcf_import_success(*args, **kwargs):
    vcf = kwargs["vcf"]

    for sample in vcf.sample_set.all():
        auto_launch_analysis_templates_for_sample(vcf.user, sample)


def handle_active_sample_gene_list_created(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    # At the moment, VCFs are sent up by API or found via sequencing scan BEFORE QCGeneLists
    # which become SampleGeneList/ActiveSampleGeneList
    # So the 1st time we called auto_launch it would have skipped the templates that requires_sample_gene_list
    # As they would have failed. Now we have them, try again to now run previously skipped

    if created:
        sample = instance.sample
        user = sample.vcf.user
        auto_launch_analysis_templates_for_sample(user, sample, skip_already_analysed=True)
