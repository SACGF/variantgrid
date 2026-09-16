"""Report where gene-level classification targets have got to, and why.

A gene-level ImportedAlleleInfo goes through the same VCF insert pipeline as any other coordinate
(@see snpdb.gene_level_variants), so "stuck in Processing" is nearly always a pipeline that died
rather than anything wrong with the record. This gathers the record states and the pipelines behind
them into one place, since chasing it by hand means joining ImportedAlleleInfo to ClassificationImport
to UploadedClassificationImport to UploadPipeline.
"""
from collections import Counter

from django.core.management.base import BaseCommand

from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from upload.models import UploadedClassificationImport, UploadPipeline


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--pipelines', action='store_true',
                            help="Also show the upload pipeline and step status behind each import")
        parser.add_argument('--list', action='store_true',
                            help="List every record rather than just the counts")

    def handle(self, *args, **options):
        qs = ImportedAlleleInfo.objects.filter(variant_coordinate__startswith=f"{GENE_LEVEL_CONTIG_NAME}:")
        total = qs.count()
        self.stdout.write(f"{total} gene-level ImportedAlleleInfo records")
        if not total:
            self.stdout.write("Nothing gene-level has a coordinate - an import that died before resolving one "
                              "shows up in 'classification_rematch_stuck --gene-level --dry-run' instead")
            return

        statuses = Counter(qs.values_list('status', flat=True))
        for status, count in sorted(statuses.items()):
            self.stdout.write(f"  {ImportedAlleleInfoStatus(status).label}: {count}")
        self.stdout.write(f"  with a matched variant: {qs.filter(matched_variant__isnull=False).count()}")
        self.stdout.write(f"  includable in exports: {qs.filter(latest_validation__include=True).count()}")

        tags = Counter()
        for validation_tags in qs.values_list('latest_validation__validation_tags', flat=True):
            for category, sub_dict in (validation_tags or {}).items():
                for field, severity in sub_dict.items():
                    tags[f"{category}.{field} ({severity})"] += 1
        if tags:
            self.stdout.write("validation tags:")
            for tag, count in sorted(tags.items()):
                self.stdout.write(f"  {tag}: {count}")

        if options['list']:
            self.stdout.write("records:")
            for allele_info in qs.order_by('pk'):
                include = allele_info.latest_validation.include if allele_info.latest_validation else None
                self.stdout.write(f"  {allele_info.pk}\t{allele_info.get_status_display()}\t"
                                  f"include={include}\t{allele_info.imported_hgvs}\t{allele_info.message}")

        if options['pipelines']:
            import_ids = sorted(set(qs.values_list('classification_import_id', flat=True)) - {None})
            self.stdout.write(f"pipelines for {len(import_ids)} ClassificationImport(s):")
            with_uploads = set(UploadedClassificationImport.objects
                               .filter(classification_import_id__in=import_ids)
                               .values_list('classification_import_id', flat=True))
            if missing := [pk for pk in import_ids if pk not in with_uploads]:
                # process_classification_import_task is what builds the uploads, and it is dispatched with
                # apply_async - with no worker consuming it the records stay in Processing forever
                self.stdout.write(f"  ClassificationImport {missing} has no upload at all - "
                                  f"process_classification_import_task never ran (is a celery worker up?)")
            for uci in UploadedClassificationImport.objects.filter(classification_import_id__in=import_ids):
                for pipeline in UploadPipeline.objects.filter(file_upload_id=uci.file_upload_id):
                    self.stdout.write(f"  ClassificationImport {uci.classification_import_id} "
                                      f"pipeline {pipeline.pk} {pipeline.get_status_display()} "
                                      f"({pipeline.file_upload.name})")
                    for step in pipeline.uploadstep_set.order_by('sort_order', 'pk'):
                        error = step.error_message.replace("\n", " ")[:200] if step.error_message else ""
                        self.stdout.write(f"    {step.name}\t{step.get_status_display()}\t{error}")
