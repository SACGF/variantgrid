from typing import Dict

from django.conf import settings
from django.contrib.auth.models import User

from classification.models import ClassificationImport, ImportedAlleleInfo, ImportedAlleleInfoStatus
from classification.tasks.classification_import_task import process_classification_import_task
from snpdb.models import ImportSource, GenomeBuild


class VariantResolver:

    def __init__(self, user: User, limit: int = 100):
        self.user = user
        self.ci_by_gb: Dict[GenomeBuild, ClassificationImport] = {}
        self.pending_count = 0
        self.all_count = 0
        self.limit = limit

    def _classification_import_from(self, genome_build: GenomeBuild):
        if existing := self.ci_by_gb.get(genome_build):
            return existing
        new_import = ClassificationImport.objects.create(genome_build=genome_build, user=self.user)
        self.ci_by_gb[genome_build] = new_import
        return new_import

    def queue_resolve(self, imported_allele_info: ImportedAlleleInfo) -> bool:
        actually_queued = False
        if settings.VARIANT_CLASSIFICATION_MATCH_VARIANTS:
            if imported_allele_info.status == ImportedAlleleInfoStatus.PROCESSING and imported_allele_info.classification_import is None:
                if genome_build_version := imported_allele_info.imported_genome_build_patch_version:
                    self.all_count += 1
                    self.pending_count += 1
                    imported_allele_info.classification_import = self._classification_import_from(genome_build_version.genome_build)
                    imported_allele_info.save()
                    actually_queued = True

                    if self.pending_count >= self.limit:
                        self.process_queue()
        return actually_queued

    def process_queue(self):
        for vcimport in self.ci_by_gb.values():
            task = process_classification_import_task.si(vcimport.pk, ImportSource.API)
            task.apply_async()

        self.ci_by_gb.clear()
        self.pending_count = 0
        return self.all_count

