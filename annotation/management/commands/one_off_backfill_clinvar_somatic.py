from collections import Counter

from django.core.management import BaseCommand
from django.db.models import Q

from annotation.models.models import ClinVar, ClinVarVersion
from annotation.vcf_files.clinvar_significance import (
    NO_CLASSIFICATION_FOR_THE_SINGLE_VARIANT,
    highest_oncogenicity,
    somatic_tier,
)
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        ClinVar's somatic clinical impact (SCI) and oncogenicity (ONC) text has been stored since
        annotation migration 0102, but the queryable somatic_tier / highest_oncogenicity columns are
        only written by imports from 0176 on. Derives them for the versions already loaded.
    """

    def add_arguments(self, parser):
        parser.add_argument('--clinvar-version', type=int,
                            help="ClinVarVersion pk to backfill (default: latest per genome build)")

    def handle(self, *args, **options):
        if clinvar_version_id := options["clinvar_version"]:
            versions = [ClinVarVersion.objects.get(pk=clinvar_version_id)]
        else:
            versions = self._latest_versions_per_build()

        for clinvar_version in versions:
            self._backfill(clinvar_version)

    @staticmethod
    def _latest_versions_per_build() -> list[ClinVarVersion]:
        versions = []
        for genome_build in GenomeBuild.builds_with_annotation():
            qs = ClinVarVersion.objects.filter(genome_build=genome_build).order_by('annotation_date')
            if clinvar_version := qs.last():
                versions.append(clinvar_version)
        return versions

    @staticmethod
    def _is_unrecognised(clinvar_text: str) -> bool:
        """ The haplotype sentinel is expected to derive nothing - anything else that does is a
            ClinVar vocabulary value the maps need to grow """
        return bool(clinvar_text) and clinvar_text != NO_CLASSIFICATION_FOR_THE_SINGLE_VARIANT

    def _backfill(self, clinvar_version: ClinVarVersion):
        qs = ClinVar.objects.filter(version=clinvar_version).filter(
            Q(somatic_clinical_significance__isnull=False) | Q(oncogenic_classification__isnull=False))

        unrecognised = Counter()
        records = []
        selected = 0
        for cv in qs:
            selected += 1
            tier = somatic_tier(cv.somatic_clinical_significance)
            if tier is None and self._is_unrecognised(cv.somatic_clinical_significance):
                unrecognised[f"SCI: {cv.somatic_clinical_significance}"] += 1

            oncogenicity = highest_oncogenicity(cv.oncogenic_classification, cv.oncogenic_conflicting_classification)
            if not oncogenicity and self._is_unrecognised(cv.oncogenic_classification):
                unrecognised[f"ONC: {cv.oncogenic_classification}"] += 1
            oncogenicity = oncogenicity or 0

            if (cv.somatic_tier, cv.highest_oncogenicity) != (tier, oncogenicity):
                cv.somatic_tier = tier
                cv.highest_oncogenicity = oncogenicity
                records.append(cv)

        if records:
            ClinVar.objects.bulk_update(records, fields=["somatic_tier", "highest_oncogenicity"], batch_size=1000)

        print(f"{clinvar_version} ({clinvar_version.genome_build}, "
              f"{clinvar_version.annotation_date}): {selected} rows selected, {len(records)} updated")
        for text, count in unrecognised.most_common():
            print(f"  unrecognised {text} x {count}")
