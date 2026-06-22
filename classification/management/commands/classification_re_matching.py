from dataclasses import dataclass
from enum import StrEnum
from functools import cached_property
from time import sleep
from django.core.management import BaseCommand
from django.db.models import Q

from classification.classification_import import reattempt_variant_matching
from classification.models import Classification, ImportedAlleleInfo, DiscordanceReport, ClinVarExport, \
    ResolvedVariantInfo
from library.guardian_utils import admin_bot
import pandas as pd
from snpdb.models import GenomeBuildPatchVersion


class RematchLevel(StrEnum):
    SOFT = "SOFT"
    HARD = "HARD"


GENOME_BUILD_COL = "Imported Genome Build"
C_HGVS_COL = "c.HGVS (Imported)"
REMATCH_LEVEL = "Rematch"


@dataclass
class RematchRequest:
    genome_build_str: str
    c_hgvs: str
    rematch_level: RematchLevel

    @staticmethod
    def from_row(row):
        genome_build_str = row[GENOME_BUILD_COL]
        c_hgvs = row[C_HGVS_COL]
        rematch_level = RematchLevel(row[REMATCH_LEVEL])

        return RematchRequest(genome_build_str=genome_build_str, c_hgvs=c_hgvs, rematch_level=rematch_level)

    @cached_property
    def genome_build_patch_ver(self) -> GenomeBuildPatchVersion:
        genome_build_patch_ver: GenomeBuildPatchVersion
        if "." in self.genome_build_str:
            genome_build_main, genome_patch_version = self.genome_build_str.split(".")
            genome_patch_version = genome_patch_version[1:]  # drop the leading p
            genome_build_patch_ver = GenomeBuildPatchVersion.objects.filter(
                genome_build__name=genome_build_main,
                patch_version=genome_patch_version
            ).first()
        else:
            genome_build_patch_ver = GenomeBuildPatchVersion.objects.filter(
                genome_build=self.genome_build_str,
                patch_version__isnull=True
            ).first()
        if not genome_build_patch_ver:
            raise ValueError(f"Could not find Genome Patch Version {self.genome_build_str}")
        return genome_build_patch_ver

    @property
    def imported_allele_info(self) -> ImportedAlleleInfo:
        iai = ImportedAlleleInfo.objects.filter(
            Q(imported_c_hgvs=self.c_hgvs) | Q(imported_g_hgvs=self.c_hgvs),
            imported_genome_build_patch_version=self.genome_build_patch_ver
        ).first()
        if not iai:
            raise ValueError(f"Could not find Imported Allele Info for {self.genome_build_str} {self.c_hgvs}")
        return iai

    def validate(self, nuke: bool = False):
        # make sure imported allele info can be found
        _ = self.imported_allele_info
        if nuke:
            return

        # if hard rematch (where we're deleting an allele)
        if self.rematch_level == RematchLevel.HARD:
            if allele := self.imported_allele_info.allele:
                dr = DiscordanceReport.objects.filter(clinical_context__allele=allele).first()
                ces = list(ClinVarExport.objects.filter(clinvar_allele__allele=allele).all())
                if dr or ces:
                    ces_ids = " ,".join(str(ce.pk) for ce in ces)
                    raise ValueError(f"Can't rematch {self.imported_allele_info} HARD - Discordance Report {dr}, ClinVarExports - {ces_ids}")

    def prep_if_hard(self, nuke: bool = False):
        if self.rematch_level == RematchLevel.HARD:
            if nuke:
                if allele := self.imported_allele_info.allele:
                    DiscordanceReport.objects.filter(clinical_context__allele=allele).delete()
                    ClinVarExport.objects.filter(clinvar_allele__allele=allele).delete()

            ResolvedVariantInfo.objects.filter(allele_info=self.imported_allele_info).delete()
            if allele := self.imported_allele_info.allele:
                Classification.objects.filter(allele=allele).update(allele=None)
                allele.delete()


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--link_unlinked', action='store_true',
                            help="Add 'clingen_allele_id' to Classifications missing it")
        parser.add_argument('--file', type=str,
                            help="A csv file with columns identifying which Allele Infos to rematch")
        parser.add_argument('--file_commit', action='store_true',
                            help="provide if you")
        parser.add_argument("--nuke", action='store_true', help="Will delete ClinVarExports and DiscordanceReports that are linked")

    def handle(self, *args, **options):
        if options["link_unlinked"]:
            self.link_all_unlinked()
        if file := options["file"]:
            self.rematch_file(file, commit=options["file_commit"], nuke=options["nuke"])

    def rematch_file(self, file_name: str, commit: bool = False, nuke: bool = False):
        print(f"Committing changes = {commit}")
        df = pd.read_csv(file_name, sep=",", low_memory=False)

        passed_validation = True
        rematch_requests: list[RematchRequest] = []
        for idx, row in df.iterrows():
            rematch_request = RematchRequest.from_row(row)
            try:
                rematch_request.validate(nuke=nuke)
            except Exception as ex:
                print(ex)
                passed_validation = False
            rematch_requests.append(rematch_request)

        if not passed_validation:
            print("1 or more records failed validation, exiting")
            return

        print("Validation complete")

        if commit:
            for rematch_request in rematch_requests:
                rematch_request.prep_if_hard(nuke=nuke)

            batch: list[int] = []
            for rematch_request in rematch_requests:
                batch.append(rematch_request.imported_allele_info.pk)
                if len(batch) > 50:
                    print("Rematching Batch")
                    reattempt_variant_matching(admin_bot(), ImportedAlleleInfo.objects.filter(pk__in=batch), True)
                    batch.clear()
                    sleep(50)

            if batch:
                print("Rematching Batch")
                reattempt_variant_matching(admin_bot(), ImportedAlleleInfo.objects.filter(pk__in=batch), True)

        print("Completed")

# IF LINKED UNLINKED

    def rematch_classifications(self, _bulk: list[Classification]):
        allele_infos = []
        for classification in _bulk:
            _, created = classification.ensure_allele_info_with_created(force_allele_info_update_check=True)
            if created:
                classification.save()
            if allele_info := classification.allele_info:
                allele_infos.append(allele_info)

        allele_info_qs = ImportedAlleleInfo.objects.filter(pk__in={ai.pk for ai in allele_infos})
        reattempt_variant_matching(admin_bot(), allele_info_qs, False)
        sleep(50)

    def link_all_unlinked(self):
        rematch_count = 0
        total_count = Classification.objects.count()

        print(f"Will rematch classifications with no Allele Info, speed will vary based on how many classifications were imported with the same details")
        print(f"Total classification count = {total_count}")

        bulk: list[Classification] = []
        for cr in Classification.objects.filter(allele_info__isnull=True).iterator():
            bulk.append(cr)
            rematch_count += 1
            if len(bulk) >= 100:
                self.rematch_classifications(bulk)
                bulk.clear()
                print(f"Rematched {rematch_count}")

        if bulk:
            self.rematch_classifications(bulk)

        print(f"Rematched {rematch_count}")
        print("Finished")
