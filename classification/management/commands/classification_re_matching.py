from enum import StrEnum
from time import sleep
from django.core.management import BaseCommand
from classification.classification_import import reattempt_variant_matching
from classification.models import Classification, ImportedAlleleInfo, DiscordanceReport, ClinVarExport, \
    ResolvedVariantInfo
from library.guardian_utils import admin_bot
import pandas as pd
from snpdb.models import GenomeBuildPatchVersion


class RematchLevel(StrEnum):
    SOFT = "SOFT"
    HARD = "HARD"


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--link_unlinked', action='store_true',
                            help="Add 'clingen_allele_id' to Classifications missing it")
        parser.add_argument('--file', type=str,
                            help="A csv file with columns identifying which Allele Infos to rematch")
        parser.add_argument('--file_commit', action='store_true',
                            help="provide if you")

    def handle(self, *args, **options):
        if options["link_unlinked"]:
            self.link_all_unlinked()
        if file := options["file"]:
            self.rematch_file(file, commit=options["file_commit"])

    GENOME_BUILD_COL = "Imported Genome Build"
    C_HGVS_COL = "c.HGVS (Imported)"
    REMATCH_LEVEL = "Rematch"

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

    def rematch_file(self, file_name: str, commit: bool = False):
        df = pd.read_csv(file_name, sep=",", low_memory=False)
        for idx, row in df.iterrows():
            genome_build_patch_ver: GenomeBuildPatchVersion
            genome_build_str = row[Command.GENOME_BUILD_COL]
            c_hgvs = row[Command.C_HGVS_COL]

            if "." in genome_build_str:
                genome_build_main, genome_patch_version = genome_build_str.split(".")
                genome_patch_version = genome_patch_version[1:]  # drop the leading p
                genome_build_patch_ver = GenomeBuildPatchVersion.objects.filter(
                    genome_build__name=genome_build_main,
                    patch_version=genome_patch_version
                ).first()
            else:
                genome_build_patch_ver = GenomeBuildPatchVersion.objects.filter(
                    genome_build=genome_build_str,
                    patch_version__isnull=True
                ).first()

            if not genome_build_patch_ver:
                raise ValueError(f"Could not find Genome Patch Version {genome_build_str}")

            rematch_level = RematchLevel(row[Command.REMATCH_LEVEL])
            # imported_c_hgvs = TextField(null=True, blank=True)
            imported_allele_info = ImportedAlleleInfo.objects.filter(
                imported_c_hgvs=c_hgvs,
                imported_genome_build_patch_version=genome_build_patch_ver
            ).first()

            if not imported_allele_info:
                raise ValueError(f"Could not find Imported Allele Info {genome_build_str} {c_hgvs}")

            if rematch_level == RematchLevel.HARD:
                if allele := imported_allele_info.allele:
                    if DiscordanceReport.objects.filter(clinical_context__allele=allele).exists():
                        raise ValueError(f"Can't rematch {imported_allele_info} hard as existing allele has Discordance Reports")
                    if ClinVarExport.objects.filter(clinvar_allele__allele=allele).exists():
                        raise ValueError(f"Can't rematch {imported_allele_info} hard as existing allele has ClinVar Export")
                    print(f"{genome_build_patch_ver}\t{c_hgvs}\tFound - Safe for hard rematch")
            else:
                print(f"{genome_build_patch_ver}\t{c_hgvs}\tFound")

            if commit:
                if rematch_level == RematchLevel.HARD:
                    allele = imported_allele_info.allele
                    ResolvedVariantInfo.objects.filter(allele_info=imported_allele_info).delete()
                    Classification.objects.filter(allele=allele).update(allele=None)
                    allele.delete()

                    reattempt_variant_matching(admin_bot(), ImportedAlleleInfo.objects.filter(pk=imported_allele_info.pk), True)
                else:
                    reattempt_variant_matching(admin_bot(), ImportedAlleleInfo.objects.filter(pk=imported_allele_info.pk), True)
                sleep(10)
                print("Rematched")

        print("Completed")

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
