import time
from typing import Dict

from django.core.management import BaseCommand

from classification.models import Classification, ImportedAlleleInfo
from classification.models.classification_variant_info_models import ResolvedVariantInfo
from flags.models import FlagComment, FlagType
from library.guardian_utils import admin_bot
from library.utils import first
from snpdb.models import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        # parser.add_argument('--all', action='store_true', default=False, help='Attempt to rematch every single classification')
        # parser.add_argument('--file', type=str, help='If provided expects a csv where the first column is a combination of genome build and imported c.hgvs it expects to import')
        # parser.add_argument('--missing', action='store_true', default=False, help='Attempt to rematch only classifications not linked to a variant - one at a time')
        parser.add_argument('--extra', action='store_true', default=False, help='Populate the allele_info of a classification')
        parser.add_argument('--sort', action='store_true', default=False, help='Fixes sort order for non-numerical chromosones')
        parser.add_argument('--revalidate_chgvs', action='store_true', default=False,
                            help='Perform re-validation on any records without outstanding c.hgvs issues')
        parser.add_argument('--validation', action='store_true', default=False, help='Perform validation on the pre-existing normalising/liftovers')
        parser.add_argument('--non-coding', action='store_true', default=False, help='Fix issue #762 NR had c. instead of n.')

    def report_unmatched(self):
        print(f"Unmatched count = {Classification.objects.filter(variant__isnull=True).count()}")

    def handle(self, *args, **options):
        # mode_all = options.get('all')
        # mode_missing = options.get('missing')
        # mode_file = options.get('file')

        mode_extra = options.get('extra')
        mode_validation = options.get('validation')
        mode_sort = options.get('sort')
        mode_revalidation_chgvs = options.get('revalidate_chgvs')
        mode_non_coding = options.get('non_coding')

        # if mode_all and mode_missing:
        #     raise ValueError("all and missing are mutually exclusive parameters")
        #
        # if mode_file and (mode_all or mode_missing or mode_extra, mode_validation):
        #     raise ValueError("file is an exclusive parameter")
        #
        # if not any([mode_all, mode_missing, mode_file, mode_extra, mode_validation]):
        #     raise ValueError("Need one of all, file, missing, mode_extra, mode_validation")

        if mode_extra:
            self.handle_extra()

        if mode_validation:
            self.handle_validation()

        if mode_sort:
            self.handle_fix_sort_order()

        if mode_revalidation_chgvs:
            self.handle_revalidate_chgvs()
        if mode_non_coding:
            self.handle_fix_non_coding()

        # row_count = 0
        # batch_size = 50
        # import_keys: Optional[Set[str]] = None
        # if mode_missing:
        #     batch_size = 1
        # elif mode_file:
        #     header_row = True
        #     import_keys = set()
        #     with open(mode_file, 'r') as file_handle:
        #         for row in csv.reader(file_handle):
        #             if header_row:
        #                 if not row[0] == "Import Key":
        #                     raise ValueError("Expected first column of file to be - 'Import Key'")
        #                 header_row = False
        #             else:
        #                 import_keys.add(row[0].strip())
        #     print(f"Only running over import keys in file - {len(import_keys)}")
        #
        # self.report_unmatched()
        #
        # qs = Classification.objects.all().order_by('evidence__genome_build__value')
        # if mode_missing:
        #     qs = qs.filter(variant__isnull=True)
        #
        # user = admin_bot()
        #
        # # setup a temporary import so discordance notifications are not sent out
        # try:
        #     batch: List[Classification] = []
        #     for c in qs:
        #         if import_keys is not None:
        #             import_key = f"{c.get(SpecialEKeys.GENOME_BUILD) or ''}#{c.get(SpecialEKeys.C_HGVS) or ''}"
        #             if import_key not in import_keys:
        #                 continue
        #
        #         c.revalidate(user=user)
        #         batch.append(c)
        #         row_count += 1
        #
        #         if len(batch) >= batch_size:
        #             self.handle_batch(batch)
        #             batch = []
        #             print(f"Handled {row_count}")
        #
        #     self.handle_batch(batch)
        #     print(f"Handled {row_count}")
        #     self.report_unmatched()
        #     sleep(10)  # give time for variant matching to complete
        # finally:
        #     ClassificationImportRun.record_classification_import("variant_rematching", 0, is_complete=True)

    def sleep_for_delay(self):
        time.sleep(10)

    def handle_fix_sort_order(self):
        print("Updating sort order - this may take a few minutes")
        updates = 0
        for rv in ResolvedVariantInfo.objects.all():
            actual_sort_string = rv.variant.sort_string
            if rv.genomic_sort != actual_sort_string:
                rv.genomic_sort = actual_sort_string
                rv.save()
                updates += 1
        print(f"Updated {updates} variant info")

    def handle_extra(self):
        i = 0
        for i, c in enumerate(Classification.objects.all()):
            if i % 100 == 0:
                print(f"Processed {i} classifications")
            # use force_update so we can be sure that the validation objects have been made
            if not c.allele_info:
                c.ensure_allele_info()
                c.save(update_modified=False)

            if c.update_allele_info_from_classification(force_update=False):
                c.save(update_modified=False)
        print(f"Finished {i} classifications")

    def handle_revalidate_chgvs(self):
        def c_hgvs_validation(c: Classification) -> str:
            if c_hgvs := c.evidence.get('c_hgvs'):
                if validation := c_hgvs.get('validation'):
                    return ", ".join(vm.get('message') for vm in validation)
            return "<no issues - but check imported allele info>"

        user = admin_bot()
        for c in Classification.objects.filter(evidence__c_hgvs__validation__isnull=False):
            before = c_hgvs_validation(c)
            # revalidate calls save
            c.revalidate(user=user)
            after = c_hgvs_validation(c)
            print(f"{c.id} from {before} to {after}")

    def handle_validation(self):
        def get_flag_comments(flag_type: FlagType, resolution_id: str) -> Dict[int, FlagComment]:
            flag_dict: Dict[int, FlagComment] = {}
            manual_closed_flag_comments = FlagComment.objects.filter(
                flag__flag_type=flag_type,
                resolution__id=resolution_id).exclude(user=admin_bot()) \
                .select_related('flag', 'flag__collection')
            for flag_comment in manual_closed_flag_comments:
                flag_collection = flag_comment.flag.collection
                # have to make sure we don't have an open flag of the type, only closed
                if not flag_collection.get_open_flag_of_type(
                        flag_type=flag_type):
                    flag_dict[flag_collection.pk] = flag_comment
            return flag_dict

        matching_variant_warning_flag_type = FlagType.objects.get(id='classification_matching_variant_warning')
        classification_transcript_version_change_flag_type = FlagType.objects.get(id='classification_transcript_version_change')
        flag_type_37_not_38 = FlagType.objects.get(pk='allele_37_not_38')

        manually_closed_variant_warning = get_flag_comments(flag_type=matching_variant_warning_flag_type, resolution_id='vm_confirmed')
        manually_closed_37_not_38 = get_flag_comments(flag_type=flag_type_37_not_38, resolution_id='closed')
        manually_closed_transcript_ver_changed = get_flag_comments(flag_type=classification_transcript_version_change_flag_type, resolution_id='tv_accepted')

        print(f"variant matching count = {len(manually_closed_variant_warning)}")
        print(f"37 not 38 flag count = {len(manually_closed_37_not_38)}")

        for i, allele_info in enumerate(ImportedAlleleInfo.objects.all()):
            if i % 100 == 0:
                print(f"Processed {i} allele infos")
            allele_info.apply_validation(force_update=True)
            allele_info.save()

            if allele := allele_info.allele:
                latest_validation = allele_info.latest_validation
                if not latest_validation.include:

                    has_normal_issues = False
                    has_liftover_issues = False
                    has_general_issues = False
                    has_builds_issues = False

                    if normalize_issues := latest_validation.validation_tags_typed.get("normalize"):
                        has_normal_issues = [True for severity in normalize_issues.values() if severity == "E"]
                    if liftover_issues := latest_validation.validation_tags_typed.get("liftover"):
                        has_liftover_issues = [True for severity in liftover_issues.values() if severity == "E"]
                    if general_issues := latest_validation.validation_tags_typed.get("general"):
                        has_general_issues = [True for severity in general_issues.values() if severity == "E"]
                    if build_issues := latest_validation.validation_tags_typed.get("builds"):
                        has_builds_issues = [True for severity in build_issues.values() if severity == "E"]

                    if not has_normal_issues and not has_liftover_issues and not has_general_issues and not has_builds_issues:
                        print("Why was this excluded in the first place if no E?")
                        print("Dev should investigate")
                        print(latest_validation.validation_tags_typed)
                        return

                    comments = set()
                    users = set()
                    if allele_flag_pk := allele.flag_collection_id:
                        if approved_flag := manually_closed_37_not_38.get(allele_flag_pk):
                            has_liftover_issues = False
                            if text := approved_flag.text:
                                comments.add(text)
                            users.add(approved_flag.user)

                    classifications = Classification.objects.filter(allele_info=allele_info)
                    for classification in classifications:
                        if approved_flag := manually_closed_variant_warning.get(classification.flag_collection_id):
                            has_normal_issues = False
                            if text := approved_flag.text:
                                comments.add(text)
                            users.add(approved_flag.user)

                    if not has_normal_issues and not has_liftover_issues and not has_general_issues and not has_builds_issues:
                        latest_validation.include = True
                        latest_validation.confirmed = True
                        latest_validation.confirmed_by = first(users)
                        latest_validation.confirmed_by_note = "\n".join(comments) if comments else ""
                        latest_validation.save()
                        print(f"Confirmed something to true on allele {allele}")

                elif tags := latest_validation.validation_tags_list:
                    transcript_ver_changes = [tag for tag in tags if tag.field == 'transcript_version_change']
                    if transcript_ver_changes and len(transcript_ver_changes) == len(tags):
                        classifications = Classification.objects.filter(allele_info=allele_info)
                        was_confirmed = False
                        comments = set()
                        users = set()
                        for classification in classifications:
                            if approved_flag := manually_closed_transcript_ver_changed.get(classification.flag_collection_id):
                                was_confirmed = True
                                if text := approved_flag.text:
                                    comments.add(text)
                                users.add(approved_flag.user)

                        if was_confirmed:
                            print("Confirmed transcript version change")
                            latest_validation.confirmed = True
                            latest_validation.confirmed_by = first(users)
                            latest_validation.confirmed_by_note = "\n".join(comments) if comments else ""
                            latest_validation.save()

    @staticmethod
    def handle_fix_non_coding():
        rvi_qs = ResolvedVariantInfo.objects.filter(allele_info__imported_c_hgvs__icontains='n.',
                                                    c_hgvs__icontains='c.')
        iai_ids_to_fix = list(rvi_qs.values_list("allele_info_id", flat=True).distinct())
        print(f"Fixing ResolvedVariantInfo for {len(iai_ids_to_fix)} ImportedAlleleInfo")
        rvi_qs.delete()

        for iai in ImportedAlleleInfo.objects.filter(pk__in=iai_ids_to_fix):
            try:
                variant_37 = iai.allele.variant_for_build(GenomeBuild.grch37())
                iai.grch37 = ResolvedVariantInfo.get_or_create(allele_info=iai, genome_build=GenomeBuild.grch37(),
                                                               variant=variant_37)
            except ValueError:
                pass

            try:
                variant_38 = iai.allele.variant_for_build(GenomeBuild.grch38())
                iai.grch38 = ResolvedVariantInfo.get_or_create(allele_info=iai, genome_build=GenomeBuild.grch38(),
                                                               variant=variant_38)
            except ValueError:
                pass

            iai.apply_validation(force_update=True)
            iai.save()


    # def handle_batch(self, batch: List[Classification]):
    #     ClassificationImportRun.record_classification_import("variant_rematching", len(batch))
    #     user = admin_bot()
    #     if batch:
    #         imports_by_genome: Dict[int, ClassificationImport] = {}
    #         for vc in batch:
    #             try:
    #                 genome_build = vc.get_genome_build()
    #                 if genome_build.pk not in imports_by_genome:
    #                     imports_by_genome[genome_build.pk] = ClassificationImport.objects.create(user=user,
    #                                                                                              genome_build=genome_build)
    #                 vc_import = imports_by_genome[genome_build.pk]
    #                 vc.set_variant(variant=None, message='Admin has re-triggered variant matching')
    #                 vc.classification_import = vc_import
    #                 vc.save()
    #             except ValueError as ve:
    #                 print(f"Couldn't revalidate {vc.id} due to bad genome build {ve}")
    #
    #         for vc_import in imports_by_genome.values():
    #             process_classification_import(vc_import, ImportSource.API)
    #
    #         self.sleep_for_delay()
