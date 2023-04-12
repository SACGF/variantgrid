import itertools
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Dict, List, Optional, Set

from classification.models import Classification, ImportedAlleleInfo
from flags.models import FlagComment, FlagType, FlagStatus
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from snpdb.models import Allele, GenomeBuild

FLAG_KEY_CLASSIFICATION_MATCHING_VARIANT_WARNING = 'classification_matching_variant_warning'
FLAG_KEY_CLASSIFICATION_TRANSCRIPT_VERSION_CHANGE = 'classification_transcript_version_change'
FLAG_KEY_37_NOT_38 = 'allele_37_not_38'


@dataclass(frozen=True)
class CHGVSIdentifier:

    genome_build: Optional[GenomeBuild] = None
    imported_c_hgvs: Optional[str] = None
    transcript: Optional[str] = None

    @property
    def effective_transcript(self) -> str:
        if self.transcript:
            return self.transcript
        return CHGVS(self.imported_c_hgvs).transcript

    def fuzzy_matches(self, other: 'CHGVSIdentifier') -> bool:
        if other.genome_build and self.genome_build and other.genome_build != self.genome_build:
            return False
        if other.imported_c_hgvs and self.imported_c_hgvs and other.imported_c_hgvs != self.imported_c_hgvs:
            return False
        if other.effective_transcript and self.effective_transcript and other.effective_transcript != self.effective_transcript:
            return False
        return True


@dataclass(frozen=True)
class AlleleIdentifier:

    allele_id: int
    identifier: Optional[CHGVSIdentifier] = None

    def with_transcript(self, transcript: str) -> 'AlleleIdentifier':
        return AlleleIdentifier(allele_id=self.allele_id, identifier=CHGVSIdentifier(transcript=transcript))


class CHGVSFlagData:

    def __init__(self, identifier: CHGVSIdentifier):
        self.identifier = identifier

        self.open_flags: Set[str] = set()
        self.manually_closed_flags: Set[str] = set()
        self.comments: List[FlagComment] = []

    def add_comment(self, flag_comment: FlagComment):
        flag_type = flag_comment.flag.flag_type.pk
        if not flag_comment.user == admin_bot():
            self.comments.append(flag_comment)
            if flag_comment.flag.resolution.status == FlagStatus.CLOSED:
                self.manually_closed_flags.add(flag_type)
        if flag_comment.flag.resolution.status == FlagStatus.OPEN:
            self.open_flags.add(flag_type)

    @property
    def closed_flags(self) -> Set[str]:
        return self.manually_closed_flags - self.open_flags


class AlleleData:

    def __init__(self):
        self.flag_datas = list()

    def flag_data_for_identifier(self, identifier: CHGVSIdentifier) -> CHGVSFlagData:
        for flag_data in self.flag_datas:
            if flag_data.identifier == identifier:
                return flag_data
        new_flag_data = CHGVSFlagData(identifier=identifier)
        self.flag_datas.append(new_flag_data)
        return new_flag_data

    def all_flag_data_for_c_hgvs(self, identifier: CHGVSIdentifier) -> List[CHGVSFlagData]:
        return [fd for fd in self.flag_datas if fd.identifier.fuzzy_matches(identifier)]


class FlagDatabase:

    def __init__(self):
        self.allele_to_flags: Dict[int, AlleleData] = defaultdict(AlleleData)

    @cached_property
    def flag_collection_to_identifier(self) -> Dict[int, AlleleIdentifier]:
        # flag collection to identifier
        flag_collection_to_allele_identifier = dict()

        # map classification flags to identifiers
        for classification_data in Classification.objects.values_list('flag_collection__id', 'allele_id',
                                                                                 'evidence__c_hgvs__value', 'evidence__genome_build__value').iterator():
            flag_collection_id, allele_id, c_hgvs, genome_build_str = classification_data
            genome_build: Optional[GenomeBuild] = None
            if genome_build_str:
                try:
                    genome_build = GenomeBuild.get_name_or_alias(genome_build_str)
                except GenomeBuild.DoesNotExist:
                    pass

            if genome_build:
                flag_collection_to_allele_identifier[flag_collection_id] = AlleleIdentifier(allele_id=allele_id, identifier=CHGVSIdentifier(genome_build=genome_build, imported_c_hgvs=c_hgvs))

        # map allele flags to identifiers
        for allele_data in Allele.objects.filter(flag_collection__id__isnull=False).values_list('flag_collection__id', 'pk').iterator():
            flag_collection_id, allele_id = allele_data
            flag_collection_to_allele_identifier[flag_collection_id] = AlleleIdentifier(allele_id=allele_id)

        return flag_collection_to_allele_identifier

    def flag_data_for_identifier(self, allele_identifier: AlleleIdentifier) -> CHGVSFlagData:
        return self.allele_to_flags[allele_identifier.allele_id].flag_data_for_identifier(identifier=allele_identifier.identifier)

    def all_flag_data_for_c_hgvs(self, allele_identifier: AlleleIdentifier) -> List[CHGVSFlagData]:
        return self.allele_to_flags[allele_identifier.allele_id].all_flag_data_for_c_hgvs(allele_identifier.identifier)

    def populate(self):
        matching_variant_warning_flag_type = FlagType.objects.get(pk=FLAG_KEY_CLASSIFICATION_MATCHING_VARIANT_WARNING)
        classification_transcript_version_change_flag_type = FlagType.objects.get(pk=FLAG_KEY_CLASSIFICATION_TRANSCRIPT_VERSION_CHANGE)
        flag_type_37_not_38 = FlagType.objects.get(pk=FLAG_KEY_37_NOT_38)

        # classification flags - comments
        for comment in FlagComment.objects.select_related('flag').filter(
            flag__flag_type__in={matching_variant_warning_flag_type, classification_transcript_version_change_flag_type}
        ).order_by('-created'):
            if allele_identifier := self.flag_collection_to_identifier.get(comment.flag.collection_id):
                self.flag_data_for_identifier(allele_identifier).add_comment(comment)

        # allele flags
        for comment in FlagComment.objects.select_related('flag').filter(
            flag__flag_type__in={flag_type_37_not_38}
        ).order_by('-created'):
            if flag_data := comment.flag.data:
                if transcript := flag_data.get('transcript'):
                    if allele_identifier := self.flag_collection_to_identifier.get(comment.flag.collection_id):
                        self.flag_data_for_identifier(allele_identifier.with_transcript(transcript)).add_comment(comment)

    def apply_to_imported_allele_infos(self):
        if not self.allele_to_flags:
            print("No relevant flags detected, nothing to do")

        for imported_allele_info in ImportedAlleleInfo.objects.select_related('latest_validation', 'imported_genome_build_patch_version__genome_build').iterator():
            latest_validation = imported_allele_info.latest_validation
            if c_hgvs := imported_allele_info.imported_c_hgvs:
                allele_identifier = AlleleIdentifier(imported_allele_info.allele_id, identifier=CHGVSIdentifier(imported_c_hgvs=c_hgvs, genome_build=imported_allele_info.imported_genome_build))

                if matches := self.all_flag_data_for_c_hgvs(allele_identifier):
                    comment_arrays = [match.comments for match in matches]
                    all_comments = list(itertools.chain(*comment_arrays))
                    all_closed_flags = set()
                    all_closed_flags.update(*(match.closed_flags for match in matches))

                    all_comments = list(sorted(all_comments, key=lambda c: c.created))

                    def convert_flag_comment(comment: FlagComment) -> str:
                        return f"({comment.created:%Y-%m-%d %H:%M} - {comment.flag.flag_type.label} @ {comment.user}) - {comment.text}"

                    comment_line = "\n\n".join(convert_flag_comment(comment) for comment in all_comments)

                    latest_validation.confirmed_by_note = comment_line
                    outstanding_issues = set()
                    for tag in latest_validation.validation_tags_list:
                        if tag.category == "normalize" and (tag.field == "c_nomen_change" or tag.field == "gene_symbol_change"):
                            outstanding_issues.add(FLAG_KEY_CLASSIFICATION_MATCHING_VARIANT_WARNING)
                        elif tag.field == "liftover" and (tag.field == "c_nomen_change" or tag.field == "gene_symbol_change"):
                            outstanding_issues.add(FLAG_KEY_37_NOT_38)
                        elif tag.severity == "E":
                            outstanding_issues.add("OTHER_ERROR_CANT_CONFIRM")

                    if latest_validation.validation_tags_list and outstanding_issues.issubset(all_closed_flags):
                        latest_validation.include = True
                        latest_validation.confirmed = True
                        latest_validation.confirmed_by = all_comments[0].user
                        latest_validation.save()
                    else:
                        # undo previous times where it was set
                        was_confirmed = latest_validation.confirmed

                        latest_validation.confirmed = False
                        latest_validation.confirmed_by = None
                        latest_validation.save()
                        if was_confirmed:
                            print(f"Unconfirmed AlleleInfo {imported_allele_info.pk}")
                            imported_allele_info.apply_validation(force_update=True)
                            imported_allele_info.save()

    @staticmethod
    def run():
        flag_database = FlagDatabase()
        print("Populating flag data into memory")
        flag_database.populate()
        print("Done")
        print("Applying flag data to imported allele info validation")
        flag_database.apply_to_imported_allele_infos()
        print("Done")