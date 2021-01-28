from datetime import datetime
from typing import List

import bs4
from django.contrib.auth.models import User
from django.db import models
from django.db.models import PROTECT, CASCADE
from django.db.models.signals import post_save
from django.dispatch import receiver
from guardian.shortcuts import assign_perm
from model_utils.models import TimeStampedModel

from classification.enums import ShareLevel, SpecialEKeys
from classification.models import ClassificationModification, EvidenceKeyMap, classification_post_publish_signal, \
    Classification, flag_types
from classification.models.evidence_mixin import VCDbRefDict
from classification.regexes import DbRegexes
from flags.models import flag_comment_action, Flag, FlagComment, FlagResolution, FlagStatus
from genes.hgvs import CHGVS
from genes.models import GeneSymbol
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from snpdb.models import Allele, Lab, GenomeBuild


class ClinVarExportStatus(models.TextChoices):
    SUBMIT_WHEN_READY = 'S', 'Submit When Ready'
    PENDING = 'P', 'Pending'
    REJECT = 'R', 'Reject'


class MultiCondition(models.TextChoices):
    NOT_DECIDED = 'N', 'Not decided'
    UNCERTAIN = 'U', 'Uncertain'  # aka uncertain
    CO_OCCURRING = 'C', 'Co-occurring'  # aka combined


class ClinVarExport(TimeStampedModel, GuardianPermissionsMixin):
    allele = models.ForeignKey(Allele, on_delete=PROTECT)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    classification_based_on = models.ForeignKey(ClassificationModification, on_delete=PROTECT)
    transcript = models.TextField()
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=PROTECT)

    review_date = models.DateTimeField(null=True, blank=True)
    review_status = models.CharField(max_length=1, choices=ClinVarExportStatus.choices, default=ClinVarExportStatus.PENDING)
    dirty_date = models.DateTimeField(null=True, blank=True)

    submit_when_possible = models.BooleanField(default=False)
    withdrawn = models.BooleanField(default=False)

    @staticmethod
    def best_clinvar_candidate(cm1: ClassificationModification, cm2: ClassificationModification) -> ClassificationModification:
        # FIXME created date isn't the best, last curated would be more accurate
        if cm1.classification.created > cm2.classification.created:
            return cm1
        return cm2

    def update_with(self, cm: ClassificationModification) -> bool:
        is_update = self.id is None or self.classification_based_on != cm
        self.classification_based_on = cm
        # self.dirty_date = datetime.now()
        # only update dirty date if json that we will send changes
        return is_update

    @staticmethod
    def chgvs_for(cm: ClassificationModification):
        # FIXME, make 38 only
        if chgvs_str := cm.classification.chgvs_grch37:
            return CHGVS(chgvs_str)
        elif chgvs_str := cm.get(SpecialEKeys.C_HGVS):
            return CHGVS(chgvs_str)
        return None

    @staticmethod
    def sync_allele(allele: Allele):

        new_count = 0
        updated_count = 0
        orphan_count = 0

        transcript_lab_to_cm = dict()

        cm: ClassificationModification
        # TODO - do we reduce this to just records shared globally??!
        # or do we accept both, and treat the share level as a to be confirmed globally
        for cm in ClassificationModification.objects.filter(
                classification__withdrawn=False,
                classification__variant__variantallele__allele=allele,
                is_last_published=True,
                share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).select_related('classification'):

            classification = cm.classification
            # decision, grabbing the 38 representation
            if c_parts := ClinVarExport.chgvs_for(cm):
                transcript_no_version = c_parts.transcript_parts.identifier
                lab_str = str(classification.lab_id)
                db_key = f"{lab_str}*{transcript_no_version}"
                use_cm = cm
                if existing := transcript_lab_to_cm.get(db_key):
                    use_cm = ClinVarExport.best_clinvar_candidate(cm, existing)
                transcript_lab_to_cm[db_key] = use_cm

        cve: ClinVarExport
        existing_records = {f"{cve.lab.id}*{cve.transcript}": cve for cve in ClinVarExport.objects.filter(allele=allele).select_related('lab')}
        for lab_transcript, cm in transcript_lab_to_cm.items():
            cve: ClinVarExport
            if cve := existing_records.pop(lab_transcript, None):
                if cve.update_with(cm) or cve.withdrawn:
                    cve.withdrawn = False
                    cve.save()
                    updated_count += 1
            else:
                c_parts = ClinVarExport.chgvs_for(cm)
                # TODO put some saftey around this?
                gs: GeneSymbol = GeneSymbol.objects.get(pk=c_parts.gene)
                cve = ClinVarExport(
                    allele=allele,
                    lab=cm.classification.lab,
                    transcript=lab_transcript.split("*", maxsplit=1)[1],
                    gene_symbol=gs
                )
                cve.update_with(cm)
                cve.save()
                new_count += 1

        for orphan in existing_records.values():
            if not orphan.withdrawn:
                # TODO allow user to withdraw
                orphan.withdrawn = True
                orphan.save()
                orphan_count += 1

        return {
            "new": new_count,
            "updated": updated_count,
            "withdrawn": orphan_count
        }

    # convinience method for forms
    def evidence_key(self, key: str):
        return EvidenceKeyMap.cached().get(key).pretty_value(self.classification_based_on.get(key), dash_for_none=True)

    @property
    def allele_origin(self):
        return self.evidence_key(SpecialEKeys.ALLELE_ORIGIN)

    @property
    def mode_of_inheritance(self):
        return self.evidence_key(SpecialEKeys.MODE_OF_INHERITANCE)

    @property
    def affected_status(self):
        return self.evidence_key(SpecialEKeys.AFFECTED_STATUS)

    @property
    def genome_build(self):
        return self.evidence_key(SpecialEKeys.GENOME_BUILD)

    @property
    def interpretation_summary(self):
        # strip out any XML
        interpret = self.evidence_key(SpecialEKeys.INTERPRETATION_SUMMARY)
        soup = bs4.BeautifulSoup(interpret, 'lxml')
        return soup.text

    @property
    def curation_context(self):
        return self.evidence_key(SpecialEKeys.CURATION_CONTEXT)

    @property
    def assertion_method(self):
        return self.evidence_key(SpecialEKeys.ASSERTION_METHOD)

    # probably wont use this
    @property
    def patient_phenotype(self):
        return self.evidence_key(SpecialEKeys.PATIENT_PHENOTYPE)

    @property
    def clinical_significance(self):
        return self.evidence_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @property
    def citation_refs(self) -> List[VCDbRefDict]:
        pubmed_refs = [ref for ref in self.classification_based_on.db_refs if ref.get('db') == DbRegexes.PUBMED.db]
        unique_refs = set()

    @property
    def condition(self) -> str:
        """
        Note condition isn't actually sent, but the value from ConditionTextMatching
        this is just here for reference
        """
        return self.evidence_key(SpecialEKeys.CONDITION)

    @property
    def curated_date(self):
        if date_str := self.evidence_key(SpecialEKeys.CURATION_DATE):
            try:
                if date := datetime.datetime.strptime(date_str, '%Y-%m-%d'):
                    return date
            except:
                pass
        return self.classification_based_on.created

    @property
    def c_hgvs(self):
        gb = GenomeBuild.get_from_fuzzy_string(self.genome_build)
        return self.classification_based_on.classification.get_c_hgvs(genome_build=gb)

    @property
    def as_json(self):
        data = dict()
        data["genome_build"] = self.genome_build
        data["c_hgvs"] = self.c_hgvs
        data["Note"] = "This ClinVar API has not been published yet, so this is just placeholder"
        return data


@receiver(post_save, sender=ClinVarExport)
def set_condition_alias_permissions(sender, created: bool, instance: ClinVarExport, **kwargs):
    if created:
        group = instance.lab.group
        assign_perm(ClinVarExport.get_read_perm(), group, instance)
        assign_perm(ClinVarExport.get_write_perm(), group, instance)


class ClinVarExportSubmission(TimeStampedModel, GuardianPermissionsMixin):
    clinvar_export = models.ForeignKey(ClinVarExport, on_delete=CASCADE)
    evidence = models.JSONField()
    submission_status = models.TextField()


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              **kwargs):

    cve: ClinVarExport
    if cve := ClinVarExport.objects.filter(classification_based_on__classification=classification).first():
        cve.update_with(newly_published)
        cve.save()


@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):
    """
    Keeps condition_text_match in sync with the classifications when withdraws happen/finish
    """
    flag = flag_comment.flag
    if flag.flag_type == flag_types.classification_flag_types.classification_withdrawn:
        cl: Classification
        if cl := Classification.objects.filter(flag_collection=flag.collection.id).first():
            cve: ClinVarExport
            if cve := ClinVarExport.objects.filter(classification_based_on__classification=cl).first():
                cve.withdrawn = flag_comment.resolution.status == FlagStatus.OPEN
                cve.save()
