import json
import operator
import re
from dataclasses import dataclass
from functools import cached_property, reduce
from operator import attrgetter
from typing import Optional, Iterable

import django
from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.db import models, transaction
from django.db.models import CASCADE, QuerySet, SET_NULL, Q
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.urls import reverse
from guardian.shortcuts import assign_perm
from model_utils.models import TimeStampedModel

from annotation.regexes import db_ref_regexes
from classification.enums import SpecialEKeys, ShareLevel
from classification.models import Classification, ClassificationModification, classification_post_publish_signal, \
    flag_types, EvidenceKeyMap, ConditionResolvedDict, ConditionResolved, classification_flag_types
from classification.models.condition_text_search import condition_text_search
from flags.models import flag_comment_action, Flag, FlagComment, FlagResolution
from genes.models import GeneSymbol, GeneSymbolAlias
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import admin_bot
from library.log_utils import report_exc_info
from library.utils import ArrayLength, get_timer
from ontology.models import OntologyTerm, OntologyService, OntologySnake, OntologyTermRelation, OntologyRelation
from ontology.ontology_matching import normalize_condition_text, \
    OPRPHAN_OMIM_TERMS, SearchText, pretty_set, PREFIX_SKIP_TERMS, IGNORE_TERMS, NON_PR_TERMS
from snpdb.models import Lab


condition_set_signal = django.dispatch.Signal()  # args: "classification", "resolved_condition"


class ConditionText(TimeStampedModel, GuardianPermissionsMixin):
    """
    For each normalized text/lab combo there'll be one ConditionText.
    Then there will be a hierarchy of ConditionTextMatches for the ConditionText which can be assigned standard
    terms and then apply that to the corresponding classification.
    """
    normalized_text = models.TextField()
    lab = models.ForeignKey(Lab, on_delete=CASCADE, null=True, blank=True)

    classifications_count = models.IntegerField(default=0)
    classifications_count_outstanding = models.IntegerField(default=0)

    class Meta:
        unique_together = ("normalized_text", "lab")

    @staticmethod
    def normalize(text: str):
        return normalize_condition_text(text)

    def __str__(self):
        return self.normalized_text

    # TODO is this better done as an before save hook?
    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        assign_perm(self.get_read_perm(), self.lab.group, self)
        assign_perm(self.get_write_perm(), self.lab.group, self)

    @transaction.atomic
    def clear(self):
        """
        Clears any condition matching at any level for this condition text
        """
        self.conditiontextmatch_set.update(condition_xrefs=[], condition_multi_operation=MultiCondition.NOT_DECIDED, last_edited_by=None)
        self.classifications_count_outstanding = self.classifications_count
        self.save()

    @property
    def classification_match_count(self) -> int:
        """
        How many classifications with this text have valid terms assigned to them
        """
        return self.classifications_count - self.classifications_count_outstanding

    @property
    def user_edited(self) -> bool:
        """
        Has a user (other than admin bot) edited this
        """
        return self.conditiontextmatch_set.annotate(condition_xrefs_length=ArrayLength('condition_xrefs')).filter(condition_xrefs_length__gt=0).exclude(last_edited_by=admin_bot()).exists()

    @property
    def root(self) -> 'ConditionTextMatch':
        return self.conditiontextmatch_set.filter(gene_symbol__isnull=True, mode_of_inheritance__isnull=True, classification=None).first()

    @property
    def gene_levels(self) -> QuerySet['ConditionTextMatch']:
        return self.conditiontextmatch_set.filter(condition_text=self, gene_symbol__isnull=False, mode_of_inheritance__isnull=True, classification=None)

    def get_absolute_url(self):
        return reverse("condition_matching", kwargs={"pk": self.pk})


class MultiCondition(models.TextChoices):
    """
    For when multiple terms exist, are we uncertain which one it is, or are we curating for multiple.
    NOT_DECIDED means a choice is still required
    """
    NOT_DECIDED = 'N', 'Not decided'
    UNCERTAIN = 'U', 'Uncertain'  # aka uncertain
    CO_OCCURRING = 'C', 'Co-occurring'  # aka combined

    @property
    def clinvar_label(self):
        return self.label


class ConditionTextMatch(TimeStampedModel, GuardianPermissionsMixin):
    """
    Represents a level of hierarchy for given ConditionText
    - root
    -- gene (e.g. BRCA1)
    --- mode_of_inheritance (e.g. [autosomal_dominant] (multiple modes are also supported)
    ---- classification X
    Where each item on the hierarchy should have 1 to many direct children
    To resolve classification X, go up the chain to find the first valid.

    To determine position in hierarchy see which elements are null.
    e.g. gene_symbol null means root level, gene_symbol with value but mode of inheritance null means gene level, etc
    """
    condition_text = models.ForeignKey(ConditionText, on_delete=CASCADE)
    last_edited_by = models.ForeignKey(User, on_delete=SET_NULL, null=True, blank=True)

    parent = models.ForeignKey('ConditionTextMatch', on_delete=CASCADE, null=True, blank=True)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE, null=True, blank=True)
    mode_of_inheritance = ArrayField(models.TextField(blank=False), default=None, null=True, blank=True)
    classification = models.OneToOneField(Classification, on_delete=CASCADE, null=True, blank=True)

    @property
    def name(self):
        """
        name of the item as it pertains to the hierarchy
        """
        if self.classification:
            return self.classification.friendly_label
        if self.mode_of_inheritance:
            e_key = EvidenceKeyMap.cached_key(SpecialEKeys.MODE_OF_INHERITANCE)
            return e_key.pretty_value(self.mode_of_inheritance)
        if self.gene_symbol:
            return self.gene_symbol.symbol
        return "Default"

    def get_permission_object(self):
        return self.condition_text

    @classmethod
    def get_permission_class(cls):
        return ConditionText

    @property
    def is_root(self):
        return not self.gene_symbol

    @property
    def is_gene_level(self):
        # just used for templates
        return not self.classification_id and self.mode_of_inheritance is None and self.gene_symbol_id is not None

    @cached_property
    def children(self) -> QuerySet['ConditionTextMatch']:
        order_by: str
        if self.is_root:
            order_by = 'gene_symbol__symbol'
        elif self.is_gene_level:
            order_by = 'mode_of_inheritance'
        else:
            order_by = 'classification__id'

        return self.conditiontextmatch_set.all().order_by(order_by)

    # resolved to this,
    condition_xrefs = ArrayField(models.TextField(blank=False), default=list)
    condition_multi_operation = models.CharField(max_length=1, choices=MultiCondition.choices, default=MultiCondition.NOT_DECIDED)

    @property
    def condition_matching_str(self) -> str:
        """
        Produce a single efficient line on what this condition text match is, e.g.
        MONDO:00001233; MONDO:0553533; co-occurring
        """
        result = ''
        if self.condition_xrefs:
            result = ", ".join(self.condition_xrefs)
        if len(self.condition_xrefs) >= 2:
            if self.condition_multi_operation == MultiCondition.NOT_DECIDED:
                result += "; uncertain/co-occurring"
            elif self.condition_multi_operation == MultiCondition.CO_OCCURRING:
                result += "; co-occurring"
            elif self.condition_multi_operation == MultiCondition.UNCERTAIN:
                result += "; uncertain"
        return result

    class Meta:
        unique_together = ("condition_text", "gene_symbol", "mode_of_inheritance", "classification")

    @property
    def is_valid(self):
        """
        Could this be used as a valid resolution of condition terms, if blank or significant errors then false
        Used to determine what parts of the hierarchy are valid
        """
        if len(self.condition_xrefs) == 0:
            return False
        if len(self.condition_xrefs) == 1:
            return True
        # we have 2+ terms from this point
        if self.condition_multi_operation == MultiCondition.NOT_DECIDED:
            # requires a joiner
            return False
        ontology_services = set()
        for xref in self.condition_xrefs:
            parts = xref.split(":")
            if len(parts) == 2:
                ontology_services.add(parts[0])
        if len(ontology_services) >= 2:
            # can't mix ontology terms
            return False

        return True

    @cached_property
    def condition_xref_terms(self) -> list[OntologyTerm]:
        terms = []
        for term_str in self.condition_xrefs:
            try:
                terms.append(OntologyTerm.get_or_stub(term_str))
            except ValueError:
                pass
        return terms

    @property
    def is_blank(self):
        return len(self.condition_xrefs) == 0

    @staticmethod
    def sync_all():
        """
        syncs ConditionTextMatch's to Classifications
        """
        cms = ClassificationModification.objects.filter(is_last_published=True, classification__withdrawn=False).select_related("classification", "classification__lab")
        cm: ClassificationModification
        for cm in cms:
            ConditionTextMatch.sync_condition_text_classification(cm=cm, update_counts=False)

        for ct in ConditionText.objects.all():
            ConditionTextMatch.attempt_automatch(condition_text=ct)
            update_condition_text_match_counts(ct)
            if ct.classifications_count == 0:
                # classifications have moved around
                ct.delete()
            else:
                ct.save()

    @staticmethod
    def attempt_automatch(condition_text: ConditionText, gene_symbol: Optional[str] = None):
        """
        Set terms that we're effectively certain of, do not override what's already there
        """
        try:
            if root := condition_text.root:
                if match := top_level_suggestion(condition_text.normalized_text):
                    if match.is_auto_assignable():
                        if not root.condition_xrefs:
                            root.condition_xrefs = match.term_str_array
                            root.condition_multi_operation = match.condition_multi_operation
                            root.last_edited_by = admin_bot()
                            root.save()
                    else:
                        gene_levels_qs = condition_text.gene_levels
                        if gene_symbol:
                            gene_levels_qs = gene_levels_qs.filter(gene_symbol=gene_symbol)

                        for gene_symbol_level in gene_levels_qs:
                            if match.is_auto_assignable(gene_symbol=gene_symbol_level.gene_symbol) and not gene_symbol_level.condition_xrefs:
                                gene_symbol_level.condition_xrefs = match.term_str_array
                                gene_symbol_level.condition_multi_operation = match.condition_multi_operation
                                gene_symbol_level.last_edited_by = admin_bot()
                                gene_symbol_level.save()
            update_condition_text_match_counts(condition_text)
            condition_text.save()
        except Exception:
            report_exc_info()

    @staticmethod
    def sync_condition_text_classification(cm: ClassificationModification, update_counts=True, attempt_automatch=False):
        """
        @param cm ClassificationModification - verify or create ConditionTextMatches for this classification
        @param update_counts - Update the counts on the ConditionText (set to False if you intend update counts of all ConditionTexts after)
        @param attempt_automatch - If true, attempt to automatch on the resulting ConditionText
        This method will save any created ConditionTexts
        """
        debug_timer = get_timer()
        debug_timer.tick("Condition Text Matching - PRE")

        classification = cm.classification
        existing: ConditionTextMatch = ConditionTextMatch.objects.filter(classification=classification).first()

        if classification.withdrawn:
            # don't worry about withdrawn records
            if existing:
                ct = existing.condition_text
                existing.delete()
                if update_counts:
                    update_condition_text_match_counts(ct)
            debug_timer.tick("Condition Text Matching - withdrawn")
            return

        lab = classification.lab
        gene_symbol: Optional[GeneSymbol] = None
        gene_symbol_str = cm.get(SpecialEKeys.GENE_SYMBOL)
        if gene_symbol_str:
            gene_symbol = GeneSymbol.objects.filter(symbol=gene_symbol_str).first()

        # see if there's a single gene symbol alias that this links to, if it's not a gene symbol itself
        if not gene_symbol:
            if alias_qs := GeneSymbolAlias.objects.filter(alias__iexact=gene_symbol_str):
                if alias_qs.count() == 1:
                    gene_symbol = alias_qs.first().gene_symbol

        # if the imported gene symbol doesn't work, see what gene symbols we can extract from the c.hgvs
        if not gene_symbol:
            try:
                genome_build = cm.get_genome_build()
                resolved_gene_symbol_str = cm.c_hgvs_best(genome_build).gene_symbol
                gene_symbol = GeneSymbol.objects.filter(symbol=resolved_gene_symbol_str).first()
            except ValueError:
                # couldn't extract genome build
                pass

        if not gene_symbol:
            # if gene_symbol_str:
            #     report_message("Classification has unrecognised gene symbol, cannot link it to condition text",
            #                    extra_data={
            #                        "target": gene_symbol_str or "<blank>",
            #                        "classification_id": classification.id,
            #                        "gene_symbol": gene_symbol_str
            #                    })
            debug_timer.tick("Condition Text Matching - no valid gene symbol")
            return

        raw_condition_text = cm.get(SpecialEKeys.CONDITION) or ""
        normalized = ConditionText.normalize(raw_condition_text)

        # need mode_of_inheritance to be not null
        mode_of_inheritance = cm.get(SpecialEKeys.MODE_OF_INHERITANCE) or []

        ct: ConditionText
        ct, ct_is_new = ConditionText.objects.get_or_create(normalized_text=normalized, lab=lab)

        debug_timer.tick("Condition Text Matching - gene symbol resolution")

        with transaction.atomic():
            ct = ConditionText.objects.select_for_update().filter(pk=ct.pk).first()

            # ensure each step of the hierarchy is present
            root, new_root = ConditionTextMatch.objects.get_or_create(
                condition_text=ct,
                gene_symbol=None,
                mode_of_inheritance=None,
                classification=None
            )

            gene_level, new_gene_level = ConditionTextMatch.objects.get_or_create(
                parent=root,
                condition_text=ct,
                gene_symbol=gene_symbol,
                mode_of_inheritance=None,
                classification=None
            )

            mode_of_inheritance_level, _ = ConditionTextMatch.objects.get_or_create(
                parent=gene_level,
                condition_text=ct,
                gene_symbol=gene_symbol,
                mode_of_inheritance=mode_of_inheritance,
                classification=None
            )

            debug_timer.tick("Condition Text Matching - ensure hierarchy")

            if existing:
                if existing.parent != mode_of_inheritance_level or \
                        existing.condition_text != ct or \
                        existing.gene_symbol != gene_symbol or \
                        existing.mode_of_inheritance != mode_of_inheritance:

                    # update existing to new hierarchy
                    # assume if a condition has been set for this classification specifically that it's
                    # still valid
                    old_text = existing.condition_text
                    existing.parent = mode_of_inheritance_level
                    existing.condition_text = ct
                    existing.mode_of_inheritance = mode_of_inheritance
                    existing.save()
                    if update_counts:
                        update_condition_text_match_counts(old_text)
                        old_text.save()
                    debug_timer.tick("Condition Text Matching - change existing")
                else:
                    debug_timer.tick("Condition Text Matching - no change")
                    # nothing has changed, no need to update anything
                    return
            else:
                ConditionTextMatch.objects.create(
                    parent=mode_of_inheritance_level,
                    condition_text=ct,
                    gene_symbol=gene_symbol,
                    mode_of_inheritance=mode_of_inheritance,
                    classification=classification
                )
                debug_timer.tick("Condition Text Matching - create new entry")

            if attempt_automatch and (new_root or new_gene_level):
                ConditionTextMatch.attempt_automatch(ct, gene_symbol=gene_symbol)
                debug_timer.tick("Condition Text Matching - auto match")
            elif update_counts:
                ct.classifications_count += 1
                is_valid = root.is_valid or gene_level.is_valid or mode_of_inheritance_level.is_valid or (existing and existing.is_valid)
                if not is_valid:
                    ct.classifications_count_outstanding += 1
                ct.save()
                debug_timer.tick("Condition Text Matching - update count quick")

    def as_resolved_condition(self) -> Optional[ConditionResolvedDict]:
        """
        Converts this specific ConditionTextMatch to a ConditionResolvedDict (suitable for saving on a classification).
        will not look up the hierarchy for a valid term, do that beforehand
        """
        if terms := self.condition_xref_terms:

            condition_resolved_obj = ConditionResolved(terms=_sort_terms(terms), join=None if len(terms) <= 1 else MultiCondition(self.condition_multi_operation))
            return condition_resolved_obj.to_json()

        return None


def update_condition_text_match_counts(ct: ConditionText):
    """
    condition count texts are cached in DB, call this method to make sure they're up to date
    does not save after
    """
    by_id: dict[int, ConditionTextMatch] = {}
    classification_related: list[ConditionTextMatch] = []
    for ctm in ct.conditiontextmatch_set.all():
        by_id[ctm.id] = ctm
        if ctm.classification_id:
            classification_related.append(ctm)

    def check_hierarchy(ctm: ConditionTextMatch) -> bool:
        if ctm.is_valid:
            return True
        if not ctm.is_blank:
            return False
        if parent := by_id.get(ctm.parent_id):
            return check_hierarchy(parent)
        return False

    invalid_count = 0
    for ctm in classification_related:
        if not check_hierarchy(ctm):
            invalid_count += 1

    ct.classifications_count = len(classification_related)
    ct.classifications_count_outstanding = invalid_count


@dataclass
class ConditionMatchingMessage:
    severity: str
    text: str

    def as_json(self):
        return {
            "severity": self.severity,  # todo should have severity enum
            "text": self.text
        }


class ConditionMatchingSuggestion:

    def __init__(self, condition_text_match: Optional[ConditionTextMatch] = None, ignore_existing: bool = False):
        """
        @param condition_text_match The ConditionTextMatch this suggestion is for (Optional)
        @param ignore_existing If true, ignore any existing assignments
        (useful for producing debug files of how we think all condition matching should go)
        """
        self.condition_text_match = condition_text_match
        self.terms: list[OntologyTerm] = _sort_terms(condition_text_match.condition_xref_terms) if condition_text_match and not ignore_existing else []
        """ A list of terms that are being suggested """

        self.is_applied = bool(self.terms)
        """ if true, this has been saved and retrieved from the DB, otherwise we're making it as a suggestion """

        self.info_only = False
        """ is this suggestion just validation on how a parent term applies to this level of the hierarchy, e.g. has messages but no terms """

        self.condition_multi_operation: MultiCondition = MultiCondition.NOT_DECIDED
        """ the way to deal with multiple terms, only set by users at this stage """
        if self.is_applied:
            self.condition_multi_operation = condition_text_match.condition_multi_operation

        self.messages: list[ConditionMatchingMessage] = []
        """ errors, warnings, info to inform the user of - or use to determine if this auto-assign appropriate """

        self.validated = False
        """ has validate been run yet"""

        self.ids_found_in_text: Optional[bool] = None
        """ Are the terms we found just straight from the condition text e.g. 'OMIM:32433', useful to determine auto-assign """

        self.alias_index: Optional[int] = None
        """ If not null, the index of the alias of the term we matched via text, useful to determine auto-assign """

        self.merged = False
        """ Was this suggestion merged, e.g. if there was an condition text 'BAM' that matched Best At Motoneuron and Blood Attacked Myliver """

    @property
    def term_str_array(self) -> list[str]:
        """ For saving back to ConditionTextMatches """
        return [term.id for term in self.terms]

    def add_term(self, term: OntologyTerm):
        if term not in self.terms:
            self.terms.append(term)

    def add_message(self, message: ConditionMatchingMessage):
        self.messages.append(message)

    def as_json(self):
        user_json = None
        if self.terms and self.is_applied:  # only report on user who filled in values
            if condition_text_match := self.condition_text_match:
                if user := condition_text_match.last_edited_by:
                    user_json = {"username": user.username}

        return {
            "id": self.condition_text_match.id if self.condition_text_match else None,
            "is_applied": self.is_applied,
            "info_only": self.info_only,
            "terms": [{"id": term.id, "name": term.name, "definition": '???' if term.is_stub else term.definition} for term in self.terms],
            "joiner": self.condition_multi_operation,
            "messages": [message.as_json() for message in self.messages if message.severity != 'debug'],
            "user": user_json
        }

    def is_auto_assignable(self, gene_symbol: Optional[GeneSymbol] = None):
        """ Is this suggestion so certain we can just assign it
        Has to be a single term (since we don't support uncertain/co-occurring in text yet).
        Can't be found via an alias.
        Can't be any warnings or errors.
        If we're a gene symbol, the term has to be a leaf term and have a relationship to the gene.
        Or if we're top level, the ID has to be found within text.
        """
        if terms := self.terms:
            if len(terms) != 1 and self.condition_multi_operation == MultiCondition.NOT_DECIDED:
                return False

            if self.alias_index is not None:
                return False

            for message in self.messages:
                if message.severity not in {"success", "info", "debug"}:
                    return False

            if gene_symbol:
                if self.is_all_leafs():
                    # if we're at a gene level, and we have a relationship and we're leafs
                    term = terms[0]
                    if OntologySnake.has_gene_relationship(term, gene_symbol):
                        return True
            else:
                # embedded ID is the only thing that will give you a top level assignment
                if self.ids_found_in_text:
                    return True
        return False

    def validate(self):
        """
        Perform validation to populate messages
        """
        if self.validated:
            return
        self.validated = True
        if terms := self.terms:
            # validate if we have multiple terms without a valid joiner
            if len(terms) > 1 and self.condition_multi_operation not in {MultiCondition.UNCERTAIN,
                                                                         MultiCondition.CO_OCCURRING}:
                self.add_message(ConditionMatchingMessage(severity="error",
                                                          text="Multiple terms provided, requires co-occurring/uncertain"))

            ontology_services: set[str] = set()
            for term in terms:
                ontology_services.add(term.ontology_service)
            if len(ontology_services) > 1:
                self.add_message(ConditionMatchingMessage(severity="error",
                                                          text=f"Only one ontology type is supported per level, {' and '.join(ontology_services)} found"))

            if valid_terms := [term for term in terms if not term.is_stub]:
                # validate that the terms have a known gene association if we're at gene level
                if ctm := self.condition_text_match:
                    if ctm.is_gene_level:
                        gene_symbol = self.condition_text_match.gene_symbol
                        for term in valid_terms:
                            if not OntologySnake.has_gene_relationship(term, gene_symbol):
                                self.add_message(ConditionMatchingMessage(severity="warning",
                                                                          text=f"{term.id} : no direct association on file to {gene_symbol.symbol}"))
                            else:
                                self.add_message(ConditionMatchingMessage(severity="success",
                                                                          text=f"{term.id} : is associated to {gene_symbol.symbol}"))

            # validate that we have the terms being referenced (if we don't big chance that they're not valid)
            for term in terms:
                if term.is_stub:
                    if term.ontology_service in OntologyService.LOCAL_ONTOLOGY_PREFIXES:
                        text = f"{term.id} : no copy of this term in our system"

                        self.add_message(ConditionMatchingMessage(severity="warning", text=text))
                    else:
                        self.add_message(ConditionMatchingMessage(severity="info",
                                                                  text=f"We do not store {term.ontology_service}, please verify externally"))
                elif term.is_obsolete:
                    self.add_message(
                        ConditionMatchingMessage(severity="error", text=f"{term.id} : \"{term.warning_text}\""))

    def is_all_leafs(self):
        if terms := self.terms:
            for term in terms:
                if not term.is_leaf:
                    return False
            return True
        return None

    def __bool__(self):
        return bool(self.terms) or bool(self.messages)


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              **kwargs):
    """
    Keeps condition_text_match in sync with the classifications when evidence changes
    """
    get_timer().tick("Condition Text Matching - post publish")
    ConditionTextMatch.sync_condition_text_classification(newly_published, attempt_automatch=True, update_counts=True)


@receiver(flag_comment_action, sender=Flag)
def check_for_withdrawn(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):
    """
    Keeps condition_text_match in sync with the classifications when withdraws start/finish
    """
    flag = flag_comment.flag
    if flag.flag_type == flag_types.classification_flag_types.classification_withdrawn:
        cl: Classification
        if cl := Classification.objects.filter(flag_collection=flag.collection.id).first():
            ConditionTextMatch.sync_condition_text_classification(cl.last_published_version, attempt_automatch=True, update_counts=True)


# @timed_cache(size_limit=2)
def top_level_suggestion(text: str) -> ConditionMatchingSuggestion:
    """ Make a suggestion at the root level for the given normalised text """
    if suggestion := embedded_ids_check(text):
        return suggestion
    return search_suggestion(text)


def embedded_ids_check(text: str) -> ConditionMatchingSuggestion:
    """
    Check the condition text for IDs e.g. 'OMIM:343233' this makes things easy
    Will also validate if there's an embedded ID with text that doesn't seem to match the term
    """
    cms = ConditionMatchingSuggestion()

    # Check if text ends with "uncertain" or "cooccurring"
    condition_check_text = text.lower().replace('-', '').strip()
    if condition_check_text.endswith("uncertain"):
        cms.condition_multi_operation = MultiCondition.UNCERTAIN
    elif condition_check_text.endswith("cooccurring"):
        cms.condition_multi_operation = MultiCondition.CO_OCCURRING

    db_matches = db_ref_regexes.search(text)
    detected_any_ids = bool(db_matches)  # see if we found any prefix suffix, if we do,
    db_matches = [match for match in db_matches if match.db in OntologyService.CONDITION_ONTOLOGIES]

    for match in db_matches:
        cms.add_term(OntologyTerm.get_or_stub(match.id_fixed))

    found_stray_omim = False

    # fall back to looking for stray OMIM terms if we haven't found any ids e.g. PMID:123456 should stop this code
    if not detected_any_ids and ':' not in text:
        stray_omim_matches = OPRPHAN_OMIM_TERMS.findall(text)
        stray_omim_matches = [term for term in stray_omim_matches if len(term) == 6]
        if stray_omim_matches:
            for omim_index in stray_omim_matches:
                omim = OntologyTerm.get_or_stub(f"OMIM:{omim_index}")
                cms.add_term(omim)
                cms.add_message(ConditionMatchingMessage(severity="warning", text=f"OMIM:{omim_index} was found without a prefix"))

    if cms.terms:
        cms.ids_found_in_text = True
        all_text_tokens = []
        all_term_tokens = []
        all_extra_words = []
        all_ignore_words = []
        condition_matching_terms = {'uncertain', 'cooccurring', 'un-certain', 'co-occurring'}
        # Split the input string by spaces to get individual terms
        input_terms = text.split()

        # Remove terms that match any in the set
        filtered_terms = [term for term in input_terms if term.lower() not in condition_matching_terms]

        # Join the filtered terms back into a string
        text = ' '.join(filtered_terms)
        for matched_term in cms.terms:
            if not matched_term.is_stub:
                ontology_prefixes = {ontology.lower() for ontology in OntologyService.CONDITION_ONTOLOGIES}
                ontology_prefixes.add("mim")
                found = False
                result_text = text
                for prefix in ontology_prefixes:
                    # fetch text between the matched term and the next ontology ID
                    pattern = fr"{matched_term.id.lower()}(.*?)({prefix}:\d+)"
                    match = re.search(pattern, text)
                    if match:
                        found = True
                        text1 = match.group(1).strip()
                        if len(text1) > 0:
                            result_text = text1
                        else:
                            result_text = None
                        break
                if not found:
                    # Find all text following the Find_id if no matching ontology found at end
                    match2 = re.search(f"{matched_term.id.lower()}(.*)", text)
                    if match2:
                        text2 = match2.group(1).strip()
                        if len(text2) > 0:
                            result_text = text2

                # work out all the words in the text, subtract any ontology IDs from the text
                text_tokens = SearchText.tokenize_condition_text(normalize_condition_text(result_text),
                                                                 deplural=True, deroman=True)
                origin_text_tokens = set(text_tokens)
                # remove all numbers that are the same as the OMIM or MONDO number
                ontology_terms = set()
                number_part = matched_term.id.split(":")[1]
                ontology_terms.add(number_part)
                while number_part.startswith("0"):
                    number_part = number_part[1:]
                    ontology_terms.add(number_part)

                ontology_prefixes = {ontology.lower() for ontology in OntologyService.CONDITION_ONTOLOGIES}
                ontology_prefixes.add("mim")
                # remove any terms with a : or # that look like they might be an ontology ID
                ontology_terms |= ontology_prefixes
                for term in text_tokens:
                    if ":" in term or "#" in term:
                        for ontology_prefix in ontology_prefixes:
                            if ontology_prefix in term:
                                ontology_terms.add(term)
                                break
                text_tokens -= ontology_terms
                text_tokens.discard(matched_term.id.lower())

                # now find all the terms that should be expected in the term
                term_tokens: set[str] = set()
                term_tokens = SearchText.tokenize_condition_text(normalize_condition_text(matched_term.name),
                                                                 deplural=True, deroman=True)
                term_tokens.discard(matched_term.id.lower())

                def check_has_non_pr_term(tokens: set[str]):
                    for term in term_tokens:
                        for non_pr in NON_PR_TERMS:
                            if non_pr in term:
                                return True
                    return False

                non_pr_terms = check_has_non_pr_term(term_tokens)

                if aliases := matched_term.aliases:
                    for alias in aliases:
                        term_tokens = term_tokens.union(
                            SearchText.tokenize_condition_text(normalize_condition_text(alias),
                                                               deplural=True, deroman=True))
                if extra := matched_term.extra:
                    if included_titles := extra.get('included_titles'):
                        for included_title in included_titles.split(';;'):
                            term_tokens = term_tokens.union(
                                SearchText.tokenize_condition_text(normalize_condition_text(included_title),
                                                                   deplural=True, deroman=True))

                if len(text_tokens) > 0:
                    all_text_tokens.append(sorted(text_tokens))
                if len(term_tokens) > 0:
                    all_term_tokens.append(sorted(term_tokens))
                extra_words = text_tokens.difference(term_tokens)
                raw_extra_words = set(extra_words)
                if len(extra_words) > 0:
                    all_extra_words.append(sorted(extra_words))
                extra_words = extra_words - PREFIX_SKIP_TERMS - IGNORE_TERMS
                gene_symbols_in_text = {word for word in extra_words if
                                        GeneSymbol.objects.filter(symbol=word).exists()}
                ignored_words = raw_extra_words - extra_words
                if len(ignored_words) > 0:
                    all_ignore_words.append(sorted(ignored_words))

                passed = False

                def includes_gene_symbol_related() -> str:
                    for gene_symbol in gene_symbols_in_text:
                        gene_symbol = gene_symbol.lower()
                        if f'{gene_symbol}-related' in text.lower():
                            return gene_symbol

                if gene_symbol_related := includes_gene_symbol_related():
                    cms.add_message(ConditionMatchingMessage(severity="info", text=f"Found {matched_term.id} in text, marked as related to \"{gene_symbol_related}\""))
                    passed = True
                elif len(extra_words) >= 3:  # 3 extra words and for words that are in common aren't longer than 9 letters combined
                    cms.add_message(ConditionMatchingMessage(severity="warning", text=f"Found {matched_term.id} in text, but also apparently unrelated words : {pretty_set(extra_words)}"))
                else:
                    passed = True

                if not passed and non_pr_terms:
                    cms.add_message(ConditionMatchingMessage(severity="info", text="Different terminology likely due to avoiding outdated wording in official term."))
        count_all_extra_words = 0
        for sublist in all_extra_words:
            if sublist:  # checks if sublist is not empty
                count_all_extra_words += len(sublist)

        cms.add_message(
            ConditionMatchingMessage(severity="debug", text=f"TEXT words (RAW) = {json.dumps(all_text_tokens)}"))
        cms.add_message(
            ConditionMatchingMessage(severity="debug", text=f"TERM words (RAW) = {json.dumps(all_term_tokens)}"))
        cms.add_message(
            ConditionMatchingMessage(severity="debug", text=f"Diff words (RAW) = {json.dumps(all_extra_words)}"))
        cms.add_message(ConditionMatchingMessage(severity="debug",
                                                 text=f"Ignored common words/gene symbols = {json.dumps(all_ignore_words)}"))
        cms.add_message(ConditionMatchingMessage(severity="debug",
                                                 text=f"Diff words (Processed) = {json.dumps(all_extra_words)}, count = {count_all_extra_words}"))
    cms.validate()
    return cms


def search_text_to_suggestion(search_text: SearchText, term: OntologyTerm) -> ConditionMatchingSuggestion:
    """
    Does this term match the search_text, if so, add it as a suggestion (possibly with some validation)
    """
    cms = ConditionMatchingSuggestion()
    if match_info := search_text.matches(term):
        if match_info.alias_index is not None:  # 0 alias is complete with acronym, 1 alias without acronym
            safe_alias = False
            if term.ontology_service == OntologyService.OMIM:
                alias = term.aliases[match_info.alias_index]
                if alias in [name_part.strip() for name_part in term.name.split(";")]:
                    # alias is part one of the name parts, would barely refer to it as an alias
                    safe_alias = True

            if term.ontology_service == OntologyService.MONDO:
                if omim_term := OntologyTermRelation.as_omim(term):
                    omim_name = omim_term.name
                    parts = [p.strip() for p in omim_name.split(';')]
                    full_name = parts[0]
                    # we could loop through all the name parts, but don't want to as the OMIM acronym is a little too ambitious
                    # if search_text.effective_equals(SearchText(full_name)):
                    #     safe_alias = True  # still mark it as True so we don't have a validation message
                    #     cms.alias_index = match_info.alias_index
                    #     cms.matched_via_alias_and_exact_term = omim_term.id

            if not safe_alias:
                cms.add_message(ConditionMatchingMessage(severity="info", text=f"Text matched on alias of {term.id}"))
                cms.alias_index = match_info.alias_index
        if term.ontology_service == OntologyService.OMIM:
            if mondo := OntologyTermRelation.as_mondo(term):
                cms.add_message(ConditionMatchingMessage(severity="warning", text=f"Converted from OMIM term {term.id}"))
                term = mondo
            else:
                # slowly direct users to MONDO
                cms.add_message(ConditionMatchingMessage(severity="warning", text="Matched on OMIM, please attempt to find MONDO term if possible"))

        cms.add_term(term)
        cms.validate()
        return cms
    return cms


def merge_matches(matches: list[ConditionMatchingSuggestion]) -> Optional[ConditionMatchingSuggestion]:
    if not matches:
        return None
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        if name_matches := [match for match in matches if match.alias_index is None]:
            if len(name_matches) == 1:
                return name_matches[0]

        if len(matches) == 2:
            # if we matched 2 terms, and one is the child of the other, just return the parent most item
            if is_descendant(terms={matches[0].terms[0]}, ancestors={matches[1].terms[0]}, check_levels=3):
                return matches[1]
            if is_descendant(terms={matches[1].terms[0]}, ancestors={matches[0].terms[0]}, check_levels=3):
                return matches[0]

        # we had multiple matches and couldn't pick a best option
        suggestion = ConditionMatchingSuggestion()
        for match in matches:
            for term in match.terms:
                suggestion.add_term(term)
        suggestion.condition_multi_operation = MultiCondition.UNCERTAIN
        suggestion.add_message(ConditionMatchingMessage(severity="error", text="Text matched multiple terms"))
        suggestion.merged = True
        return suggestion


def find_local_term(match_text: SearchText, service: OntologyService) -> Optional[ConditionMatchingSuggestion]:
    """
    Note that local search is pretty dumb, uses SearchText which already simplifies a bunch of words
    Code looks through local database, erring of false positives, then uses search_text_to_suggestion to
    see if returned result is a genuine positive
    """
    q = []
    # TODO, can we leverage phenotype matching?
    if match_text.prefix_terms:
        term_list = list(match_text.prefix_terms)
        if len(term_list) == 1 and len(term_list[0]) <= 4:
            term_str: str = term_list[0]
            # check array contains (and hope we don't have any mixed case aliases)
            q.append(Q(name__iexact=term_str) | Q(aliases__contains=[term_str.upper()]) | Q(
                aliases__contains=[term_str.lower()]))
        else:
            for term_str in term_list:
                if len(term_str) > 1 and not term_str.isnumeric():
                    # exclude numeric because the value might be stored as roman or arabic
                    # problem with icontains in aliases is it converts array list to a string, and then finds text in there
                    # so "hamper,laundry" would be returned for icontains="ham"
                    q.append(Q(name__icontains=term_str) | Q(aliases__icontains=term_str))

    matches = []
    if q:
        qs = OntologyTerm.objects.filter(ontology_service=service).filter(reduce(operator.and_, q))
        for term in qs[0:200]:
            if not term.is_obsolete:
                if cms := search_text_to_suggestion(match_text, term):
                    matches.append(cms)

    return merge_matches(matches)


def search_suggestion(text: str) -> ConditionMatchingSuggestion:
    """
    embedded ids have already been searched for, check the MONDO search database
    """
    match_text = SearchText(text)
    if local_mondo := find_local_term(match_text, OntologyService.MONDO):
        return local_mondo

    try:
        matches: list[ConditionMatchingSuggestion] = []
        for term in condition_text_search(text):
            if cms := search_text_to_suggestion(match_text, term):
                matches.append(cms)
        if search_match := merge_matches(matches):
            return search_match
    except:
        report_exc_info()

    if local_omim := find_local_term(match_text, OntologyService.OMIM):
        return local_omim

    return ConditionMatchingSuggestion()


def _is_descendant_ids(term_ids: set[int], ancestors_ids: set[int], seen_ids: set[int], check_levels: int):
    # TODO move this to OntologyTerm class
    for term in term_ids:
        if term in ancestors_ids:
            return True
    if check_levels == 0:
        return False

    all_parent_terms = set()
    if term_ids:
        parent_qs = OntologyTermRelation.objects.filter(source_term__pk__in=term_ids, relation=OntologyRelation.IS_A)
        for parent in parent_qs.values_list("dest_term", flat=True):
            if parent not in seen_ids:
                all_parent_terms.add(parent)

    if all_parent_terms:
        seen_ids = seen_ids.union(all_parent_terms)
        return _is_descendant_ids(all_parent_terms, ancestors_ids, seen_ids, check_levels - 1)


def is_descendant(terms: set[OntologyTerm], ancestors: set[OntologyTerm], check_levels: int = 10):
    return _is_descendant_ids({term.id for term in terms}, {term.id for term in ancestors}, set(), check_levels)


def _sort_terms(terms: Iterable[OntologyTerm]) -> list[OntologyTerm]:
    # sort by name (so OMIM, MONDO etc will be sorted together)
    # if we don't have a name, then fall back to ontology service, and then index, so DOID:00034 will appear before DOID:400
    return sorted(list(terms), key=attrgetter("name", "ontology_service", "index"))


def condition_matching_suggestions(ct: ConditionText, ignore_existing=False) -> list[ConditionMatchingSuggestion]:
    suggestions = []

    root_level = ct.root
    root_cms: Optional[ConditionMatchingSuggestion]

    root_cms = ConditionMatchingSuggestion(root_level, ignore_existing=ignore_existing)
    display_root_cms = root_cms
    if not root_cms.terms:
        root_cms = top_level_suggestion(ct.normalized_text)
        root_cms.condition_text_match = root_level
        if root_cms.ids_found_in_text:
            display_root_cms = root_cms
        else:
            display_root_cms = root_cms
            display_root_cms.info_only = True

    root_cms.validate()
    suggestions.append(display_root_cms)

    # filled in and gene level, exclude root as we take care of that before-hand
    filled_in = ct.conditiontextmatch_set.annotate(condition_xrefs_length=ArrayLength('condition_xrefs'))
    filled_in = filled_in.filter(Q(condition_xrefs_length__gt=0) | Q(parent=root_level)).exclude(gene_symbol=None)
    root_level_terms = root_cms.terms
    root_level_mondo: set[OntologyTerm] = {term for term in root_level_terms if term.ontology_service == OntologyService.MONDO}

    for ctm in filled_in.order_by('gene_symbol'):
        if ctm.condition_xref_terms and not ignore_existing:
            cms = ConditionMatchingSuggestion(ctm)
            cms.validate()
            suggestions.append(cms)

        elif ctm.is_gene_level:  # should be the only other option
            # chances are we'll have some suggestions or warnings for gene level if we have root level
            # if we don't no foul, we just send down an empty context
            # and if we got rid of root level, might need to blank out gene level suggestions/warnings
            cms = ConditionMatchingSuggestion(condition_text_match=ctm, ignore_existing=ignore_existing)
            suggestions.append(cms)

            gene_symbol = ctm.gene_symbol
            if root_level_mondo:
                gene_level_terms: set[OntologyTerm] = OntologySnake.mondo_terms_for_gene_symbol(gene_symbol=gene_symbol)
                matches_gene_level: set[OntologyTerm] = set()
                for gene_level in gene_level_terms:
                    if is_descendant({gene_level}, root_level_mondo):
                        matches_gene_level.add(gene_level)

                not_root_gene_terms: list[OntologyTerm] = []
                for term in matches_gene_level:
                    if term not in root_level_terms:
                        not_root_gene_terms.append(term)
                not_root_gene_terms = _sort_terms(not_root_gene_terms)

                matches_gene_level_leafs = [term for term in matches_gene_level if term.is_leaf]
                root_level_str = ', '.join([term.id for term in _sort_terms(root_level_mondo)])

                if not matches_gene_level:
                    cms.add_message(ConditionMatchingMessage(severity="warning", text=f"{root_level_str} : Could not find association to {gene_symbol}"))
                elif len(matches_gene_level) == 1:
                    term = list(matches_gene_level)[0]
                    if term in root_level_terms:
                        cms.add_message(ConditionMatchingMessage(severity="success",
                                                                 text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                        if len(root_level_terms) != 1 and not root_cms.ids_found_in_text:  # if multiple terms from root level
                            # just leave them the same (don't make a suggestion, just validate)
                            cms.add_term(term)
                    else:
                        cms.add_message(ConditionMatchingMessage(severity="success",
                                                                 text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                        cms.add_term(term)  # not guaranteed to be a leaf, but no associations on child terms
                elif len(matches_gene_level_leafs) == 1:
                    term = list(matches_gene_level_leafs)[0]
                    cms.add_message(ConditionMatchingMessage(severity="success",
                                                             text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                    cms.add_term(matches_gene_level_leafs[0])
                    # code that used to tell you which other leaf terms were related, but got too messy
                    # for term in sorted(list(not_root_gene_terms), key=attrgetter("name")):
                    #    if term != matches_gene_level_leafs[0]:
                    #        cms.add_message(ConditionMatchingMessage(severity="info", text=f"{term.id} {term.name} is also associated to {gene_symbol}"))

                elif len(not_root_gene_terms) == 1:
                    term = not_root_gene_terms[0]
                    cms.add_term(term)
                    cms.add_message(ConditionMatchingMessage(severity="success",
                                                             text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                else:
                    cms.add_message(ConditionMatchingMessage(severity="warning", text=f"{root_level_str} : Multiple descendants of this term are associated to {gene_symbol}"))
                    for term in not_root_gene_terms:
                        cms.add_message(ConditionMatchingMessage(severity="info", text=f"{term.id} {term.name} is associated to {gene_symbol}"))

            else:
                # if not MONDO term, see if this term has a known relationship directly
                parent_term_missing_gene = []
                parent_term_has_gene = []
                for term in root_level_terms:
                    if not term.is_stub:
                        if not OntologySnake.has_gene_relationship(term, gene_symbol):
                            parent_term_missing_gene.append(term)
                        else:
                            parent_term_has_gene.append(term)
                parent_term_missing_gene = _sort_terms(parent_term_missing_gene)
                parent_term_has_gene = _sort_terms(parent_term_has_gene)

                for term in parent_term_missing_gene:
                    cms.add_message(ConditionMatchingMessage(severity="warning", text=f"{term.id} : no relationship on file to {gene_symbol.symbol}"))
                for term in parent_term_has_gene:
                    cms.add_message(ConditionMatchingMessage(severity="success",
                                                             text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))

                if cms.terms == root_cms.terms and root_cms.ids_found_in_text:
                    cms.terms = []  # no need to duplicate when ids found in text

            if not cms.terms:
                if root_cms.is_applied:
                    # if parent was applied, and all we have are warnings
                    # put them in the applied column not suggestion
                    cms.is_applied = True
                else:
                    # if root was a suggestion, but we couldn't come up with a more specific suggestion
                    # suggest the root at each gene level anyway (along with any warnings we may have generated)
                    if not root_cms.ids_found_in_text and not root_cms.merged:
                        cms.terms = root_cms.terms  # just copy parent term if we couldn't use the child term
                        # for message in root_cms.messages:
                        #    cms.add_message(message)
    return suggestions


def apply_condition_resolution(classification: Classification, new_condition_resolution: Optional[ConditionResolvedDict] = None):
    with transaction.atomic():
        classification_data = Classification.objects.filter(id=classification.id
                                                            ).select_for_update().values('condition_resolution').first()
        old_condition_resolution = classification_data['condition_resolution']
        condition_text_old = "No Matched Condition"
        condition_text = "No Matched Condition"

        if old_condition_resolution:
            condition_text_old = ConditionResolved.from_dict(old_condition_resolution).summary

        new_condition = None
        if new_condition_resolution:
            new_condition = ConditionResolved.from_dict(new_condition_resolution)
            condition_text = new_condition.summary

        if old_condition_resolution != new_condition_resolution:
            condition_text = f"{condition_text_old} --> {condition_text}"
            classification.condition_resolution = new_condition_resolution

            # TODO, handle this via signal?
            classification.flag_collection_safe.add_flag(
                comment=condition_text,
                data=new_condition_resolution,
                flag_type=classification_flag_types.condition_resolution,
            )
            classification.save(update_fields=['condition_resolution'])

            condition_set_signal.send(sender=Classification, classification=classification,
                                      condition=new_condition)


def apply_condition_resolution_to_classifications(ctm: ConditionTextMatch):
    """
    Updates classification's cached condition_resolution (where appropriate) for the ConditionTextMatch
    Call when the ConditionTextMatch has changed
    """
    if ctm.classification_id:
        classification = ctm.classification
        if matching := condition_text_match_for_classification_id(ctm.classification_id):
            apply_condition_resolution(classification, matching.as_resolved_condition())
        else:
            apply_condition_resolution(classification, None)
    else:
        # we're higher level than a classification now
        lowest_valid = ctm
        while lowest_valid and not lowest_valid.is_valid:
            lowest_valid = lowest_valid.parent
        condition_resolution = lowest_valid.as_resolved_condition() if lowest_valid else None
        # now to find all condition text matches under this one, that don't have an override
        children = ConditionTextMatch.objects.filter(parent=ctm).select_related("classification")
        while children:
            new_children: list[ConditionTextMatch] = []
            for child in children:
                if not child.is_valid:  # if child is valid then it already has a value unaffected by this
                    if classification := child.classification:
                        apply_condition_resolution(classification, condition_resolution)
                    else:
                        new_children.append(child)
            children = ConditionTextMatch.objects.filter(parent__in=new_children)


@receiver(post_save, sender=ConditionTextMatch)
def condition_text_saved(sender, instance: ConditionTextMatch, created: bool, raw, using, update_fields, **kwargs):
    if not created or instance.classification_id:
        # only worth doing this for new classification links or already existing ctms
        # e.g. a newly created root, gene level, inheritance level will not have any terms to assign
        apply_condition_resolution_to_classifications(instance)


def sync_all_condition_resolutions_to_classifications():
    Classification.objects.update(condition_resolution=None)
    ctms = ConditionTextMatch.objects.annotate(condition_xrefs_length=ArrayLength('condition_xrefs')).filter(condition_xrefs_length__gt=0)
    for ctm in ctms:
        if ctm.is_valid:
            apply_condition_resolution_to_classifications(ctm)


def condition_text_match_for_classification_id(cid: int) -> Optional[ConditionTextMatch]:
    try:
        ctm: ConditionTextMatch = ConditionTextMatch.objects.get(classification__id=cid)
        while ctm and not ctm.is_valid:
            ctm = ctm.parent
        return ctm
    except ConditionTextMatch.DoesNotExist:
        return None
