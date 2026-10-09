"""
Phenotype text matched to ontology terms. A PhenotypeDescription is one owner's text (a Patient's or a Cohort's,
through its owner field), split into TextPhenotypeSentences; each points at the TextPhenotype for that sentence,
shared by every description containing it and matched once per PhenotypeMatchVersion into TextPhenotypeMatches.

Entry points: PhenotypeDescription.get_results (the JSON pages highlight from), get_ontology_term_ids,
patient_phenotype_terms / patient_phenotypes_for_samples (many patients in one query). Matching is
patients/phenotype_matching.py.
"""
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Optional

from cache_memoize import cache_memoize
from django.contrib.auth.models import User
from django.db import models
from django.db.models import F, Prefetch, Q, QuerySet
from django.db.models.deletion import CASCADE, SET_NULL
from model_utils.models import TimeStampedModel

from library.constants import DAY_SECS
from ontology.models import OntologyService, OntologyTerm, OntologyVersion
from patients.models.has_phenotype_description_mixin import PHENOTYPE_DESCRIPTION
from patients.models.models_patient import Patient
from patients.phenotype_matcher import PHENOTYPE_MATCHER_VERSION, get_ambiguous_acronym_denylist

# TextPhenotypeMatch back to the patient it was matched for
TPM_PATIENT_PATH = "text_phenotype__textphenotypesentence__phenotype_description__patient"
TPM_DESCRIPTION_TEXT_PATH = "text_phenotype__textphenotypesentence__phenotype_description__original_text"

# The 3 ontologies a phenotype description matches to, labelled as OntologyTerm.split_hpo_omim_mondo_as_dict has them
PHENOTYPE_ONTOLOGY_SERVICE_LABELS = {
    OntologyService.HPO: "HPO",
    OntologyService.OMIM: "OMIM",
    OntologyService.MONDO: "MONDO",
}

AmbiguousAcronymDenylist = Mapping[str, tuple[tuple[str, str], ...]]


def ambiguous_acronym_pattern(denylist: AmbiguousAcronymDenylist) -> Optional[re.Pattern]:
    """ Finds the denylisted acronyms in a sentence - compile once per description, there are a lot of them """
    if not denylist:
        return None
    return re.compile(r"\b(" + "|".join(re.escape(k) for k in denylist) + r")\b", re.IGNORECASE)


def without_ambiguous_acronyms(matches: Iterable['TextPhenotypeMatch'],
                               denylist: AmbiguousAcronymDenylist) -> list['TextPhenotypeMatch']:
    """ Matches (saved or not) whose text is not an ambiguous acronym - one short string mapping to multiple
        distinct concept clusters, which downstream ontology-term queries can't use without picking the wrong one """
    return [m for m in matches if not m.is_ambiguous_acronym(denylist)]


def _ambiguous_alias_result(match_text: str, offset_start: int, offset_end: int,
                            denylist: AmbiguousAcronymDenylist) -> dict:
    """ A warning-only entry: the UI lists the concepts the acronym could be rather than picking one """
    entry = {
        "ambiguous_alias": match_text,
        "offset_start": offset_start,
        "offset_end": offset_end,
    }
    if candidates := denylist.get(TextPhenotypeMatch.ambiguous_acronym_key(match_text)):
        entry["ambiguous_alias_candidates"] = [{"accession": acc, "name": name} for acc, name in candidates]
    return entry


class PhenotypeMatchVersion(TimeStampedModel):
    """ What a sentence's matches were produced with: the matcher code (PHENOTYPE_MATCHER_VERSION) and the ontology
        its lookups were built from. One row per pair, made the first time that pair matches a sentence (#2131).
        Deleting an OntologyVersion cascades here and leaves its sentences unstamped, so awaiting a rematch. """
    matcher_version = models.IntegerField()
    ontology_version = models.ForeignKey(OntologyVersion, null=True, blank=True, on_delete=CASCADE)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["matcher_version", "ontology_version"], nulls_distinct=False,
                                    name="phenotype_match_version_unique"),
        ]

    def __str__(self):
        return f"matcher v{self.matcher_version} / {self.ontology_version}"

    @staticmethod
    def _current_kwargs() -> dict:
        return {"matcher_version": PHENOTYPE_MATCHER_VERSION, "ontology_version": OntologyVersion.latest(validate=False)}

    @staticmethod
    def get_or_create_current() -> 'PhenotypeMatchVersion':
        """ The current matcher against the latest ontology - what a PhenotypeMatcher built now stamps """
        phenotype_match_version, _ = PhenotypeMatchVersion.objects.get_or_create(**PhenotypeMatchVersion._current_kwargs())
        return phenotype_match_version

    @staticmethod
    def current_qs() -> QuerySet['PhenotypeMatchVersion']:
        """ Empty until something has been matched with the current pair - every matched sentence is stale then """
        return PhenotypeMatchVersion.objects.filter(**PhenotypeMatchVersion._current_kwargs())


class TextPhenotype(models.Model):
    """ One row per distinct sentence - the unit of matching, shared by every description containing it """
    text = models.TextField(unique=True)
    # Null = awaiting matching; stamped when matched, set back to null to requeue
    match_version = models.ForeignKey(PhenotypeMatchVersion, null=True, blank=True, on_delete=SET_NULL)

    def __str__(self):
        return f"{self.text} (matched with: {self.match_version})"

    def mark_matched(self, match_version: PhenotypeMatchVersion):
        self.match_version = match_version
        self.save(update_fields=["match_version"])

    @staticmethod
    def awaiting_qs() -> QuerySet['TextPhenotype']:
        """ Never matched, or requeued - bulk_patient_phenotype_matching matches these """
        return TextPhenotype.objects.filter(match_version__isnull=True)

    @staticmethod
    def stale_qs() -> QuerySet['TextPhenotype']:
        """ Matched with an older matcher or ontology - `match_patient_phenotypes --stale` rematches these """
        return (TextPhenotype.objects.filter(match_version__isnull=False)
                .exclude(match_version__in=PhenotypeMatchVersion.current_qs()))


class PhenotypeDescription(models.Model):
    """ A patient's or cohort's phenotype text, split into sentences. Owned by its patient or its cohort, so it
        goes when they do; owned by neither only for the moment a live preview needs one.
        The owner field is named after the owner's model - @see HasPhenotypeDescriptionMixin """
    original_text = models.TextField()
    patient = models.OneToOneField("patients.Patient", null=True, blank=True, related_name=PHENOTYPE_DESCRIPTION,
                                   on_delete=CASCADE)
    cohort = models.OneToOneField("snpdb.Cohort", null=True, blank=True, related_name=PHENOTYPE_DESCRIPTION,
                                  on_delete=CASCADE)
    approved_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)

    class Meta:
        constraints = [
            models.CheckConstraint(condition=Q(patient__isnull=True) | Q(cohort__isnull=True),
                                   name="phenotype_description_one_owner"),
        ]

    def __str__(self):
        owner = self.patient or self.cohort or "Preview"
        return f"{owner}: {self.original_text[:50]}"

    def get_results(self) -> list[dict]:
        """ The matches of every sentence, offset into the whole text - what pages highlight from """
        match_qs = TextPhenotypeMatch.objects.select_related("ontology_term")
        sentences = list(self.textphenotypesentence_set.select_related("text_phenotype")
                         .prefetch_related(Prefetch("text_phenotype__textphenotypematch_set", queryset=match_qs)))
        ontology_version = None
        if any(sentence.text_phenotype.textphenotypematch_set.all() for sentence in sentences):
            ontology_version = OntologyVersion.latest()
        denylist = get_ambiguous_acronym_denylist()
        acronym_pattern = ambiguous_acronym_pattern(denylist)

        results = []
        for sentence in sentences:
            results.extend(sentence.get_results(ontology_version, denylist, acronym_pattern))
        return results

    @cache_memoize(timeout=DAY_SECS, args_rewrite=lambda s: (s.pk, ))
    def get_ontology_term_ids(self) -> list[str]:
        denylist = get_ambiguous_acronym_denylist()
        tpm_qs = (TextPhenotypeMatch.objects
                  .filter(text_phenotype__textphenotypesentence__phenotype_description=self)
                  .select_related("text_phenotype"))
        return sorted({tpm.ontology_term_id for tpm in without_ambiguous_acronyms(tpm_qs, denylist)})


class TextPhenotypeSentence(models.Model):
    """ A description broken up into sentences which are represented by TextPhenotypes and matches """
    phenotype_description = models.ForeignKey(PhenotypeDescription, on_delete=CASCADE)
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    sentence_offset = models.IntegerField()

    def __str__(self):
        return f"{self.sentence_offset}: {self.text_phenotype_id}"

    def get_results(self, ontology_version: Optional[OntologyVersion], denylist: AmbiguousAcronymDenylist,
                    acronym_pattern: Optional[re.Pattern]) -> list[dict]:
        """ ontology_version may be None only when the sentence has no matches.
            Ambiguous acronyms come back as warning-only entries rather than terms """
        all_matches = self.text_phenotype.textphenotypematch_set.all()
        matches = without_ambiguous_acronyms(all_matches, denylist)
        # Where an Ontology Service has multiple matches to the exact same text
        span_counts = Counter(tpm.service_span for tpm in matches)

        results = []
        for tpm in matches:
            data = tpm.to_dict(ontology_version)
            data["offset_start"] += self.sentence_offset
            data["offset_end"] += self.sentence_offset
            if span_counts[tpm.service_span] >= 2:
                data["ambiguous"] = tpm.match_text
            results.append(data)

        # Ambiguous acronyms aren't saved as matches, so find them by rescanning the sentence text
        warned_spans = set()
        if acronym_pattern:
            for m in acronym_pattern.finditer(self.text_phenotype.text):
                warned_spans.add((m.start(), m.end()))
                results.append(_ambiguous_alias_result(m.group(0), m.start() + self.sentence_offset,
                                                       m.end() + self.sentence_offset, denylist))
        # A row matched before its text joined the denylist
        for tpm in all_matches:
            if (tpm.offset_start, tpm.offset_end) not in warned_spans and tpm.is_ambiguous_acronym(denylist):
                warned_spans.add((tpm.offset_start, tpm.offset_end))
                results.append(_ambiguous_alias_result(tpm.match_text, tpm.offset_start + self.sentence_offset,
                                                       tpm.offset_end + self.sentence_offset, denylist))
        return results


class TextPhenotypeMatch(models.Model):
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    offset_start = models.IntegerField()
    offset_end = models.IntegerField()

    def __str__(self):
        return f"{self.ontology_term} from ({self.offset_start}-{self.offset_end})"

    @property
    def match_text(self) -> str:
        txt = self.text_phenotype.text
        return txt[self.offset_start:self.offset_end]

    @property
    def service_span(self) -> tuple:
        return self.ontology_term.ontology_service, self.offset_start, self.offset_end

    @staticmethod
    def ambiguous_acronym_key(match_text: str) -> str:
        """ How match text is looked up in get_ambiguous_acronym_denylist() """
        return match_text.lower().replace(",", "")

    def is_ambiguous_acronym(self, denylist: AmbiguousAcronymDenylist) -> bool:
        return bool(denylist) and self.ambiguous_acronym_key(self.match_text) in denylist

    def to_dict(self, ontology_version: OntologyVersion) -> dict:
        """ This is what's sent as JSON back to client for highlighting and grids """
        # Iterate the memoized QuerySet (its pickle holds the rows) rather than querying it again
        gene_symbols = [gene_symbol.pk for gene_symbol in
                        ontology_version.cached_gene_symbols_for_terms_tuple((self.ontology_term.pk,))]
        accession = str(self.ontology_term)
        return {
            "accession": accession,
            "gene_symbols": gene_symbols,
            "match": accession,
            "ontology_service": self.ontology_term.get_ontology_service_display(),
            "name": self.ontology_term.name,
            "offset_start": self.offset_start,
            "offset_end": self.offset_end,
            "pk": self.ontology_term.pk,
        }


@dataclass(frozen=True)
class PatientPhenotypeTerms:
    """ A patient's phenotype text and the ontology terms matched in it. match_texts is the text each term
        matched (term pk -> text), so a page can collapse the terms one phrase matched into one chip """
    text: str
    terms: dict[str, list[OntologyTerm]]
    match_texts: dict[str, str] = field(default_factory=dict)

    def to_json(self) -> dict:
        terms = {}
        for service_label in PHENOTYPE_ONTOLOGY_SERVICE_LABELS.values():
            terms[service_label] = [{"id": term.pk, "name": term.name, "match_text": self.match_texts.get(term.pk)}
                                    for term in self.terms.get(service_label, [])]
        return {"text": self.text, "terms": terms}


def patient_phenotype_terms(patients: Iterable[Patient]) -> dict[int, PatientPhenotypeTerms]:
    """ Patient pk -> matched terms, for every patient that has phenotype text, in one query.
        The single object path is Patient.get_ontology_term_ids() - the two must agree. """
    denylist = get_ambiguous_acronym_denylist()
    tpm_qs = (TextPhenotypeMatch.objects
              .filter(**{TPM_PATIENT_PATH + "__in": patients})
              .select_related("text_phenotype", "ontology_term")
              .annotate(patient_id=F(TPM_PATIENT_PATH), phenotype_text=F(TPM_DESCRIPTION_TEXT_PATH))
              .order_by("pk"))

    text_by_patient_id = {}
    terms_by_patient_id = defaultdict(set)
    match_texts_by_patient_id = defaultdict(dict)
    for tpm in tpm_qs:
        text_by_patient_id[tpm.patient_id] = tpm.phenotype_text
        if tpm.is_ambiguous_acronym(denylist):
            continue
        terms_by_patient_id[tpm.patient_id].add(tpm.ontology_term)
        match_texts_by_patient_id[tpm.patient_id].setdefault(tpm.ontology_term_id, tpm.match_text)

    phenotype_terms = {}
    for patient_id, text in text_by_patient_id.items():
        terms_by_service = defaultdict(list)
        for term in sorted(terms_by_patient_id[patient_id]):
            if service_label := PHENOTYPE_ONTOLOGY_SERVICE_LABELS.get(term.ontology_service):
                terms_by_service[service_label].append(term)
        phenotype_terms[patient_id] = PatientPhenotypeTerms(text=text, terms=dict(terms_by_service),
                                                            match_texts=match_texts_by_patient_id[patient_id])
    return phenotype_terms


def patient_phenotypes_for_samples(user, samples: Iterable) -> dict[int, dict]:
    """ Patient pk -> the JSON a page draws phenotype chips from, for the patients of these samples
        the user can view - one query for the whole page. A patient the user cannot view is absent. """
    patient_ids = {sample.patient_id for sample in samples if sample.patient_id}
    if not patient_ids:
        return {}
    patients = Patient.filter_for_user(user).filter(pk__in=patient_ids)
    return {patient_id: phenotype_terms.to_json()
            for patient_id, phenotype_terms in patient_phenotype_terms(patients).items()}
