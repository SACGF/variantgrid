"""
A series of models that currently stores the combination of MONDO, OMIM, HPO & HGNC.
(Note that HGNC is included just to model relationships, please use GeneSymbol for all your GeneSymbol needs).
"""
import functools
import logging
import operator
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Optional, Union, Iterable, Any, Iterator

from cache_memoize import cache_memoize
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.db import models, connection
from django.db.models import PROTECT, CASCADE, QuerySet, Q, Max, TextChoices
from django.urls import reverse
from model_utils.models import TimeStampedModel, now
from psqlextra.models import PostgresPartitionedModel
from psqlextra.types import PostgresPartitioningMethod

from genes.models import GeneSymbol
from library.cache import timed_cache
from library.constants import DAY_SECS, WEEK_SECS
from library.log_utils import report_exc_info
from library.preview_request import PreviewData, PreviewModelMixin
from library.utils import Constant, first


class OntologyImportSource:
    PANEL_APP_AU = "PAAU"
    MONDO = "MONDO"
    OMIM = "OMIM"
    HPO = "HP"
    HGNC = "HGNC"
    GENCC = 'gencc'


class OntologyService(models.TextChoices):
    MONDO = "MONDO", "MONDO"
    OMIM = "OMIM", "OMIM"
    HPO = "HP", "HP"
    HGNC = "HGNC", "HGNC"

    #Not stored locally
    DOID = "DOID", "DOID"
    ORPHANET = "ORPHA", "ORPHA"
    MEDGEN = "MedGen", "MedGen"
    MeSH = "MeSH", "MeSH"

    EXPECTED_LENGTHS: dict[str, int] = Constant({
        MONDO[0]: 7,
        OMIM[0]: 6,
        HPO[0]: 7,
        HGNC[0]: 1,  # HGNC ids aren't typically 0 padded, because they're not monsters
        DOID[0]: None,  # variable length with padded 0s
        ORPHANET[0]: 1,  # ORPHANET ids aren't typically 0 padded
        MEDGEN[0]: None,
        MeSH[0]: None
    })

    IMPORTANCE: dict[str, int] = Constant({
        MONDO[0]: 2,
        OMIM[0]: 3,
        HPO[0]: 4,  # put HPO relationships last as they occasionally spam OMIM
        DOID[0]: 5,
        ORPHANET[0]: 6,
        HGNC[0]: 1,  # show gene relationships first
        MEDGEN[0]: 7,
        MeSH[0]: 8
    })

    URLS: dict[str, str] = Constant({
        MONDO[0]: "https://monarchinitiative.org/disease/MONDO:${1}",
        OMIM[0]: "http://www.omim.org/entry/${1}",
        HPO[0]: "https://hpo.jax.org/app/browse/term/HP:${1}",
        HGNC[0]: "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:${1}",
        DOID[0]: "https://www.ebi.ac.uk/ols/ontologies/doid/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FDOID_${1}",
        ORPHANET[0]: "https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=EN&Expert=${1}",
        MEDGEN[0]: "https://www.ncbi.nlm.nih.gov/medgen/${1}",
        MeSH[0]: "https://meshb.nlm.nih.gov/record/ui?ui=${1}"
    })

    LOCAL_ONTOLOGY_PREFIXES: set[str] = Constant({
        MONDO[0],
        OMIM[0],
        HPO[0],
        HGNC[0]
    })

    CONDITION_ONTOLOGIES: set[str] = Constant({
        MONDO[0],
        OMIM[0],
        HPO[0],
        DOID[0],
        ORPHANET[0],
        MEDGEN[0],
        MeSH[0]
    })

    @staticmethod
    def index_to_id(ontology_service: 'OntologyService', index: Union[int|str]):
        num_part = str(index)
        if expected_length := OntologyService.EXPECTED_LENGTHS.get(ontology_service):
            num_part = str(index).rjust(expected_length, '0')
        return f"{ontology_service}:{num_part}"


class OntologyRelation:
    """
    Common relationships, relationship is free text.
    Note it's best to look at the import source for these.
    """
    IS_A = "is_a"
    EXACT = "exact"  # defined by HPO and MONDO
    EXACT_SYNONYM = "exact_synonym"
    RELATED = "related"  # defined by HPO and MONDO (also use related Synonym from mondo to populate this)
    RELATED_SYNONYM = "related_synonym"
    CLOSE = "close"  # defined by HPO and MONDO
    BROAD = "broad"  # defined by HPO and MONDO
    NARROW = "narrow"  # defined by HPO and MONDO
    ALTERNATIVE = "alternative"  # HPO has alternative ID in mondo file
    XREF = "xref"  # listed in MONDO xrefs, probably is the same term
    REPLACED = "replaced"

    CONSIDER = "consider"

    FREQUENCY = "frequency"
    PANEL_APP_AU = "panelappau"
    ASSOCIATED = "associated"  # used by GenCC, phenotypes_to_genes
    ALL_FREQUENCY = "frequency"  # used by OMIM_ALL_FREQUENCIES
    ENTREZ_ASSOCIATION = "associated condition"
    MIM_2_GENE = "mim2gene"

    DISPLAY_NAMES = {
        IS_A: "is a",
        EXACT_SYNONYM: "has an exact synonym from",
        RELATED_SYNONYM: "has a related synonym from",
        ALL_FREQUENCY: "frequently occurs with",
        ENTREZ_ASSOCIATION: "has an associated gene of",
        PANEL_APP_AU: "PanelApp AU association"
    }

    """
    # MONDO associations
    "http://purl.obolibrary.org/obo/RO_0004025": "disease causes dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004001": "has material basis in gain of function germline mutation in",
    "http://purl.obolibrary.org/obo/RO_0004021": "disease has basis in disruption of",
    "http://purl.obolibrary.org/obo/RO_0004020": "disease has basis in dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004026": "disease has location",
    "http://purl.obolibrary.org/obo/RO_0004027": "disease has inflammation site",
    "http://purl.obolibrary.org/obo/RO_0004030": "disease arises from structure"
    """


class PanelAppClassification(models.TextChoices):
    RED = "1", "Expert Review Red"
    AMBER = "2", "Expert Review Amber"
    GREEN = "3", "Expert Review Green"

    @property
    def is_strong_enough(self) -> bool:
        return self == PanelAppClassification.GREEN

    def __lt__(self, other) -> bool:
        return int(self.value) < int(other.value)

    @staticmethod
    def get_by_label_pac(label: str) -> 'PanelAppClassification':
        for pac in PanelAppClassification:
            if pac.label == label:
                return pac
        raise ValueError(f"No PanelAppClassification for {label}")

    @staticmethod
    def get_above_min(min_classification: str) -> set[str]:
        classifications = set()
        for e in reversed(PanelAppClassification):
            classifications.add(e.label)
            if e.value == min_classification:
                break
        return classifications


class GeneDiseaseClassification(models.TextChoices):
    # @see https://thegencc.org/faq.html#validity-termsdelphi-survey - where sort order comes from
    # Not using e.g. "GENCC:100009" (= Supportive) as that has different sort order than page above
    REFUTED = "1", "Refuted Evidence"
    NO_KNOWN = "2", "No Known Disease Relationship"
    ANIMAL = "3", "Animal Model Only"
    DISPUTED = "4", "Disputed Evidence"
    LIMITED = "5", "Limited"
    SUPPORTIVE = "6", "Supportive"
    MODERATE = "7", "Moderate"
    STRONG = "8", "Strong"
    DEFINITIVE = "9", "Definitive"

    @property
    def is_strong_enough(self) -> bool:
        return self in {GeneDiseaseClassification.STRONG, GeneDiseaseClassification.DEFINITIVE}

    @staticmethod
    def get_by_label(label: str) -> 'GeneDiseaseClassification':
        for gdc in GeneDiseaseClassification:
            if gdc.label == label:
                return gdc
        raise ValueError(f"No GeneDiseaseClassification for {label}")

    @staticmethod
    def get_above_min(min_classification: str) -> set[str]:
        classifications = set()
        for e in reversed(GeneDiseaseClassification):
            classifications.add(e.label)
            if e.value == min_classification:
                break
        return classifications


class OntologyImport(TimeStampedModel):
    """
    Keeps track of when data was imported, typically used to see how old the data is and if it needs
    to be imported again
    """
    import_source = models.TextField()
    filename = models.TextField()
    version = models.IntegerField()
    context = models.TextField()
    hash = models.TextField()
    processor_version = models.IntegerField(default=1)
    processed_date = models.DateTimeField(auto_created=True)
    completed = models.BooleanField(default=False)

    class Meta:
        unique_together = ('import_source', 'filename', 'version')

    def save(self, *args, **kwargs):
        created = not self.pk
        if created and not self.version:
            # Assign version to be next highest
            existing_ontology = OntologyImport.objects.filter(import_source=self.import_source,
                                                              filename=self.filename)
            data = existing_ontology.aggregate(Max("version"))
            self.version = (data.get("version__max") or 0) + 1

        super().save(*args, **kwargs)
        if created and OntologyVersion.in_ontology_version(self):
            logging.info("Creating new OntologyTermRelation partition")
            import_source_id = self.pk
            connection.schema_editor().add_list_partition(
                model=OntologyTermRelation,
                name=f"import_source_{import_source_id}",
                values=[import_source_id],
            )

    def __str__(self):
        name = f"OntologyImport ({self.pk}) - {self.import_source}: {self.filename} ({self.created.date()})"
        if not self.completed:
            name += " (incomplete)"
        return name


class OntologyTermStatus(TextChoices):
    CONDITION = 'C'  # Also phenotypes are mixed up in there right now
    DEPRECATED = 'D'
    NON_CONDITION = 'N'
    STUB = 'S'


@dataclass
class OntologyIdNormalized:
    prefix: str
    postfix: str
    full_id: str
    clean: bool

    @property
    def num_part(self) -> int:
        return int(self.postfix)

    def num_part_safe(self) -> int:
        try:
            return self.num_part
        except:
            return 0

    @staticmethod
    def normalize(dirty_id: str) -> 'OntologyIdNormalized':
        parts = re.split("[:|_]", dirty_id)
        if len(parts) != 2:
            raise ValueError(f"Can not convert {dirty_id} to a proper id")

        prefix = parts[0].strip().upper()
        postfix = parts[1].strip()

        if prefix in ("ORPHA", "ORPHANET"):  # Orphanet is the one ontology (so far) where the standard is sentance case
            prefix = "ORPHA"
        elif prefix.upper() == "MIM":
            prefix = "OMIM"
        elif prefix.upper() == "MEDGEN":
            prefix = "MedGen"
            postfix = postfix.upper()
        elif prefix.upper() == "MESH":
            prefix = "MeSH"
            postfix = postfix.upper()
        prefix = OntologyService(prefix)

        try:
            if expected_length := OntologyService.EXPECTED_LENGTHS[prefix]:
                postfix = str(int(postfix))  # strip leading 0s, so we can then add the correct number
                postfix = postfix.rjust(expected_length, '0')
            clean_id = f"{prefix}:{postfix}"
            return OntologyIdNormalized(prefix=prefix, postfix=postfix, full_id=clean_id, clean=True)

        except ValueError:
            return OntologyIdNormalized(prefix=prefix, postfix=postfix, full_id=dirty_id, clean=False)

    def __str__(self):
        return self.full_id


class OntologyTerm(TimeStampedModel, PreviewModelMixin):

    """
    id is Term as it should be referenced <prefix>:<zero padded index> e.g.
    MONDO:0000043, OMIM:0000343
    """
    id = models.TextField(primary_key=True)
    ontology_service = models.CharField(max_length=10, choices=OntologyService.choices)
    index = models.IntegerField()
    name = models.TextField(null=True, blank=True)  # should only be null if we're using it as a placeholder reference
    definition = models.TextField(null=True, blank=True)
    extra = models.JSONField(null=True, blank=True)
    aliases = ArrayField(models.TextField(blank=False), null=False, blank=True, default=list)
    status = models.CharField(max_length=1, default=OntologyTermStatus.CONDITION, choices=OntologyTermStatus.choices)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)

    move_to_re = re.compile(r"MOVED TO (\d+)")

    @property
    def oxo_suffix(self):
        if self.ontology_service == OntologyService.MEDGEN:
            lower_id = self.id.lower()
            if ":c" in lower_id:
                return "UMLS:" + lower_id[lower_id.index(":") + 1:].upper()
        elif self.ontology_service == OntologyService.ORPHANET:
            return f"orphanet:{self.index}"
        return self.id

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-disease"

    @property
    def preview(self) -> PreviewData:
        name: str
        if self.is_stub:
            name = ""
        elif not self.name:
            name = "No Name Provided"
        else:
            name = self.name
        return self.preview_with(
            title=name,
            summary=self.definition
        )

    @property
    def is_locally_stored(self):
        return self.ontology_service in OntologyService.LOCAL_ONTOLOGY_PREFIXES

    def __str__(self):
        return f"{self.id} {self.name}"

    @property
    def short(self) -> str:
        if self.ontology_service == OntologyService.HGNC:
            return self.name
        else:
            return self.id

    class Meta:
        unique_together = ("ontology_service", "index")

    def __lt__(self, other):
        if self.ontology_service != other.ontology_service:
            return self.ontology_service < other.ontology_service
        if self.ontology_service == OntologyService.HGNC:
            return self.name < other.name
        return self.index < other.index

    @property
    def is_stub(self):
        return self._state.adding

    @property
    def is_obsolete(self) -> bool:
        # obsolete covers deprecated, moved, and gene (not the most accurate title)
        # more so "is not valid for condition"
        return self.status != OntologyTermStatus.CONDITION

    @property
    def is_valid_for_condition(self) -> bool:
        return self.status == OntologyTermStatus.CONDITION

    @property
    def warning_text(self) -> Optional[str]:
        if self.status == OntologyTermStatus.DEPRECATED:
            return "Term is Deprecated"
        elif self.status == OntologyTermStatus.NON_CONDITION:
            term_type = (self.extra or {}).get('type', 'Unknown')
            return f"Term is of type - {term_type}"
        elif self.status == OntologyTermStatus.STUB:
            return "Term was referenced by 3rd party but not yet from our authoritative source"
        else:
            return None

    @cached_property
    def is_leaf(self) -> bool:
        # Warning, just meant to be called on MONDO terms
        if not self.is_stub and self.ontology_service == OntologyService.MONDO:
            return not OntologyVersion.get_latest_and_live_ontology_qs().filter(dest_term=self, relation=OntologyRelation.IS_A).exists()
        return True

    @property
    def url_safe_id(self):
        return self.id.replace(":", "_")

    @staticmethod
    def get_from_slug(slug_pk):
        pk = slug_pk.replace("_", ":")
        return OntologyTerm.objects.get(pk=pk)

    @staticmethod
    def get_gene_symbol(gene_symbol: Union[str, GeneSymbol]) -> 'OntologyTerm':
        """
        Returns the OntologyTerm for a GeneSymbol (which exists just to make relationship nodes)
        This should be used as little as possible and only to bridge between GeneSymbols and OntologyTerms
        aka if you're talking about GeneSymbols, pass around the GeneSymbol object, not the OntologyTerm
        """
        from genes.gene_matching import HGNCMatcher

        if isinstance(gene_symbol, GeneSymbol):
            gene_symbol = gene_symbol.symbol

        if gene_ontology := OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC, name=gene_symbol).first():
            return gene_ontology
        if gene_ontology := OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC, aliases__contains=[gene_symbol]).first():
            return gene_ontology

        hgnc_matcher = HGNCMatcher.instance()
        if hgnc := hgnc_matcher.match_hgnc(gene_symbol):
            if hgnc_term := OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC, index=hgnc.id).first():
                # we found an Ontology Term for HGNC, but the ID is already in use for another sybmol
                # prioritise the name as defined by HGNCGeneNames but keep the other names in the aliases
                # Also only lazily load aliases, we shouldn't have too many sources of data that request based on
                # deprecated gene names
                hgnc_term.aliases.append(gene_symbol)
                hgnc_term.save()
                return hgnc_term
            # every term needs an import
            o_import = OntologyImport.objects.create(
                import_source=OntologyService.HGNC,
                filename="HGNC Aliases",
                context="adhoc_hgnc",
                hash="N/A",
                processor_version=1,
                processed_date=now,
                completed=True)

            term = OntologyTerm.objects.create(
                id=f"HGNC:{hgnc.id}",
                ontology_service=OntologyService.HGNC,
                index=hgnc.id,
                name=hgnc.gene_symbol_id,
                definition=hgnc.approved_name,
                from_import=o_import
            )
            return term

        raise ValueError(f"Cannot find HGNC for {gene_symbol}")

    @staticmethod
    @timed_cache(size_limit=30, ttl=60)
    def get_or_stub_cached(id_str: str) -> 'OntologyTerm':
        """
        Call this when you're in a context that's not currently loading OntologyTerms into the database
        """
        return OntologyTerm.get_or_stub(id_str)

    @staticmethod
    def get_or_stub(id_str: Union[str, OntologyIdNormalized]) -> 'OntologyTerm':
        """
        Returns an OntologyTerm for the given ID.
        If the OntologyTerm doesn't exist in the database, will create an OntologyTerm but
        WILL NOT persist it to the database
        """
        if isinstance(id_str, OntologyIdNormalized):
            normal_id = id_str
        else:
            normal_id = OntologyIdNormalized.normalize(id_str)

        if normal_id.clean:
            if existing := OntologyTerm.objects.filter(id=normal_id.full_id).first():
                return existing
            try:
                index_num_part_value = normal_id.num_part
            except:
                index_num_part_value = normal_id.num_part_safe  # Ontologies like MedGen can have alpha characters in the "index", providing an index of 0 until we update the model
            return OntologyTerm(
                id=normal_id.full_id,
                ontology_service=normal_id.prefix,
                index=index_num_part_value,
                name=""
            )
        else:
            if existing := OntologyTerm.objects.filter(ontology_service=normal_id.prefix, name__iexact=normal_id.postfix).first():
                return existing
            else:
                raise ValueError(f"Can not convert {id_str} to a proper id")

    @property
    def padded_index(self) -> str:
        # ID should already be padded to right number of indexes
        return self.id.split(":")[1]

    @property
    def external_url(self):
        return OntologyService.URLS[self.ontology_service].replace("${1}", self.padded_index)

    def gene_symbol_url(self) -> Optional[str]:
        """ Only for HGNC, get link to gene_symbol """
        url = None
        if self.ontology_service == OntologyService.HGNC:
            if gene_symbol := GeneSymbol.cast(self.name):
                url = gene_symbol.get_absolute_url()
        return url

    def get_absolute_url(self):
        if url := self.gene_symbol_url():
            return url
        return reverse('ontology_term', kwargs={"term": self.url_safe_id})

    @property
    def best_url(self):
        if self.is_locally_stored:
            return self.get_absolute_url()
        else:
            return self.external_url

    @property
    def url(self):
        if settings.ONTOLOGY_EXTERNAL_LINKS:
            return self.external_url
        else:
            return self.get_absolute_url()

    @property
    def moved_to(self) -> Optional[str]:
        if self.name:
            if match := OntologyTerm.move_to_re.match(self.name):
                return match.group(1)
        return None

    @staticmethod
    def split_hpo_omim_mondo(ontology_term_ids: Iterable[str]) -> tuple[QuerySet, QuerySet, QuerySet]:
        hpo_qs = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.HPO)
        omim_qs = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.OMIM)
        mondo_qs = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.MONDO)
        return hpo_qs, omim_qs, mondo_qs

    @staticmethod
    def split_hpo_omim_mondo_as_dict(ontology_term_ids: Iterable[str]) -> dict[str, QuerySet]:
        hpo_qs, omim_qs, mondo_qs = OntologyTerm.split_hpo_omim_mondo(ontology_term_ids)
        return {"HPO": hpo_qs, "OMIM": omim_qs, "MONDO": mondo_qs}


class OntologyTermRelationManager(models.Manager):

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.select_related("source_term", "dest_term", "from_import")


class OntologyTermRelation(PostgresPartitionedModel, TimeStampedModel):
    """
    Relationship between two terms, is generally considered to be bi-directional (or at the very least
    code typically checks relationships in both directions)

    I haven't elected to use django_dag node_factory here as it only allows one kind of relationship,
    and we have quite a lot.
    """
    objects = OntologyTermRelationManager()
    source_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE, related_name="subject")
    dest_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    relation = models.TextField()
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)

    class PartitioningMeta:
        method = PostgresPartitioningMethod.LIST
        key = ["from_import_id"]

    class Meta:
        unique_together = ("from_import", "source_term", "dest_term", "relation")

    def __str__(self):
        name = f"{self.source_term} -> ({self.relation}) -> {self.dest_term}"
        # Add extra info for gene/disease
        if self.relation == OntologyRelation.RELATED and "strongest_classification" in self.extra:
            moi_classifications = self.get_gene_disease_moi_classifications()
            all_classifications = reversed(GeneDiseaseClassification.labels)
            name += ". Gene/Disease: " + ", ".join(self.get_moi_summary(moi_classifications, all_classifications))
        return name

    def __lt__(self, other):
        return self.source_term < other.source_term

    def other_end(self, term: OntologyTerm) -> OntologyTerm:
        """
        Given a relationship is between two terms, return the other one
        If term is neither of the ends, throw an Error
        """
        if term == self.source_term:
            return self.dest_term
        if term == self.dest_term:
            return self.source_term
        raise ValueError("Term was neither a source or dest term, can't find the other end")

    @property
    def relation_display(self):
        return OntologyRelation.DISPLAY_NAMES.get(self.relation, self.relation)

    @property
    def relationship_quality(self) -> GeneDiseaseClassification | PanelAppClassification | None:
        if extra := self.extra:
            if strongest := extra.get('strongest_classification'):
                try:
                    label = GeneDiseaseClassification.get_by_label(strongest)
                except ValueError:
                    label = PanelAppClassification.get_by_label_pac(strongest)
                return label
        return None

    @staticmethod
    def _as_ontology(term: OntologyTerm, service: OntologyService) -> Optional[OntologyTerm]:
        if term.ontology_service == service:
            return term

        q_dest_service = Q(source_term=term) & Q(dest_term__ontology_service=service)
        q_source_service = Q(dest_term=term) & Q(source_term__ontology_service=service)
        otr_qs = OntologyVersion.get_latest_and_live_ontology_qs().filter(q_dest_service | q_source_service, relation=OntologyRelation.EXACT)
        if mondo_rel := otr_qs.first():
            return mondo_rel.other_end(term)
        return None

    @staticmethod
    def as_mondo(term: OntologyTerm) -> Optional[OntologyTerm]:
        return OntologyTermRelation._as_ontology(term, OntologyService.MONDO)

    @staticmethod
    def as_omim(term: OntologyTerm) -> Optional[OntologyTerm]:
        return OntologyTermRelation._as_ontology(term, OntologyService.OMIM)

    @staticmethod
    def relations_of(term: OntologyTerm, otr_qs: Optional[QuerySet['OntologyTermRelation']] = None) -> list['OntologyTermRelation']:
        def sort_relationships(rel1, rel2):
            other1 = rel1.other_end(term)
            other2 = rel2.other_end(term)
            if other1.ontology_service != other2.ontology_service:
                return OntologyService.IMPORTANCE[other1.ontology_service] - OntologyService.IMPORTANCE[other2.ontology_service]
            rel1source = rel1.source_term_id == term.id
            rel2source = rel2.source_term_id == term.id
            if rel1source != rel2source:
                return -1 if rel1source else 1
            return -1 if other1.index < other2.index else 1

        if otr_qs is None:
            otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()

        items = list(otr_qs.filter(Q(source_term=term) | Q(dest_term=term)))
        items.sort(key=functools.cmp_to_key(sort_relationships))
        return items

    # Gene / Disease classifications
    def get_gene_disease_moi_classifications(self):
        sources = self.extra.get("sources")
        if not sources:
            raise ValueError("Extra does not contain 'sources' - only call this on gene/disease classifications")

        moi_classifications = defaultdict(lambda: defaultdict(set))
        for source in sources:
            moi = source["mode_of_inheritance"]
            classification = source["gencc_classification"]
            submitter = source["submitter"]
            moi_classifications[moi][classification].add(submitter)
        return moi_classifications

    @staticmethod
    def get_moi_summary(moi_classifications, valid_classifications) -> list[str]:
        moi_summary = []
        for moi, classifications in moi_classifications.items():
            classification_submitters = []
            for classification in valid_classifications:
                if submitters := classifications.get(classification):
                    classification_submitters.append(f"{classification}: {'/'.join(sorted(submitters))}")
            if classification_submitters:
                moi_summary.append(f"{moi} ({' '.join(classification_submitters)})")
        return moi_summary


OntologyList = Optional[Union[QuerySet, list[OntologyTerm]]]


@dataclass
class OntologyRelationshipQualityFilter:

    min_gencc_strength: GeneDiseaseClassification
    min_panelapp_strength: PanelAppClassification

    @cached_property
    def filter_q(self) -> Q:
        gencc_q = Q(from_import__import_source=OntologyImportSource.GENCC) & Q(extra__strongest_classification__in=GeneDiseaseClassification.get_above_min(self.min_gencc_strength))
        panel_app_q = Q(from_import__import_source=OntologyImportSource.PANEL_APP_AU) & Q(extra__strongest_classification__in=PanelAppClassification.get_above_min(self.min_panelapp_strength))
        other_import_source_q = ~Q(from_import__import_source__in={OntologyImportSource.GENCC, OntologyImportSource.PANEL_APP_AU})

        return gencc_q | panel_app_q | other_import_source_q


ONTOLOGY_RELATIONSHIP_NO_QUALITY_FILTER = OntologyRelationshipQualityFilter(min_gencc_strength=GeneDiseaseClassification.DISPUTED, min_panelapp_strength=PanelAppClassification.RED)
ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER = OntologyRelationshipQualityFilter(min_gencc_strength=GeneDiseaseClassification.ANIMAL, min_panelapp_strength=PanelAppClassification.RED)
ONTOLOGY_RELATIONSHIP_MEDIUM_QUALITY_FILTER = OntologyRelationshipQualityFilter(min_gencc_strength=GeneDiseaseClassification.LIMITED, min_panelapp_strength=PanelAppClassification.AMBER)
ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER = OntologyRelationshipQualityFilter(min_gencc_strength=GeneDiseaseClassification.STRONG, min_panelapp_strength=PanelAppClassification.GREEN)


class OntologyVersion(TimeStampedModel):
    """ This is used by annotation.AnnotationVersion to keep track of different OntologyImports,
        so we can load historical versions of OntologyTermRelation """

    gencc_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="gencc_ontology_version")
    mondo_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="mondo_ontology_version")
    hp_owl_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="hp_owl_ontology_version")
    hp_phenotype_to_genes_import = models.ForeignKey(OntologyImport, on_delete=PROTECT,
                                                     related_name="hp_phenotype_to_genes_ontology_version")
    omim_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="omim_ontology_version", null=True, blank=True)

    class Meta:
        unique_together = (
            'gencc_import',
            'mondo_import',
            'hp_owl_import',
            'hp_phenotype_to_genes_import',
            'omim_import'  # warning, nulls screw up unique together
        )

    OPTIONAL_IMPORTS = {
        'omim_import',
    }

    ONTOLOGY_IMPORTS = {
        "gencc_import": (OntologyImportSource.GENCC,
                         ['https://search.thegencc.org/download/action/submissions-export-csv',
                          'gencc-submissions.csv']),
        "mondo_import": (OntologyImportSource.MONDO, ['mondo.json']),
        "hp_owl_import": (OntologyImportSource.HPO, ['hp.owl']),
        "hp_phenotype_to_genes_import": (OntologyImportSource.HPO,
                                         ['phenotype_to_genes.txt',
                                          'OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt']),
        "omim_import": (OntologyImportSource.OMIM, ['mimTitles.txt'])
    }

    @staticmethod
    def in_ontology_version(ontology_import: OntologyImport) -> bool:
        versioned = defaultdict(set)
        for import_source, filenames in OntologyVersion.ONTOLOGY_IMPORTS.values():
            versioned[import_source].update(filenames)
        return ontology_import.filename in versioned[ontology_import.import_source]

    @staticmethod
    def latest(validate=True) -> Optional['OntologyVersion']:
        oi_qs = OntologyImport.objects.all()
        kwargs = {}
        missing_fields = set()
        for field, (import_source, filenames) in OntologyVersion.ONTOLOGY_IMPORTS.items():
            if ont_import := oi_qs.filter(import_source=import_source, filename__in=filenames).order_by("pk").last():
                kwargs[field] = ont_import
            elif field not in OntologyVersion.OPTIONAL_IMPORTS:
                missing_fields.add(field)

        if not missing_fields:
            values = list(kwargs.values())
            last_date = max(oi.created for oi in values)
            ontology_version, created = OntologyVersion.objects.get_or_create(**kwargs,
                                                                              defaults={"created": last_date})
            if created:
                # Avoid circular import
                from annotation.models import AnnotationVersion
                AnnotationVersion.new_sub_version(None)
        else:
            if validate:
                msg = "OntologyVersion.latest() - missing fields: %s", ", ".join(missing_fields)
                raise OntologyVersion.DoesNotExist(msg)
            else:
                ontology_version = None
        return ontology_version

    def get_ontology_imports(self):
        return [ont_import for ont_import in [
            self.gencc_import,
            self.mondo_import,
            self.hp_owl_import,
            self.hp_phenotype_to_genes_import,
            self.omim_import
        ] if ont_import is not None]

    def get_ontology_term_relations(self):
        return OntologyTermRelation.objects.filter(from_import__in=self.get_ontology_imports())

    @staticmethod
    @timed_cache(ttl=60, quick_key_access=True)
    def get_latest_and_live_ontology_qs() -> QuerySet[OntologyTermRelation]:
        latest = OntologyVersion.latest()
        # live relationships of panelappau aren't versioned
        # TODO could restrict only if we have live enabled in settings
        return OntologyTermRelation.objects.filter(Q(from_import__in=latest.get_ontology_imports()) | Q(relation='panelappau'))

    def get_gene_disease_relations_qs(self) -> QuerySet:
        return self.get_ontology_term_relations().filter(relation=OntologyRelation.RELATED,
                                                         extra__strongest_classification__isnull=False)

    @cache_memoize(WEEK_SECS)
    def moi_and_submitters(self) -> tuple[list[str], list[str]]:
        """ Cached lists of MOI/Submitters from GenCC gene/disease extra JSON """
        moi = set()
        submitters = set()
        for extra in self.get_gene_disease_relations_qs().values_list("extra", flat=True):
            for source in extra["sources"]:
                moi.add(source["mode_of_inheritance"])
                submitters.add(source["submitter"])
        return list(sorted(moi)), list(sorted(submitters))

    @cache_memoize(DAY_SECS)
    def cached_gene_symbols_for_terms_tuple(self, terms_tuple: tuple[int]) -> QuerySet:
        """ Slightly restricted signature so we can cache it """
        return self.gene_symbols_for_terms(terms_tuple)

    def gene_symbols_for_terms(self, terms: OntologyList) -> QuerySet:
        """ This is uncached, see also: cached_gene_symbols_for_terms """
        gene_symbol_names = set()
        otr_qs = self.get_ontology_term_relations()
        for term in terms:
            if isinstance(term, str):
                term = OntologyTerm.get_or_stub(term)
                if term.is_stub:
                    return GeneSymbol.objects.none()
            gene_symbol_snakes = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC, otr_qs=otr_qs)
            gene_symbol_names.update([snake.leaf_term.name for snake in gene_symbol_snakes])
        return GeneSymbol.objects.filter(symbol__in=gene_symbol_names)

    def terms_for_gene_symbol(self, gene_symbol: Union[str, GeneSymbol], desired_ontology: OntologyService,
                              max_depth=1, quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER) -> 'OntologySnakes':
        otr_qs = self.get_ontology_term_relations()
        return OntologySnake.terms_for_gene_symbol(gene_symbol, desired_ontology, max_depth=max_depth,
                                                   quality_filter=quality_filter, otr_qs=otr_qs)

    def gene_disease_relations(self, gene_symbol: Union[str, GeneSymbol],
                               quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER) -> list[OntologyTermRelation]:
        snake = self.terms_for_gene_symbol(gene_symbol, OntologyService.MONDO,
                                           max_depth=0, quality_filter=quality_filter)
        return snake.leaf_relations(ontology_relation=OntologyRelation.RELATED)

    def __str__(self):
        date_str = self.created.strftime("%d %B %Y")
        return f"v{self.pk}. ({date_str})"


@dataclass
class OntologySnakeStep:
    """
    A step in the snake, used for showing the path
    """
    relation: OntologyTermRelation
    dest_term: OntologyTerm
    reversed: bool = False

    @property
    def relationship(self) -> OntologyTermRelation:
        # relationship is the preferred term, add this property to help migrate over to better wording
        return self.relation

    @property
    def source_term(self) -> OntologyTerm:
        return self.relation.other_end(self.dest_term)


@dataclass
class OntologyTermDescendant:
    term: OntologyTerm
    depth: int

    @property
    def sort_key(self) -> Any:
        return self.depth, self.term

    def __lt__(self, other):
        return self.sort_key < other.sort_key


class OntologySnake:
    """
    Use to "Snake" through Ontology nodes, typically to resolve to/from gene symbols.
    An OntologySnake is meant to be immutable, creating new snake_steps each time snake_step is called
    """

    def __init__(self, source_term: OntologyTerm, leaf_term: Optional[OntologyTerm] = None,
                 paths: Optional[list[OntologyTermRelation]] = None):
        self.source_term = source_term
        self.leaf_term = leaf_term or source_term
        self.paths = paths or []

    @property
    def is_strong_enough(self) -> bool:
        for path in self.paths:
            if gencc_quality := path.relationship_quality:
                if not gencc_quality.is_strong_enough:
                    return False
        return True

    def snake_step(self, relationship: OntologyTermRelation) -> 'OntologySnake':
        """
        Creates a new OntologySnake with this extra relationship
        """
        new_leaf = relationship.other_end(self.leaf_term)
        new_paths = list(self.paths)
        new_paths.append(relationship)
        return OntologySnake(source_term=self.source_term, leaf_term=new_leaf, paths=new_paths)

    def show_steps(self) -> list[OntologySnakeStep]:
        steps: list[OntologySnakeStep] = []
        node = self.source_term
        for path in self.paths:
            node = path.other_end(node)
            relationship_reversed = path.source_term == node
            steps.append(OntologySnakeStep(relation=path, dest_term=node, reversed=relationship_reversed))
        return steps

    def reverse(self) -> 'OntologySnake':
        # simply reversing leaf and source will reverse the direction of all the relationships inside
        return OntologySnake(source_term=self.leaf_term, leaf_term=self.source_term, paths=list(reversed(self.paths)))

    def __repr__(self):
        text = f"{self.source_term}"
        for step in self.show_steps():
            forwards = step.relation.dest_term == step.dest_term
            text += f" {'<' if not forwards else ''}-{step.relation.relation}-{'>' if forwards else ''} {step.dest_term}"
        return text

    @cached_property
    def _sort_key(self):
        if steps := self.show_steps():
            last_step = steps[0]
            return (last_step.relation.from_import.import_source or '').lower(), len(steps)
        else:
            return "", 0

    def __lt__(self, other: 'OntologySnake'):
        return self._sort_key < other._sort_key

    @property
    def leaf_relationship(self) -> OntologyTermRelation:
        return self.paths[-1]

    @cached_property
    def start_source(self) -> OntologyImportSource:
        return self.show_steps()[0].relation.from_import.import_source

    @cached_property
    def get_import_relations(self) -> Optional[OntologyTermRelation]:
        for step in self.show_steps():
            if step.relation.from_import.import_source in {OntologyImportSource.PANEL_APP_AU,
                                                           OntologyImportSource.GENCC,
                                                           OntologyImportSource.MONDO}:
                return step.relation

    @staticmethod
    def check_if_ancestor(descendant: OntologyTerm, ancestor: OntologyTerm, max_levels=4) -> list['OntologySnake']:
        if ancestor == descendant:
            return OntologySnake(source_term=ancestor, leaf_term=descendant)

        if descendant.ontology_service != ancestor.ontology_service:
            raise ValueError(f"Can only check for ancestry within the same ontology service, not {descendant.ontology_service} vs {ancestor.ontology_service}")

        seen: set[OntologyTerm] = {descendant}
        new_snakes: list[OntologySnake] = list([OntologySnake(source_term=descendant)])
        valid_snakes: list[OntologySnake] = []
        level = 0
        while new_snakes:
            level += 1
            snakes_by_leaf: dict[OntologyTerm, OntologySnake] = {}
            for snake in new_snakes:
                snakes_by_leaf[snake.leaf_term] = snake

            new_snakes: list[OntologySnake] = []
            for relationship in OntologyTermRelation.objects\
                    .filter(source_term__in=snakes_by_leaf.keys(), relation=OntologyRelation.IS_A)\
                    .exclude(dest_term__in=seen):
                snake = snakes_by_leaf[relationship.source_term]
                new_snake = snake.snake_step(relationship)
                leaf_term = new_snake.leaf_term
                seen.add(new_snake.leaf_term)
                if leaf_term == ancestor:
                    valid_snakes.append(new_snake)
                else:
                    new_snakes.append(new_snake)
            if valid_snakes:
                return valid_snakes
            if level >= max_levels:
                return []
        return []

    @staticmethod
    def all_descendants_of(term: OntologyTerm, limit: int = 100, otr_qs: QuerySet[OntologyTermRelation] = None) -> tuple[list[OntologyTermDescendant], bool]:
        """
        :param term: The term to find all descendants of
        :param limit: The maximum number of results returned
        :param otr_qs: A query limiting relationships, if not provided will use latest OntologyVersion
        :return: A list of unique descendants and a boolean indicating True = More results but truncated
        """
        if otr_qs is None:
            otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()

        review_terms: set[OntologyTerm] = {term}
        reviewed_terms: set[OntologyTerm] = {term}
        results: list[OntologyTermDescendant] = []

        depth = 1
        while review_terms:
            next_level = OntologyVersion.get_latest_and_live_ontology_qs().filter(dest_term__in=review_terms, relation=OntologyRelation.IS_A)
            review_terms = set()

            for child_relationship in next_level:
                child = child_relationship.source_term
                if child not in reviewed_terms:
                    results.append(OntologyTermDescendant(child, depth))
                    if len(results) >= limit:
                        return results, True

                    reviewed_terms.add(child)  # don't look at this term again in case of multiple inheritance (if some other branch has it as a child again)
                    review_terms.add(child)  # but this should be the first time we're looking at it, so see if it has any children

            depth += 1

        results.sort()
        return results, False

    # TODO only allow EXACT between two anythings that aren't Gene Symbols
    @staticmethod
    def snake_from(term: OntologyTerm, to_ontology: OntologyService,
                   quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
                   max_depth: int = 1, otr_qs: QuerySet[OntologyTermRelation] = None) -> 'OntologySnakes':
        """
        Returns the smallest snake/paths from source term to the desired OntologyService
        Ignores IS_A paths
        """
        if term.ontology_service == to_ontology:
            return OntologySnakes([OntologySnake(source_term=term)])

        if otr_qs is None:
            otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()

        seen: set[OntologyTerm] = set()
        seen.add(term)
        new_snakes: list[OntologySnake] = list([OntologySnake(source_term=term)])
        valid_snakes: list[OntologySnake] = []

        relation_q_list = [
            # the list of relationships below is hardly complete for stopping MONDO <-> OMIM, that's done as an extra step
            # but filter out the most common ones here (and IS_A as we don't want to go up/down the hierarchy)
            ~Q(relation__in={OntologyRelation.IS_A, OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM, OntologyRelation.XREF}),
            quality_filter.filter_q
        ]
        q_relation = functools.reduce(operator.and_, relation_q_list)

        iteration = -1
        while new_snakes:
            iteration += 1
            snakes: list[OntologySnake] = list(new_snakes)
            new_snakes: list[OntologySnake] = []
            by_leafs: dict[OntologyTerm, OntologySnake] = {}
            for snake in snakes:
                if existing := by_leafs.get(snake.leaf_term):
                    if len(snake.paths) < len(existing.paths):
                        by_leafs[snake.leaf_term] = snake
                else:
                    by_leafs[snake.leaf_term] = snake
            all_leafs = by_leafs.keys()

            outgoing = otr_qs.filter(source_term__in=all_leafs).exclude(dest_term__in=seen).filter(q_relation)
            incoming = otr_qs.filter(dest_term__in=all_leafs).exclude(source_term__in=seen).filter(q_relation)
            if to_ontology == OntologyService.HGNC:
                outgoing = outgoing.exclude(dest_term__ontology_service=OntologyService.HPO)
                incoming = incoming.exclude(source_term__ontology_service=OntologyService.HPO)

            all_relations = list(outgoing) + list(incoming)

            for relation in all_relations:
                snake = by_leafs.get(relation.source_term) or by_leafs.get(relation.dest_term)

                if snake.leaf_term in (relation.source_term, relation.dest_term):
                    other_term = relation.other_end(snake.leaf_term)

                    ontology_services = {snake.leaf_term.ontology_service, other_term.ontology_service}
                    # Possibly Narrow or Broad would also be valid??
                    if OntologyService.MONDO in ontology_services and OntologyService.OMIM in ontology_services and relation.relation != OntologyRelation.EXACT:
                        continue

                    new_snake = snake.snake_step(relation)
                    if other_term.ontology_service == to_ontology:
                        valid_snakes.append(new_snake)
                        continue
                    if len(new_snake.paths) <= max_depth:
                        new_snakes.append(new_snake)
                    seen.add(other_term)

        valid_snakes.sort()

        return OntologySnakes(valid_snakes)

    @staticmethod
    def direct_relationships_for_gene_symbol(
            gene_symbol: Union[str, GeneSymbol],
            quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER
        ) -> list[OntologyTermRelation]:
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()
        q_relation = quality_filter.filter_q
        return list(otr_qs.filter(
            q_relation,
            dest_term=gene_ontology,
            source_term__ontology_service__in={OntologyService.MONDO, OntologyService.OMIM}
        ))

    @staticmethod
    def mondo_terms_for_gene_symbol(gene_symbol: Union[str, GeneSymbol]) -> set[OntologyTerm]:
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        terms = set()

        otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()
        q_relation = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER.filter_q

        mondos = otr_qs.filter(
            q_relation,
            dest_term=gene_ontology,
            source_term__ontology_service=OntologyService.MONDO
        )
        terms = terms.union(set(mondos.values_list("source_term_id", flat=True)))
        omim_ids = otr_qs.filter(
            q_relation,
            dest_term=gene_ontology,
            source_term__ontology_service=OntologyService.OMIM
        ).values_list("source_term_id", flat=True)

        if omim_ids:
            # relationships are always MONDO -> OMIM, and MONDO -> HGNC, OMIM -> HGNC
            via_omim_mondos = OntologyVersion.get_latest_and_live_ontology_qs().filter(
                source_term__ontology_service=OntologyService.MONDO, dest_term_id__in=omim_ids).\
                exclude(relation__in={OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM}).\
                values_list("source_term_id", flat=True)

            terms = terms.union(set(via_omim_mondos))
        if terms:
            return set(OntologyTerm.objects.filter(pk__in=terms))
        return set()

    @staticmethod
    def terms_for_gene_symbol(gene_symbol: Union[str, GeneSymbol], desired_ontology: OntologyService,
                              max_depth=1,
                              quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
                              otr_qs: QuerySet[OntologyTermRelation] = None) -> 'OntologySnakes':
        # FIXME, can the min_classification default to STRONG and other code can filter it out?
        """ max_depth: How many steps in snake path to go through """
        # TODO, do this with hooks
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        return OntologySnake.snake_from(term=gene_ontology, to_ontology=desired_ontology,
                                        max_depth=max_depth, quality_filter=quality_filter,
                                        otr_qs=otr_qs)

    @staticmethod
    def has_gene_relationship(term: Union[OntologyTerm, str], gene_symbol: Union[GeneSymbol, str], quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER) -> bool:
        # TODO, have this run off get_all_term_to_gene_relationships
        # just need to filter through results until we reach one of a high enough quality
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        if isinstance(term, str):
            term = OntologyTerm.get_or_stub(term)
            if term.is_stub:
                return False
        try:
            gene_term = OntologyTerm.get_gene_symbol(gene_symbol)
            # try direct link first
            quality_q = quality_filter.filter_q
            otr_qs = OntologyVersion.get_latest_and_live_ontology_qs().filter(dest_term=gene_term).filter(quality_q)
            if otr_qs.filter(source_term=term).exists():
                return True

            if term.ontology_service == OntologyService.OMIM:
                if mondo_term := OntologyTermRelation.as_mondo(term):
                    if otr_qs.filter(source_term=mondo_term).exists():
                        return True
            elif term.ontology_service == OntologyService.MONDO:
                if omim_term := OntologyTermRelation.as_omim(term):
                    if otr_qs.filter(source_term=omim_term).exists():
                        return True

            return False
        except ValueError:
            report_exc_info()
            return False

    @staticmethod
    def get_all_term_to_gene_relationships(term: Union[OntologyTerm, str], gene_symbol: Union[GeneSymbol, str], try_related_terms: bool = True) -> Iterator['OntologySnake']:
        # iterates all ontology term relationships between the term and the gene symbol (as well as any relationships to the equiv MONDO/OMIM)
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        if isinstance(term, str):
            term = OntologyTerm.get_or_stub(term)
            if term.is_stub:
                return None
        try:
            gene_term = OntologyTerm.get_gene_symbol(gene_symbol)
            # try direct link first
            otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()
            for relationship in otr_qs.filter(source_term=term, dest_term=gene_term):
                yield OntologySnake(source_term=term, leaf_term=gene_term, paths=[relationship])

            if not try_related_terms:
                return None

            # optimisations for OMIM/MONDO
            if term.ontology_service in {OntologyService.MONDO, OntologyService.OMIM}:
                if term.ontology_service == OntologyService.MONDO:
                    if omim := OntologyTermRelation.as_omim(term):
                        yield from OntologySnake.get_all_term_to_gene_relationships(omim, gene_symbol, try_related_terms=False)
                elif term.ontology_service == OntologyService.OMIM:
                    if mondo := OntologyTermRelation.as_mondo(term):
                        yield from OntologySnake.get_all_term_to_gene_relationships(mondo, gene_symbol, try_related_terms=False)
        except ValueError:
            report_exc_info()
            return None

    @staticmethod
    def get_children(term: OntologyTerm) -> set[OntologyTerm]:
        relationships = OntologyVersion.get_latest_and_live_ontology_qs().filter(dest_term=term, relation=OntologyRelation.IS_A).select_related("source_term")
        return set(relationship.source_term for relationship in relationships)

    @staticmethod
    def get_parents(term: OntologyTerm) -> set[OntologyTerm]:
        relationships = OntologyVersion.get_latest_and_live_ontology_qs().filter(source_term=term, relation=OntologyRelation.IS_A).select_related("dest_term")
        return set(relationship.dest_term for relationship in relationships)


class OntologySnakes:

    def __bool__(self):
        return bool(self.snakes)

    def __init__(self, snakes: list[OntologySnake]):
        self.snakes = snakes

    def __iter__(self):
        return self.snakes.__iter__()

    def __len__(self):
        return len(self.snakes)

    def __getitem__(self, item):
        return self.snakes[item]

    def leafs(self) -> list[OntologyTerm]:
        return list(sorted({snake.leaf_term for snake in self}))

    def leaf_relations(self, ontology_relation: str = None) -> list[OntologyTermRelation]:
        relations = {snake.leaf_relationship for snake in self}
        if ontology_relation:
            relations = {otr for otr in relations if otr.relation == ontology_relation}
        return list(sorted(relations))


class SingleTermH:

    def __init__(self, term: OntologyTerm):
        self.term = term
        self.children: set[SingleTermH] = set()
        self.depth = 99
        self.resolved = False

    @cached_property
    def all_relevant_children(self):
        all_relevant_children = {self.term}
        for child in self.children:
            all_relevant_children.update(child.all_relevant_children)
        return all_relevant_children

    def __repr__(self):
        return f"{self.term} depth: {self.depth}"


class AncestorCalculator:

    def __init__(self, terms: Iterable[OntologyTerm]):
        terms = set(terms)
        if not terms:
            raise ValueError("Terms must have at least one entry")
        else:
            ontology_services = set()
            for term in terms:
                ontology_services.add(term.ontology_service)
            if len(ontology_services) > 1:
                raise ValueError(f"Can't find ancestors between different ontologies {ontology_services}")
            ontology_service = first(ontology_services)
            if ontology_service not in {OntologyService.MONDO, OntologyService.HPO}:
                raise ValueError(f"Can only find common ancestor of MONDO & HPO, not {ontology_service}")

        self.base_terms = terms
        self.term_hs: dict[OntologyTerm, SingleTermH] = {}
        for term in terms:
            self.get_term(term)
        self.processed = False

    def get_term(self, term: OntologyTerm) -> SingleTermH:
        if term_h := self.term_hs.get(term):
            return term_h
        else:
            term_h = SingleTermH(term=term)
            self.term_hs[term] = term_h
            return term_h

    def get_unresolved(self):
        return [sh for sh in self.term_hs.values() if not sh.resolved]

    def set_depth(self, single_h: SingleTermH, depth: int):
        if single_h.depth > depth:
            single_h.depth = depth

            for child in single_h.children:
                self.set_depth(child, depth+1)

    def process(self) -> OntologyTerm:
        if self.processed:
            raise ValueError("Can only call process once")
        self.processed = True

        if len(self.base_terms) == 1:
            return first(self.base_terms)

        root_term: Optional[SingleTermH] = None

        unprocessed = self.get_unresolved()
        while unprocessed:
            for sh in unprocessed:
                parents = OntologySnake.get_parents(sh.term)
                if parents:
                    for parent in parents:
                        self.get_term(parent).children.add(sh)
                else:
                    root_term = sh
                sh.resolved = True

            unprocessed = self.get_unresolved()

        if not root_term:
            raise ValueError("Did not find root term going up IS_A relationships")

        self.set_depth(root_term, 0)

        best_contender = root_term
        for sh in self.term_hs.values():
            if sh.depth > best_contender.depth and sh.all_relevant_children.issuperset(self.base_terms):
                best_contender = sh

        return best_contender.term

    @staticmethod
    def common_ancestor(terms: Iterable[OntologyTerm]) -> OntologyTerm:
        return AncestorCalculator(terms).process()
