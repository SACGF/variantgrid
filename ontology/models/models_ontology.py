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
from typing import Optional, List, Dict, Set, Union, Tuple, Iterable

from cache_memoize import cache_memoize
from django.contrib.postgres.fields import ArrayField
from django.db import models, connection
from django.db.models import PROTECT, CASCADE, QuerySet, Q, Max, TextChoices
from django.urls import reverse
from lazy import lazy
from model_utils.models import TimeStampedModel, now
from psqlextra.models import PostgresPartitionedModel
from psqlextra.types import PostgresPartitioningMethod

from genes.models import GeneSymbol
from library.cache import timed_cache
from library.constants import DAY_SECS, WEEK_SECS
from library.log_utils import report_exc_info
from library.utils import Constant


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

    DOID = "DOID", "DOID"
    ORPHANET = "Orphanet", "Orphanet"

    EXPECTED_LENGTHS: Dict[str, int] = Constant({
        MONDO[0]: 7,
        OMIM[0]: 6,
        HPO[0]: 7,
        HGNC[0]: 1,  # HGNC ids aren't typically 0 padded, because they're not monsters
        DOID[0]: None,  # variable length with padded 0s
        ORPHANET[0]: 1  # ORPHANET ids aren't typically 0 padded
    })

    IMPORTANCE: Dict[str, int] = Constant({
        MONDO[0]: 2,
        OMIM[0]: 3,
        HPO[0]: 4,  # put HPO relationships last as they occasionally spam OMIM
        DOID[0]: 5,
        ORPHANET[0]: 6,
        HGNC[0]: 1  # show gene relationships first
    })

    URLS: Dict[str, str] = Constant({
        MONDO[0]: "https://monarchinitiative.org/disease/MONDO:${1}",
        OMIM[0]: "http://www.omim.org/entry/${1}",
        HPO[0]: "https://hpo.jax.org/app/browse/term/HP:${1}",
        HGNC[0]: "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:${1}",
        DOID[0]: "https://www.ebi.ac.uk/ols/ontologies/doid/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FDOID_${1}",
        ORPHANET[0]: "https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=EN&Expert=${1}"
    })

    LOCAL_ONTOLOGY_PREFIXES: Set[str] = Constant({
        MONDO[0],
        OMIM[0],
        HPO[0],
        HGNC[0]
    })

    CONDITION_ONTOLOGIES: Set[str] = Constant({
        MONDO[0],
        OMIM[0],
        HPO[0],
        DOID[0],
        ORPHANET[0]
    })

    @staticmethod
    def index_to_id(ontology_service: 'OntologyService', index: int):
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
    RELATED = "related"  # defined by HPO and MONDO (also use relatedSynonymn from mondo to populate this)
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


class GeneDiseaseClassification(models.TextChoices):
    # @see https://thegencc.org/faq.html#validity-termsdelphi-survey - where sort order comes from
    # Not using eg "GENCC:100009" (= Supportive) as that has different sort order than page above
    REFUTED = "1", "Refuted Evidence"
    NO_KNOWN = "2", "No Known Disease Relationship"
    ANIMAL = "3", "Animal Model Only"
    DISPUTED = "4", "Disputed Evidence"
    LIMITED = "5", "Limited"
    SUPPORTIVE = "6", "Supportive"
    MODERATE = "7", "Moderate"
    STRONG = "8", "Strong"
    DEFINITIVE = "9", "Definitive"

    @staticmethod
    def get_above_min(min_classification: str) -> List[str]:
        classifications = []
        for e in reversed(GeneDiseaseClassification):
            classifications.append(e.label)
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

    def save(self, **kwargs):
        created = not self.pk
        if created and not self.version:
            # Assign version to be next highest
            existing_ontology = OntologyImport.objects.filter(import_source=self.import_source,
                                                              filename=self.filename)
            data = existing_ontology.aggregate(Max("version"))
            self.version = (data.get("version__max") or 0) + 1

        super().save(**kwargs)
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

    @staticmethod
    def normalize(dirty_id: str) -> 'OntologyIdNormalized':
        parts = re.split("[:|_]", dirty_id)
        if len(parts) != 2:
            raise ValueError(f"Can not convert {dirty_id} to a proper id")

        prefix = parts[0].strip().upper()
        if prefix == "ORPHANET":  # Orphanet is the one ontology (so far) where the standard is sentance case
            prefix = "Orphanet"
        prefix = OntologyService(prefix)
        postfix = parts[1].strip()
        try:
            num_part = int(postfix)
            clean_id: str
            if expected_length := OntologyService.EXPECTED_LENGTHS[prefix]:
                clean_id = OntologyService.index_to_id(prefix, num_part)
            else:
                # variable length IDs like DOID
                clean_id = f"{prefix}:{postfix}"

            return OntologyIdNormalized(prefix=prefix, postfix=postfix, full_id=clean_id, clean=True)

        except ValueError:
            return OntologyIdNormalized(prefix=prefix, postfix=postfix, full_id=dirty_id, clean=False)

    def __str__(self):
        return self.full_id


class OntologyTerm(TimeStampedModel):

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

    def __str__(self):
        return f"{self.id} {self.name}"

    class Meta:
        unique_together = ("ontology_service", "index")

    def __lt__(self, other):
        if self.ontology_service != other.ontology_service:
            return self.ontology_service < other.ontology_service
        return self.index < other.index

    def get_absolute_url(self):
        return reverse('ontology_term', kwargs={"term": self.url_safe_id})

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
            term_type = (self.extra or dict()).get('type', 'Unknown')
            return f"Term is of type - {term_type}"
        elif self.status == OntologyTermStatus.STUB:
            return f"Term was referenced by 3rd party but not yet from our authoritative source"
        else:
            return None

    @lazy
    def is_leaf(self) -> bool:
        # Warning, just meant to be called on MONDO terms
        if not self.is_stub and self.ontology_service == OntologyService.MONDO:
            return not OntologyTermRelation.objects.filter(dest_term=self, relation=OntologyRelation.IS_A).exists()
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
        WONT persist it to the database
        """
        if not isinstance(id_str, OntologyIdNormalized):
            normal_id = OntologyIdNormalized.normalize(id_str)
        if normal_id.clean:
            if existing := OntologyTerm.objects.filter(id=normal_id.full_id).first():
                return existing
            return OntologyTerm(
                id=normal_id.full_id,
                ontology_service=normal_id.prefix,
                index=normal_id.num_part,
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
    def url(self):
        return OntologyService.URLS[self.ontology_service].replace("${1}", self.padded_index)

    @staticmethod
    def split_hpo_omim_mondo(ontology_term_ids: Iterable[str]) -> Tuple[QuerySet, QuerySet, QuerySet]:
        hpo_qs = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.HPO)
        omim_qs = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.OMIM)
        mondo_qs = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.MONDO)
        return hpo_qs, omim_qs, mondo_qs

    @staticmethod
    def split_hpo_omim_mondo_as_dict(ontology_term_ids: Iterable[str]) -> Dict[str, QuerySet]:
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

    I haven't elected to use django_dag node_factory here as it only allows one kind of relationship
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

    @staticmethod
    def as_mondo(term: OntologyTerm) -> Optional[OntologyTerm]:
        if term.ontology_service == OntologyService.MONDO:
            return term

        q_dest_modo = Q(source_term=term) & Q(dest_term__ontology_service=OntologyService.MONDO)
        q_source_modo = Q(dest_term=term) & Q(source_term__ontology_service=OntologyService.MONDO)
        otr_qs = OntologyTermRelation.objects.filter(q_dest_modo | q_source_modo, relation=OntologyRelation.EXACT)
        if mondo_rel := otr_qs.first():
            return mondo_rel.other_end(term)
        return None

    @staticmethod
    def relations_of(term: OntologyTerm, otr_qs: Optional[QuerySet['OntologyTermRelation']] = None) -> List['OntologyTermRelation']:
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
    def get_moi_summary(moi_classifications, valid_classifications) -> List[str]:
        moi_summary = []
        for moi, classifications in moi_classifications.items():
            classification_submitters = []
            for classification in valid_classifications:
                if submitters := classifications.get(classification):
                    classification_submitters.append(f"{classification}: {'/'.join(sorted(submitters))}")
            if classification_submitters:
                moi_summary.append(f"{moi} ({' '.join(classification_submitters)})")
        return moi_summary


OntologyList = Optional[Union[QuerySet, List[OntologyTerm]]]


class OntologyVersion(TimeStampedModel):
    """ This is used by annotation.AnnotationVersion to keep track of different OntologyImports
        so we can load historical versions of OntologyTermRelation """

    gencc_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="gencc_ontology_version")
    mondo_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="mondo_ontology_version")
    hp_owl_import = models.ForeignKey(OntologyImport, on_delete=PROTECT, related_name="hp_owl_ontology_version")
    hp_phenotype_to_genes_import = models.ForeignKey(OntologyImport, on_delete=PROTECT,
                                                     related_name="hp_phenotype_to_genes_ontology_version")

    class Meta:
        unique_together = ('gencc_import', 'mondo_import', 'hp_owl_import', 'hp_phenotype_to_genes_import')

    ONTOLOGY_IMPORTS = {
        "gencc_import": (OntologyImportSource.GENCC,
                         ['https://search.thegencc.org/download/action/submissions-export-csv',
                          'gencc-submissions.csv']),
        "mondo_import": (OntologyImportSource.MONDO, ['mondo.json']),
        "hp_owl_import": (OntologyImportSource.HPO, ['hp.owl']),
        "hp_phenotype_to_genes_import": (OntologyImportSource.HPO,
                                         ['phenotype_to_genes.txt',
                                          'OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt']),
    }

    @staticmethod
    def in_ontology_version(ontology_import: OntologyImport) -> bool:
        versioned = defaultdict(set)
        for (import_source, filenames) in OntologyVersion.ONTOLOGY_IMPORTS.values():
            versioned[import_source].update(filenames)
        return ontology_import.filename in versioned[ontology_import.import_source]

    @staticmethod
    def latest() -> Optional['OntologyVersion']:
        oi_qs = OntologyImport.objects.all()
        kwargs = {}
        for field, (import_source, filenames) in OntologyVersion.ONTOLOGY_IMPORTS.items():
            kwargs[field] = oi_qs.filter(import_source=import_source, filename__in=filenames).order_by("pk").last()

        values = list(kwargs.values())
        if all(values):
            last_date = max([oi.created for oi in values])
            ontology_version, created = OntologyVersion.objects.get_or_create(**kwargs,
                                                                              defaults={"created": last_date})
            if created:
                # Avoid circular import
                from annotation.models import AnnotationVersion
                AnnotationVersion.new_sub_version(None)
        else:
            ontology_version = None
            missing_fields = [field for field, value in kwargs.items() if value is None]
            if missing_fields:
                msg = "OntologyVersion.latest() - missing fields: %s", ", ".join(missing_fields)
                raise OntologyVersion.DoesNotExist(msg)
        return ontology_version

    def get_ontology_imports(self):
        return [self.gencc_import, self.mondo_import, self.hp_owl_import, self.hp_phenotype_to_genes_import]

    def get_ontology_term_relations(self):
        return OntologyTermRelation.objects.filter(from_import__in=self.get_ontology_imports())

    @staticmethod
    def get_latest_and_live_ontology_qs():
        latest = OntologyVersion.latest()
        # live relationships of panelappau aren't versioned
        # TODO could restrict only if we have live enabled in settings
        return OntologyTermRelation.objects.filter(Q(from_import__in=latest.get_ontology_imports()) | Q(relation='panelappau'))

    def get_gene_disease_relations_qs(self) -> QuerySet:
        return self.get_ontology_term_relations().filter(relation=OntologyRelation.RELATED,
                                                         extra__strongest_classification__isnull=False)

    @cache_memoize(WEEK_SECS)
    def moi_and_submitters(self) -> Tuple[List[str], List[str]]:
        """ Cached lists of MOI/Submitters from GenCC gene/disease extra JSON """
        moi = set()
        submitters = set()
        for extra in self.get_gene_disease_relations_qs().values_list("extra", flat=True):
            for source in extra["sources"]:
                moi.add(source["mode_of_inheritance"])
                submitters.add(source["submitter"])
        return list(sorted(moi)), list(sorted(submitters))

    @cache_memoize(DAY_SECS)
    def cached_gene_symbols_for_terms_tuple(self, terms_tuple: Tuple[int]) -> QuerySet:
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
                              max_depth=1, min_classification: Optional[GeneDiseaseClassification] = None) -> 'OntologySnakes':
        otr_qs = self.get_ontology_term_relations()
        return OntologySnake.terms_for_gene_symbol(gene_symbol, desired_ontology, max_depth=max_depth,
                                                   min_classification=min_classification, otr_qs=otr_qs)

    def gene_disease_relations(self, gene_symbol: Union[str, GeneSymbol],
                               min_classification: GeneDiseaseClassification = None) -> List[OntologyTermRelation]:
        snake = self.terms_for_gene_symbol(gene_symbol, OntologyService.MONDO,
                                           max_depth=0, min_classification=min_classification)
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

    @property
    def relationship(self) -> OntologyTermRelation:
        # relationship is the preferred term, add this property to help migrate over to better wording
        return self.relation

    @property
    def source_term(self) -> OntologyTerm:
        return self.relation.other_end(self.dest_term)


class OntologySnake:
    """
    Use to "Snake" through Ontology nodes, typically to resolve to/from gene symbols.
    An OntologySnake is meant to be immutable, creating new snake_steps each time snake_step is called
    """

    def __init__(self, source_term: OntologyTerm, leaf_term: Optional[OntologyTerm] = None,
                 paths: Optional[List[OntologyTermRelation]] = None):
        self.source_term = source_term
        self.leaf_term = leaf_term or source_term
        self.paths = paths or []

    def snake_step(self, relationship: OntologyTermRelation) -> 'OntologySnake':
        """
        Creates a new OntologySnake with this extra relationship
        """
        new_leaf = relationship.other_end(self.leaf_term)
        new_paths = list(self.paths)
        new_paths.append(relationship)
        return OntologySnake(source_term=self.source_term, leaf_term=new_leaf, paths=new_paths)

    def show_steps(self) -> List[OntologySnakeStep]:
        steps: List[OntologySnakeStep] = []
        node = self.source_term
        for path in self.paths:
            node = path.other_end(node)
            steps.append(OntologySnakeStep(relation=path, dest_term=node))
        return steps

    def __str__(self):
        text = f"{self.source_term}"
        for step in self.show_steps():
            forwards = step.relation.dest_term == step.dest_term
            text += f" {'<' if not forwards else ''}-{step.relation.relation}-{'>' if forwards else ''} {step.dest_term}"
        return text

    @lazy
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

    @staticmethod
    def check_if_ancestor(descendant: OntologyTerm, ancestor: OntologyTerm, max_levels=4) -> List['OntologySnake']:
        if ancestor == descendant:
            return OntologySnake(source_term=ancestor, leaf_term=descendant)

        if descendant.ontology_service != ancestor.ontology_service:
            raise ValueError(f"Can only check for ancestry within the same ontology service, not {descendant.ontology_service} vs {ancestor.ontology_service}")

        seen: Set[OntologyTerm] = {descendant}
        new_snakes: List[OntologySnake] = list([OntologySnake(source_term=descendant)])
        valid_snakes: List[OntologySnake] = []
        level = 0
        while new_snakes:
            level += 1
            snakes_by_leaf: Dict[OntologyTerm, OntologySnake] = {}
            for snake in new_snakes:
                snakes_by_leaf[snake.leaf_term] = snake

            new_snakes: List[OntologySnake] = []
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

    # TODO only allow EXACT between two anythings that aren't Gene Symbols
    @staticmethod
    def snake_from(term: OntologyTerm, to_ontology: OntologyService,
                   min_classification: Optional[GeneDiseaseClassification] = GeneDiseaseClassification.STRONG,
                   max_depth: int = 1, otr_qs: QuerySet[OntologyTermRelation] = None) -> 'OntologySnakes':
        """
        Returns the smallest snake/paths from source term to the desired OntologyService
        Ignores IS_A paths
        """
        if term.ontology_service == to_ontology:
            return OntologySnakes([OntologySnake(source_term=term)])

        if otr_qs is None:
            otr_qs = OntologyVersion.get_latest_and_live_ontology_qs()
            # otr_qs = OntologyTermRelation.objects.all()

        seen: Set[OntologyTerm] = set()
        seen.add(term)
        new_snakes: List[OntologySnake] = list([OntologySnake(source_term=term)])
        valid_snakes: List[OntologySnake] = []

        relation_q_list = [
            ~Q(relation__in={OntologyRelation.IS_A, OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM}),
            OntologySnake.gencc_quality_filter(min_classification),
        ]
        q_relation = functools.reduce(operator.and_, relation_q_list)

        iteration = -1
        while new_snakes:
            iteration += 1
            snakes: List[OntologySnake] = list(new_snakes)
            new_snakes: List[OntologySnake] = []
            by_leafs: Dict[OntologyTerm, OntologySnake] = {}
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

                if relation.source_term == snake.leaf_term or relation.dest_term == snake.leaf_term:
                    other_term = relation.other_end(snake.leaf_term)

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
    def gencc_quality_filter(quality: GeneDiseaseClassification = GeneDiseaseClassification.STRONG) -> Q:
        gencc_classifications = GeneDiseaseClassification.get_above_min(GeneDiseaseClassification.STRONG)
        return ~Q(from_import__import_source='gencc') | Q(extra__strongest_classification__in=gencc_classifications)

    @staticmethod
    def mondo_terms_for_gene_symbol(gene_symbol: Union[str, GeneSymbol]) -> Set[OntologyTerm]:
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        terms = set()

        q_relation = OntologySnake.gencc_quality_filter()

        mondos = OntologyTermRelation.objects.filter(q_relation, dest_term=gene_ontology,
                                                     source_term__ontology_service=OntologyService.MONDO)
        terms = terms.union(set(mondos.values_list("source_term_id", flat=True)))
        omim_ids = OntologyTermRelation.objects.filter(q_relation, dest_term=gene_ontology,
                                                       source_term__ontology_service=OntologyService.OMIM).values_list("source_term_id", flat=True)
        if omim_ids:
            # relationships are always MONDO -> OMIM, and MONDO -> HGNC, OMIM -> HGNC
            via_omim_mondos = OntologyTermRelation.objects.filter(source_term__ontology_service=OntologyService.MONDO, dest_term_id__in=omim_ids).exclude(relation__in={OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM}).values_list("source_term_id", flat=True)
            terms = terms.union(set(via_omim_mondos))
        if terms:
            return set(OntologyTerm.objects.filter(pk__in=terms))
        return set()

    @staticmethod
    def terms_for_gene_symbol(gene_symbol: Union[str, GeneSymbol], desired_ontology: OntologyService,
                              max_depth=1, min_classification: Optional[GeneDiseaseClassification] = None,
                              otr_qs: QuerySet[OntologyTermRelation] = None) -> 'OntologySnakes':
        # FIXME, can the min_classification default to STRONG and other code can filter it out?
        """ max_depth: How many steps in snake path to go through """
        # TODO, do this with hooks
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        return OntologySnake.snake_from(term=gene_ontology, to_ontology=desired_ontology,
                                        max_depth=max_depth, min_classification=min_classification,
                                        otr_qs=otr_qs)

    @staticmethod
    def has_gene_relationship(term: Union[OntologyTerm, str], gene_symbol: Union[GeneSymbol, str], quality: Optional[GeneDiseaseClassification] = GeneDiseaseClassification.STRONG) -> bool:
        # TODO, do this with hooks
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        if isinstance(term, str):
            term = OntologyTerm.get_or_stub(term)
            if term.is_stub:
                return False
        try:
            gene_term = OntologyTerm.get_gene_symbol(gene_symbol)
            # try direct link first
            quality_q = OntologySnake.gencc_quality_filter(quality)
            if OntologyTermRelation.objects.filter(source_term=term, dest_term=gene_term).filter(quality_q).exists():
                return True
            # optimisations for OMIM/MONDO
            if term.ontology_service in {OntologyService.MONDO, OntologyService.OMIM}:
                via_ids: QuerySet = None
                exclude_mondo_omim = ~Q(relation__in={OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM})
                if term.ontology_service == OntologyService.MONDO:
                    via_ids = OntologyTermRelation.objects.filter(source_term=term, dest_term__ontology_service=OntologyService.OMIM).filter(quality_q).filter(exclude_mondo_omim).values_list("dest_term_id", flat=True)
                else:
                    via_ids = OntologyTermRelation.objects.filter(dest_term=term, source_term__ontology_service=OntologyService.MONDO).filter(quality_q).filter(exclude_mondo_omim).values_list("source_term_id", flat=True)
                return OntologyTermRelation.objects.filter(source_term_id__in=via_ids, dest_term=gene_term).exists()

            hgnc_terms = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC).leafs()
            return gene_term in hgnc_terms
        except ValueError:
            report_exc_info()
            return False


class OntologySnakes:

    def __bool__(self):
        return bool(self.snakes)

    def __init__(self, snakes: List[OntologySnake]):
        self.snakes = snakes

    def __iter__(self):
        return self.snakes.__iter__()

    def __len__(self):
        return len(self.snakes)

    def __getitem__(self, item):
        return self.snakes[item]

    def leafs(self) -> List[OntologyTerm]:
        return list(sorted({snake.leaf_term for snake in self}))

    def leaf_relations(self, ontology_relation: str = None) -> List[OntologyTermRelation]:
        relations = {snake.leaf_relationship for snake in self}
        if ontology_relation:
            relations = {otr for otr in relations if otr.relation == ontology_relation}
        return list(sorted(relations))
