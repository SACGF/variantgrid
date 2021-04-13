import functools
import re
from dataclasses import dataclass
from typing import Optional, List, Dict, Set, Union, Tuple, Iterable

from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models import PROTECT, CASCADE, QuerySet, Q
from django.urls import reverse
from lazy import lazy
from model_utils.models import TimeStampedModel, now

from genes.models import GeneSymbol
from library.log_utils import report_exc_info
from library.utils import Constant

"""
A series of models that currently stores the combination of MONDO, OMIM, HPO & HGNC.
(Note that HGNC is included just to model relationships, please use GeneSymbol for all your GeneSymbol needs).
"""
class OntologyImportSource:
    PANEL_APP_AU = "PAAU"
    MONDO = "MONDO"
    OMIM = "OMIM"
    HPO = "HP"
    HGNC = "HGNC"


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
        MONDO[0]: "https://vm-monitor.monarchinitiative.org/disease/MONDO:${1}",
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


class OntologyImport(TimeStampedModel):
    """
    Keeps track of when data was imported, typically used to see how old the data is and if it needs
    to be imported again
    """
    import_source = models.TextField()
    filename = models.TextField()
    context = models.TextField()
    hash = models.TextField()
    processor_version = models.IntegerField(default=1)
    processed_date = models.DateTimeField(auto_created=True)
    completed = models.BooleanField(default=False)


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
        if self.name:
            return "obsolete" in self.name.lower()
        return False

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
    def get_or_stub(id_str: str) -> 'OntologyTerm':
        """
        Returns an OntologyTerm for the given ID.
        If the OntologyTerm doesn't exist in the database, will create an OntologyTerm but
        WONT persist it to the database
        """
        parts = re.split("[:|_]", id_str)
        if len(parts) != 2:
            raise ValueError(f"Can not convert {id_str} to a proper id")

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

            if existing := OntologyTerm.objects.filter(id=clean_id).first():
                return existing

            return OntologyTerm(
                id=clean_id,
                ontology_service=prefix,
                index=num_part,
                name=""
            )
        except ValueError:
            if existing := OntologyTerm.objects.filter(ontology_service=prefix, name=postfix).first():
                return existing
            raise ValueError(f"Can not convert {id_str} to a proper id")

    @property
    def padded_index(self) -> str:
        # ID should already be padded to right number of indexes
        return self.id.split(":")[1]

    @property
    def url(self):
        return OntologyService.URLS[self.ontology_service].replace("${1}", self.padded_index)

    @staticmethod
    def split_hpo_and_omim(ontology_term_ids: Iterable[str]) -> Tuple[QuerySet, QuerySet]:
        hpo = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.HPO)
        omim = OntologyTerm.objects.filter(pk__in=ontology_term_ids, ontology_service=OntologyService.OMIM)
        return hpo, omim


class OntologyTermRelation(TimeStampedModel):
    """
    Relationship between two terms, is generally considered to be bi-directional (or at the very least
    code typically checks relationships in both directions)

    I haven't elected to use django_dag node_factory here as it only allows one kind of relationship
    and we have quite a lot.
    """
    source_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE, related_name="subject")
    dest_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    relation = models.TextField()
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)

    class Meta:
        unique_together = ("source_term", "dest_term", "relation")

    # as in we can have multiple copies, up to code to try to not duplicate

    def __str__(self):
        return f"{self.source_term} -> ({self.relation}) -> {self.dest_term}"

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
    def parents_of(source_term: OntologyTerm) -> QuerySet:
        """
        QuerySet of OntologyTerms representing parents of this term.
        Note that MONDO terms very often have multiple parents
        """
        parent_ids = OntologyTermRelation.objects.filter(source_term=source_term, relation=OntologyRelation.IS_A).values_list("dest_term", flat=True)
        return OntologyTerm.objects.filter(id__in=parent_ids).order_by("-id")

    @staticmethod
    def children_of(term: OntologyTerm) -> QuerySet:
        """
        QuerySet of OntologyTerms representing children of this term
        """
        children_ids = OntologyTermRelation.objects.filter(dest_term=term, relation=OntologyRelation.IS_A).values_list("source_term", flat=True)
        return OntologyTerm.objects.filter(id__in=children_ids).order_by("-id")

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
    def relations_of(term: OntologyTerm) -> List['OntologyTermRelation']:
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

        items = list(OntologyTermRelation.objects.filter(Q(source_term=term) | Q(dest_term=term)).select_related("source_term", "dest_term", "from_import"))
        items.sort(key=functools.cmp_to_key(sort_relationships))
        return items


@dataclass
class OntologySnakeStep:
    """
    A step in the snake, used for showing the path
    """
    relation: OntologyTermRelation
    dest_term: OntologyTerm


OntologyList = Optional[Union[QuerySet, List[OntologyTerm]]]


class OntologySnake:
    """
    Use to "Snake" through Ontology nodes, typically to resolve to/from gene symbols.
    An OntologySnake is meant to be immutable, creating new snake_steps each time snake_step is called
    """

    def __init__(self, source_term: OntologyTerm, leaf_term: Optional[OntologyTerm] = None,
                 paths: Optional[List[OntologyTermRelation]] = None):
        self.source_term = source_term
        self.leaf_term = leaf_term or source_term
        self.paths = paths or list()

    def snake_step(self, relationship: OntologyTermRelation) -> 'OntologySnake':
        """
        Creates a new OntologySnake with this extra relationship
        """
        new_leaf = relationship.other_end(self.leaf_term)
        new_paths = list(self.paths)
        new_paths.append(relationship)
        return OntologySnake(source_term=self.source_term, leaf_term=new_leaf, paths=new_paths)

    def show_steps(self) -> List[OntologySnakeStep]:
        steps: List[OntologySnakeStep] = list()
        node = self.source_term
        for path in self.paths:
            node = path.other_end(node)
            steps.append(OntologySnakeStep(relation=path, dest_term=node))
        return steps

    @property
    def leaf_relationship(self) -> OntologyTermRelation:
        return self.paths[-1]

    # TODO only allow EXACT between two anythings that aren't Gene Symbols
    @staticmethod
    def snake_from(term: OntologyTerm, to_ontology: OntologyService, max_depth: int = 1) -> 'OntologySnakes':
        """
        Returns the smallest snake/paths from source term to the desired OntologyService
        Ignores IS_A paths
        """
        if term.ontology_service == to_ontology:
            return OntologySnakes([OntologySnake(source_term=term)])

        seen: Set[OntologyTerm] = set()
        seen.add(term)
        new_snakes: List[OntologySnake] = list([OntologySnake(source_term=term)])
        valid_snakes: List[OntologySnake] = list()

        best_relationships = ~Q(relation__in={OntologyRelation.IS_A, OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM})

        iteration = -1
        while new_snakes:
            iteration += 1
            snakes: List[OntologySnake] = list(new_snakes)
            new_snakes: List[OntologySnake] = list()
            by_leafs: Dict[OntologyTerm, OntologySnake] = dict()
            for snake in snakes:
                if existing := by_leafs.get(snake.leaf_term):
                    if len(snake.paths) < len(existing.paths):
                        by_leafs[snake.leaf_term] = snake
                else:
                    by_leafs[snake.leaf_term] = snake
            all_leafs = by_leafs.keys()

            outgoing = OntologyTermRelation.objects.filter(source_term__in=all_leafs).exclude(dest_term__in=seen).filter(best_relationships).select_related("source_term", "dest_term")
            incoming = OntologyTermRelation.objects.filter(dest_term__in=all_leafs).exclude(source_term__in=seen).filter(best_relationships).select_related("source_term", "dest_term")
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

        return OntologySnakes(valid_snakes)

    @staticmethod
    def mondo_terms_for_gene_symbol(gene_symbol: Union[str, GeneSymbol]) -> Set[OntologyTerm]:
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        terms = set()

        mondos = OntologyTermRelation.objects.filter(dest_term=gene_ontology, source_term__ontology_service=OntologyService.MONDO).values_list("source_term_id", flat=True)
        terms = terms.union(set(mondos))
        omim_ids = OntologyTermRelation.objects.filter(dest_term=gene_ontology, source_term__ontology_service=OntologyService.OMIM).values_list("source_term_id", flat=True)
        if omim_ids:
            # relationships are always MONDO -> OMIM, and MONDO -> HGNC, OMIM -> HGNC
            via_omim_mondos = OntologyTermRelation.objects.filter(source_term__ontology_service=OntologyService.MONDO, dest_term_id__in=omim_ids).exclude(relation__in={OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM}).values_list("source_term_id", flat=True)
            terms = terms.union(set(via_omim_mondos))
        if terms:
            return set(OntologyTerm.objects.filter(pk__in=terms))
        return set()

    @staticmethod
    def terms_for_gene_symbol(gene_symbol: Union[str, GeneSymbol], desired_ontology: OntologyService, max_depth=1) -> 'OntologySnakes':
        # TODO, do this with hooks
        from ontology.panel_app_ontology import update_gene_relations
        update_gene_relations(gene_symbol)
        gene_ontology = OntologyTerm.get_gene_symbol(gene_symbol)
        return OntologySnake.snake_from(term=gene_ontology, to_ontology=desired_ontology, max_depth=max_depth)

    @staticmethod
    def gene_symbols_for_terms(terms: OntologyList) -> QuerySet:
        gene_symbol_names = set()
        for term in terms:
            if isinstance(term, str):
                term = OntologyTerm.get_or_stub(term)
                if term.is_stub:
                    return GeneSymbol.objects.none()
            gene_symbol_snakes = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC)
            gene_symbol_names.update([snake.leaf_term.name for snake in gene_symbol_snakes])
        return GeneSymbol.objects.filter(symbol__in=gene_symbol_names)

    @staticmethod
    def has_gene_relationship(term: Union[OntologyTerm, str], gene_symbol: Union[GeneSymbol, str]) -> bool:
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
            if OntologyTermRelation.objects.filter(source_term=term, dest_term=gene_term).exists():
                return True
            # optimisations for OMIM/MONDO
            if term.ontology_service in {OntologyService.MONDO, OntologyService.OMIM}:
                via_ids: QuerySet = None
                exclude_mondo_omim = ~Q(relation__in={OntologyRelation.EXACT_SYNONYM, OntologyRelation.RELATED_SYNONYM})
                if term.ontology_service == OntologyService.MONDO:
                    via_ids = OntologyTermRelation.objects.filter(source_term=term, dest_term__ontology_service=OntologyService.OMIM).filter(exclude_mondo_omim).values_list("dest_term_id", flat=True)
                else:
                    via_ids = OntologyTermRelation.objects.filter(dest_term=term, source_term__ontology_service=OntologyService.MONDO).filter(exclude_mondo_omim).values_list("source_term_id", flat=True)
                return OntologyTermRelation.objects.filter(source_term_id__in=via_ids, dest_term=gene_term).exists()

            hgnc_terms = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC).leafs()
            return gene_term in hgnc_terms
        except ValueError:
            report_exc_info()
            return False


class OntologySnakes:

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
