import operator
import re
from functools import reduce
from typing import Set, Iterable, Union, Optional, Match, List, Any

from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import QuerySet
from django.db.models.query_utils import Q
from django.urls.base import reverse
from lazy import lazy
from pyhgvs import InvalidHGVSName

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.manual_variant_entry import check_can_create_variants, CreateManualVariantForbidden
from annotation.models.models import AnnotationVersion
from genes.hgvs import HGVSMatcher
from genes.models import TranscriptVersion, Transcript, MissingTranscript, Gene, GeneSymbol, GeneSymbolAlias, \
    gene_symbol_withdrawn_str
from genes.models_enums import AnnotationConsortium
from library.genomics import format_chrom
from library.log_utils import report_exc_info
from ontology.models import OntologyTerm
from patients.models import ExternalPK, Patient
from seqauto.models import SequencingRun, Experiment
from snpdb.clingen_allele import get_clingen_allele, get_clingen_alleles_from_external_code
from snpdb.models import VARIANT_PATTERN, LOCUS_PATTERN, LOCUS_NO_REF_PATTERN, DBSNP_PATTERN, \
    ClinGenAllele, GenomeBuild, Sample, Variant, Sequence, VariantCoordinate, UserSettings, Organization, Lab, Allele
from snpdb.models.models_enums import ClinGenAlleleExternalRecordType
from upload.models import ModifiedImportedVariant
from classification.models.classification import ClassificationModification, \
    Classification, CreateNoClassificationForbidden
from variantgrid.perm_path import get_visible_url_names
from variantopedia.models import SearchTypes

DB_PREFIX_PATTERN = re.compile(fr"^(v|{settings.VARIANT_VCF_DB_PREFIX})(\d+)$")
VARIANT_VCF_PATTERN = re.compile(r"((?:chr)?\S*)\s+(\d+)\s+\.?\s*([GATC]+)\s+([GATC]+)")
VARIANT_GNOMAD_PATTERN = re.compile(r"(?:chr)?(\S*)-(\d+)-([GATC]+)-([GATC]+)")
HGVS_MINIMUM_TO_SHOW_ERROR_PATTERN = re.compile(r":(c|g|p)\..*\d+")
ONTOLOGY_PATTERN = re.compile(r"\w+:\s*.*")

class AbstractMetaVariant:

    def is_valid_for_user(self, user: User) -> bool:
        return True


class CreateManualVariant(AbstractMetaVariant):
    def __init__(self, genome_build: GenomeBuild, variant_string: str):
        self.genome_build = genome_build
        self.variant_string = variant_string

    def is_valid_for_user(self, user: User) -> bool:
        try:
            check_can_create_variants(user)
            return True
        except CreateManualVariantForbidden:
            return False

    def get_absolute_url(self):
        kwargs = {"genome_build_name": self.genome_build.pk, "variants_text": self.variant_string}
        return reverse('create_manual_variant_entry_from_text', kwargs=kwargs)

    def __str__(self):
        return "Click to create and annotate this variant"


class ClassifyVariant(AbstractMetaVariant):

    def __init__(self, variant: Variant, transcript_id: str = None):
        self.variant = variant
        self.transcript_id = transcript_id

    def get_absolute_url(self):
        kwargs = {"variant_id": self.variant.pk}
        if self.transcript_id:
            url = "create_classification_for_variant_and_transcript"
            kwargs["transcript_id"] = self.transcript_id
        else:
            url = "create_classification_for_variant"
        return reverse(url, kwargs=kwargs)

    def __str__(self):
        return f"Click to classify '{self.variant}'"


class ClassifyNoVariantHGVS(AbstractMetaVariant):

    def __init__(self, genome_build: GenomeBuild, hgvs_string: str):
        self.genome_build = genome_build
        self.hgvs_string = hgvs_string

    def is_valid_for_user(self, user: User) -> bool:
        try:
            Classification.check_can_create_no_classification_via_web_form(user)
            return True
        except CreateNoClassificationForbidden:
            return False

    def get_absolute_url(self):
        kwargs = {"genome_build_name": self.genome_build.pk, "hgvs_string": self.hgvs_string}
        return reverse('create_classification_from_hgvs', kwargs=kwargs)

    def __str__(self):
        return f"Click to classify from unvalidated HGVS: '{self.hgvs_string}'"


def search_data(user: User, search_string: str, classify: bool) -> 'SearchResults':
    """ build_results:  list of (search_type, result, genome_build)
                        result - has string representation to display
                                 get_absolute_url() link or to load if only result
        search_types: list of what we searched due to regex match
        search_errors: list of (search_type, exception, genome_build)
    """
    searcher = Searcher(user=user, search_string=search_string, classify=classify)
    return searcher.search()


class SearchResult:

    def __init__(self, record: Any,
                 genome_builds: Optional[Set[GenomeBuild]] = None,
                 annotation_consortia: Optional[str] = None,
                 message: Optional[Union[str, List[str]]] = None,
                 is_debug_data: bool = False):
        self.record = record
        self.genome_builds = genome_builds
        self.annotation_consortia = annotation_consortia
        search_message = getattr(record, "search_message", None)
        if search_message and not message:
            message = search_message
        if isinstance(message, str):
            message = [message]

        self.messages: Optional[List[str]] = message
        self.is_debug_data = is_debug_data
        self.search_type: Optional[str] = None

        self.is_preferred_build: Optional[bool] = None
        self.is_preferred_annotation: Optional[bool] = None

    def apply_search_score(self, preferred_genome_build: GenomeBuild):
        self.is_preferred_build = not self.genome_builds or preferred_genome_build in self.genome_builds
        self.is_preferred_annotation = not self.annotation_consortia or preferred_genome_build.annotation_consortium == self.annotation_consortia

    @property
    def header_string(self):
        text = self.search_type
        extras = []
        if self.annotation_consortia:
            if self.annotation_consortia == AnnotationConsortium.ENSEMBL:
                extras.append("Ensembl")
            elif self.annotation_consortia == AnnotationConsortium.REFSEQ:
                extras.append("RefSeq")
        if self.genome_builds:
            extras.extend((genome_build.name for genome_build in self.genome_builds))

        if extras:
            text = f"{text} ({', '.join(extras)})"

        return text

    @property
    def is_preferred(self):
        return self.is_preferred_build and self.is_preferred_annotation and not self.is_debug_data and not self.messages

    @lazy
    def preference_score(self):
        if self.is_debug_data:
            return 0

        score = 0
        if self.is_preferred_build:
            score += 4
        if self.is_preferred_annotation:
            score += 2
        if not self.messages:
            score += 1
        return score

    def __lt__(self, other):
        # sort by score first
        score_diff = self.preference_score - other.preference_score
        if score_diff:
            return score_diff < 0

        # then inverse by type (e.g. so big scores and alphabetically small go up the top)
        if other.search_type < self.search_type:
            return True
        if other.search_type > self.search_type:
            return False

        if str(other.record) < str(self.record):
            return True
        return False


class SearchResults:

    def __init__(self, genome_build_preferred: GenomeBuild):
        self.genome_build_preferred = genome_build_preferred
        self.results: List[SearchResult] = []
        self.search_types: Set[str] = set()
        self.search_errors: List = []

    def append_search_type(self, search_type: str):
        self.search_types.add(search_type)

    def append_result(self, result: Any):
        self.results.append(result)

    def append_error(self, error: Any):
        self.search_errors.append(error)

    def extend(self, other: 'SearchResults'):
        self.search_errors.extend(other.search_errors)
        self.results.extend(other.results)
        self.search_types = self.search_types.union(other.search_types)

    @property
    def search_types_string(self):
        return ', '.join(self.search_types)

    def complete(self):
        for result in self.results:
            result.apply_search_score(preferred_genome_build=self.genome_build_preferred)
        self.results.sort(reverse=True)

    @property
    def non_debug_results(self):
        return [result for result in self.results if not result.is_debug_data]

    def single_preferred_result(self) -> Optional[SearchResult]:
        preferred = [result for result in self.results if result.is_preferred]
        if len(preferred) == 1:
            return preferred[0]
        return None


class Searcher:

    def __init__(self, user: User, search_string: str, classify: bool):
        """ build_results:  list of (search_type, result, genome_build)
                            result - has string representation to display
                                     get_absolute_url() link or to load if only result
            search_types: list of what we searched due to regex match
            search_errors: list of (search_type, exception, genome_build)
        """
        HAS_ALPHA_PATTERN = r"[a-zA-Z]"
        NO_WHITESPACE = r"^\S+$"
        HGVS_UNCLEANED_PATTERN = r"[cmg].*\d+"  # Bare bones match - may be able to fix it...
        COSMIC_PATTERN = re.compile(r"^COS[MV](\d+)$")  # Old or new

        self.genome_build_searches = [
            (SearchTypes.DBSNP, DBSNP_PATTERN, search_dbsnp),
            (SearchTypes.HGVS, HGVS_UNCLEANED_PATTERN, search_hgvs),
            (SearchTypes.LOCUS, LOCUS_PATTERN, search_locus),
            (SearchTypes.LOCUS, LOCUS_NO_REF_PATTERN, search_locus),
            (SearchTypes.SAMPLE, HAS_ALPHA_PATTERN, search_sample),
            (SearchTypes.VARIANT, VARIANT_PATTERN, search_variant),
            (SearchTypes.CLINGEN_ALLELE_ID, r"^CA", search_clingen_allele),
            (SearchTypes.VARIANT, VARIANT_VCF_PATTERN, search_variant_vcf),
            (SearchTypes.VARIANT, VARIANT_GNOMAD_PATTERN, search_variant_gnomad)
        ]
        self.genome_agnostic_searches = [
            (SearchTypes.GENE, NO_WHITESPACE, search_gene_symbol),  # special case
            (SearchTypes.GENE, r"(ENSG|Gene.*:)\d+", search_gene),  # special case
            (SearchTypes.COSMIC, COSMIC_PATTERN, search_cosmic),
            (SearchTypes.EXPERIMENT, HAS_ALPHA_PATTERN, search_experiment),
            (SearchTypes.EXTERNAL_PK, HAS_ALPHA_PATTERN, search_external_pk),
            (SearchTypes.CLASSIFICATION, None, search_classification),
            (SearchTypes.PATIENT, HAS_ALPHA_PATTERN, search_patient),
            (SearchTypes.SEQUENCING_RUN, r"^\d{6}", search_sequencing_run),
            (SearchTypes.TRANSCRIPT, r"^(ENST|NM_)\d+\.?\d*$", search_transcript),
            (SearchTypes.VARIANT, DB_PREFIX_PATTERN, search_variant_id),
            (SearchTypes.LAB, r"[a-zA-Z]{3,}", search_lab),
            (SearchTypes.ORG, r"[a-zA-Z]{3,}", search_org),
            (SearchTypes.ONTOLOGY, ONTOLOGY_PATTERN, search_ontology)
        ]

        exclude_search_types = set()
        visible_url_parts = get_visible_url_names()
        if not visible_url_parts.get('patients'):
            exclude_search_types.add(SearchTypes.PATIENT)
            exclude_search_types.add(SearchTypes.SAMPLE)
            exclude_search_types.add(SearchTypes.EXTERNAL_PK)
        if not visible_url_parts.get('data'):
            exclude_search_types.add(SearchTypes.EXPERIMENT)
            exclude_search_types.add(SearchTypes.SEQUENCING_RUN)

        self.exclude_search_types = exclude_search_types

        # if user can only see variants with classifications, then they def shouldn't be able to classify new variants
        self.can_create = True
        self.user = user
        self.search_string = search_string.strip()
        self.classify = classify

        if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not user.is_superuser:
            self.classify = False
            self.can_create = False

        user_settings = UserSettings.get_for_user(user)
        self.genome_build_preferred = user_settings.default_genome_build
        self.genome_build_all = GenomeBuild.builds_with_annotation()

    def search(self) -> SearchResults:
        search_results = SearchResults(genome_build_preferred=self.genome_build_preferred)

        # do genome build agnostic search
        search_results.extend(self.search_within(genome_build=None))
        for genome_build in self.genome_build_all:
            # do genome build specific searches
            search_results.extend(self.search_within(genome_build=genome_build))

        search_results.complete()
        return search_results

    def search_within(self, genome_build: Optional[GenomeBuild] = None) -> SearchResults:
        """
        :param genome_build: If provided, only searches that use genome_build will be run, if None will search genome_build agnostic searches
        :param annotation_consortia: If provided, only annotation consortia searches will be run
        :return: SearchResults
        """
        print(f"Searching within genome_build {genome_build}")

        variant_qs: Optional[QuerySet]
        searches: List[Any]
        results = SearchResults(genome_build_preferred=genome_build)
        if genome_build:
            try:
                variant_qs = get_visible_variants(self.user, genome_build)
                searches = self.genome_build_searches
            except Exception as e:
                report_exc_info(extra_data={
                    'search_string': self.search_string,
                    'classify': self.classify,
                    'genome_build_id': genome_build.name if genome_build else None
                })
                results.append_error(("variants", e, genome_build))
                return results
        else:
            variant_qs = None
            searches = self.genome_agnostic_searches

        genome_build_results = set()
        if genome_build:
            genome_build_results.add(genome_build)

        for search_type, regex_pattern, get_results_func in searches:
            if search_type in self.exclude_search_types:
                continue

            if regex_pattern is None or re.search(regex_pattern, self.search_string):
                results.append_search_type(search_type)
                try:
                    build_results = get_results_func(self.search_string,
                                                     user=self.user,
                                                     genome_build=genome_build,
                                                     variant_qs=variant_qs,
                                                     classify=self.classify)
                    if build_results is not None:
                        for sr in build_results:
                            if not isinstance(sr, SearchResult):
                                sr = SearchResult(record=sr)
                            r = sr.record

                            if isinstance(r, Variant) and self.classify:
                                sr.record = ClassifyVariant(r)

                            if isinstance(r, AbstractMetaVariant) and not r.is_valid_for_user(user=self.user):
                                # don't let regular users create new variants if they can only see variants with classifications attached
                                # don't count as filtered out result, as it's not a "real" result
                                continue

                            if not sr.search_type:
                                sr.search_type = search_type
                            if not sr.genome_builds:
                                sr.genome_builds = genome_build_results

                            results.append_result(sr)

                except Exception as e:
                    report_exc_info(extra_data={
                        'search_string': self.search_string,
                        'search_type': search_type,
                        'classify': self.classify,
                        'genome_build_id': genome_build.name if genome_build else None
                    })
                    results.append_error((search_type, e, genome_build))
        return results


VARIANT_SEARCH_RESULTS = Optional[Union[Iterable[Union[ClassifyVariant, Allele, Variant, CreateManualVariant, SearchResult]], QuerySet]]


def get_visible_variants(user: User, genome_build: GenomeBuild) -> VARIANT_SEARCH_RESULTS:
    """ Shariant wants to restrict search to only classified variants """

    annotation_version = AnnotationVersion.latest(genome_build)
    variant_qs = get_variant_queryset_for_annotation_version(annotation_version)
    if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not user.is_superuser:
        variant_qs = variant_qs.filter(Classification.get_variant_q(user, genome_build))

    return variant_qs


def search_cosmic(search_string: str, user: User, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    results = None
    if search_string.startswith("COSV"):
        results = variant_qs.filter(variantannotation__cosmic_id__contains=search_string)
    elif search_string.startswith("COSM"):
        results = variant_qs.filter(variantannotation__cosmic_legacy_id__contains=search_string)
    return results


def search_dbsnp(search_string, user, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    results = variant_qs.filter(variantannotation__dbsnp_rs_id=search_string)
    if not results.exists():
        cmv_list = []
        for clingen_allele in get_clingen_alleles_from_external_code(ClinGenAlleleExternalRecordType.DBSNP_ID, search_string):
            variant_string = clingen_allele.get_variant_string(genome_build)
            search_message = f"'{search_string}' resolved to ClinGenAllele {clingen_allele} coordinate '{variant_string}'"
            cmv_list.append(SearchResult(CreateManualVariant(genome_build, variant_string), message=search_message))
        return cmv_list

    return results


def search_gene_symbol(search_string: str, **kwargs) -> Iterable[Union[GeneSymbol, GeneSymbolAlias]]:
    # itertools.chain doesn't work
    # only return a GeneSymbol alias if we're not returning the source GeneSymbol
    gene_symbols = list(GeneSymbol.objects.filter(symbol__iexact=search_string).exclude(symbol__endswith=gene_symbol_withdrawn_str))
    gene_symbol_strs = {gene_symbol.symbol for gene_symbol in gene_symbols}
    aliases = [alias for alias in GeneSymbolAlias.objects.filter(alias__iexact=search_string).all() if alias.alias not in gene_symbol_strs]
    return gene_symbols + aliases


def search_ontology(search_string: str, **kwargs) -> Optional[SearchResult]:
    try:
        return [SearchResult(OntologyTerm.get_or_stub(search_string))]
    except:
        return []


def search_gene(search_string: str, **kwargs) -> Iterable[Gene]:
    """ Symbols have been separated into search_gene_symbol - this returns Gene objects """

    CONSORTIUM_REGEX = {
        r"(ENSG\d+)": AnnotationConsortium.ENSEMBL,
        r"Gene:(\d+)": AnnotationConsortium.REFSEQ,
        r"GeneID:(\d+)": AnnotationConsortium.REFSEQ,
        r"Gene ID:(\d+)": AnnotationConsortium.REFSEQ,
    }

    for c_regex, annotation_consortium in CONSORTIUM_REGEX.items():
        if m := re.match(c_regex, search_string, re.IGNORECASE):
            gene_id = m.group(1)
            return Gene.objects.filter(identifier=gene_id, annotation_consortium=annotation_consortium)
    return []


def search_transcript(search_string: str, **kwargs) -> VARIANT_SEARCH_RESULTS:
    """ return Transcript or TranscriptVersion (build independent) """
    transcript_id, version = TranscriptVersion.get_transcript_id_and_version(search_string)
    if version:
        obj = TranscriptVersion.objects.filter(transcript_id=transcript_id, version=version).first()
    else:
        obj = Transcript.objects.filter(identifier=transcript_id).first()

    if obj:
        return [obj]
    return []


def search_hgvs(search_string: str, user: User, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    hgvs_matcher = HGVSMatcher(genome_build)
    classify = kwargs.get("classify")
    try:
        # can add on search_message to objects to (stop auto-jump and show message)
        search_messages = []
        hgvs_string = search_string
        try:
            variant_tuple = hgvs_matcher.get_variant_tuple(hgvs_string)
        except (InvalidHGVSName, NotImplementedError) as original_error:
            original_hgvs_string = hgvs_string
            try:
                hgvs_string = HGVSMatcher.clean_hgvs(hgvs_string)
                variant_tuple = hgvs_matcher.get_variant_tuple(hgvs_string)
                search_messages.append(f"Warning: Cleaned '{search_string}' => '{hgvs_string}'")
            except InvalidHGVSName:
                if classify:
                    search_message = f"Error reading HGVS: '{original_error}'"
                    return [SearchResult(ClassifyNoVariantHGVS(genome_build, original_hgvs_string), message=search_message)]

                # We want to be able to rescue dodgy HGVS but don't want to raise errors
                # for everything as that will cause lots of false positives.
                if HGVS_MINIMUM_TO_SHOW_ERROR_PATTERN.findall(search_string):
                    raise original_error  # cleaning didn't work don't tell anyone
                return None

        if not hgvs_matcher.matches_reference(hgvs_string):
            ref = variant_tuple[2]
            build_and_patch = genome_build.get_build_with_patch()
            if len(ref) > 20:
                ref = f"{Sequence.abbreviate(ref)} ({len(ref)} bases)"
            search_messages.append(f"Warning: Using reference '{ref}' from our build: {build_and_patch}")

        original_tv, used_tv = hgvs_matcher.get_original_and_used_transcript_versions(hgvs_string)
        if original_tv != used_tv:
            search_messages.append(f"Warning: Missing transcript version '{original_tv}', using best match '{used_tv}'")

        transcript_id = hgvs_matcher.get_transcript_id(hgvs_string,
                                                       transcript_version=False)
        try:
            results = get_results_from_variant_tuples(variant_qs, variant_tuple)
            variant = results.get()
            if classify:
                return [SearchResult(ClassifyVariant(variant, transcript_id), message=search_messages)]
            return [SearchResult(variant, message=search_messages)]
        except Variant.DoesNotExist:
            variant_string = Variant.format_tuple(*variant_tuple)
            variant_string_abbreviated = Variant.format_tuple(*variant_tuple, abbreviate=True)
            search_messages.append(f"'{search_string}' resolved to {variant_string_abbreviated}")
            results = [SearchResult(CreateManualVariant(genome_build, variant_string), message=search_messages)]
            results.extend(search_for_alt_alts(variant_qs, variant_tuple, search_messages))
            return results
    except MissingTranscript:
        pass


def search_for_alt_alts(variant_qs: QuerySet, variant_tuple: VariantCoordinate, messages: List[str]) -> VARIANT_SEARCH_RESULTS:
    if not messages:
        messages = list()
    variants = list()
    results = get_results_from_variant_tuples(variant_qs, variant_tuple, any_alt=True)
    for result in results:
        this_messages = list()
        this_messages.extend(messages)
        this_messages.append(f"Warning: No results for alt '{variant_tuple.alt}', but found this using alt '{result.alt}'")
        variants.append(SearchResult(record=result, message=this_messages))
    return variants


def search_locus(search_string: str, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    if m := re.match(LOCUS_PATTERN, search_string):
        chrom, position, ref = m.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
    else:
        if m := re.match(r"([^:]+):(\d+)", search_string):  # No ref supplied
            chrom, position = m.groups()
            chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
            ref = None
        else:
            return None

    kwargs = {"locus__contig__name": chrom,
              "locus__position": position}
    if ref:
        kwargs["locus__ref__seq"] = ref
    return variant_qs.filter(**kwargs)


def search_sample(search_string: str, user: User, genome_build: GenomeBuild, **kwargs) -> Iterable[Sample]:
    return Sample.filter_for_user(user).filter(name__icontains=search_string,
                                               vcf__genome_build=genome_build)


def search_variant_match(m: Match, user: User, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    if m:
        (chrom, position, ref, alt) = m.groups()
        ref = ref.upper()
        alt = alt.upper()
        if ref == alt:
            alt = Variant.REFERENCE_ALT
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        results = get_results_from_variant_tuples(variant_qs, (chrom, position, ref, alt))
        if results.exists():
            return results

        variant_string = Variant.format_tuple(chrom, position, ref, alt)
        search_message = f"The variant {variant_string} does not exist in the database"
        return [SearchResult(CreateManualVariant(genome_build, variant_string), message=search_message)]
    return None


def search_variant(search_string: str, user: User, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    m = VARIANT_PATTERN.match(search_string)
    return search_variant_match(m, user, genome_build, variant_qs, **kwargs)


def search_variant_vcf(search_string: str, user: User, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    m = VARIANT_VCF_PATTERN.match(search_string)
    return search_variant_match(m, user, genome_build, variant_qs, **kwargs)


def search_variant_gnomad(search_string: str, user: User, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    m = VARIANT_GNOMAD_PATTERN.match(search_string)
    return search_variant_match(m, user, genome_build, variant_qs, **kwargs)


def search_variant_id(search_string: str, **kwargs) -> VARIANT_SEARCH_RESULTS:
    m = re.match(DB_PREFIX_PATTERN, search_string)
    if m:
        variant_id = m.group(2)
        return Variant.objects.filter(pk=variant_id)
    return None


def search_experiment(search_string: str, **kwargs) -> Iterable[Experiment]:
    return Experiment.objects.filter(name__icontains=search_string)


def search_external_pk(search_string: str, user: User, **kwargs):
    # Returns related objects
    # TODO: Handle user permissions on this??
    RELATED_OBJECT_FIELDS = ["case", "pathologytestorder", "patient"]
    results = []
    for external_pk in ExternalPK.objects.filter(code__iexact=search_string):
        for f in RELATED_OBJECT_FIELDS:
            try:
                obj = getattr(external_pk, f)
                results.append(obj)
            except ObjectDoesNotExist:
                pass
    return results


def search_classification(search_string: str, user: User, **kwargs) -> Iterable[Classification]:
    """ Search for LabId which can be either:
        "vc1080" or "Molecular Genetics, Frome Road / vc1080" (as it appears in classification) """

    filters = [Q(classification__lab_record_id=search_string)]  # exact match
    slash_index = search_string.find("/")
    if slash_index > 0:
        lab_name = search_string[:slash_index].strip()
        lab_record_id = search_string[slash_index + 1:].strip()
        if lab_name and lab_record_id:
            q_lab_name = Q(classification__lab__name=lab_name)
            q_lab_record = Q(classification__lab_record_id=lab_record_id)
            filters.append(q_lab_name & q_lab_record)
    q_vcm = reduce(operator.or_, filters)
    vcm_qs = ClassificationModification.filter_for_user(user)
    vcm_ids = vcm_qs.filter(q_vcm, is_last_published=True).values('classification')
    # convert from modifications back to Classification so absolute_url returns the editable link
    return Classification.objects.filter(pk__in=vcm_ids)


def search_clingen_allele(search_string: str, user: User, genome_build: GenomeBuild, variant_qs: QuerySet, **kwargs) -> VARIANT_SEARCH_RESULTS:
    if ClinGenAllele.looks_like_id(search_string):
        clingen_allele = get_clingen_allele(search_string)
        if settings.PREFER_ALLELE_LINKS:
            return [clingen_allele.allele]

        variant_qs = variant_qs.filter(variantallele__allele__clingen_allele=clingen_allele,
                                       variantallele__genome_build=genome_build)
        if variant_qs.exists():
            return variant_qs
        variant_string = clingen_allele.get_variant_string(genome_build)
        variant_string_abbreviated = clingen_allele.get_variant_string(genome_build, abbreviate=True)
        search_message = f"'{clingen_allele}' resolved to '{variant_string_abbreviated}'"
        return [SearchResult(CreateManualVariant(genome_build, variant_string), message=search_message)]

    return None


def search_patient(search_string: str, user: User, **kwargs) -> Iterable[Patient]:
    parts = search_string.split(",")
    if len(parts) == 2:
        (last_name, first_name) = parts
        q_last = Q(last_name__iexact=last_name.strip())
        q_first = Q(first_name__iexact=first_name.strip())
        patient_q = q_last & q_first
    else:
        patient_q = Q(last_name__iexact=search_string) | Q(first_name__iexact=search_string)
    return Patient.filter_for_user(user).filter(patient_q)


def search_sequencing_run(search_string: str, **kwargs) -> Iterable[SequencingRun]:
    return SequencingRun.objects.filter(name__icontains=search_string)


def get_results_from_variant_tuples(qs: QuerySet, data: VariantCoordinate, any_alt: bool = False) -> Union[QuerySet, Iterable[Variant]]:
    """
    :param qs: A query set that we'll be searching inside of (except for when returning ModifiedImportVariants)
    :param data: The variant coordinate to lookup
    :param any_alt: If true, search without using alt and return all matches
    :return: A QuerySet of variants
    """
    (chrom, position, ref, alt) = data
    position = int(position)

    results = qs.filter(Variant.get_chrom_q(chrom), locus__position=position, locus__ref__seq=ref)
    if not any_alt:
        results = results.filter(alt__seq=alt)

    if not results:
        if not any_alt:
            results = ModifiedImportedVariant.get_variants_for_unnormalized_variant(chrom, position, ref, alt)
        else:
            # should we really be searching ModifiedImportVariants with any alt? or should that just happen for
            # the filter of "real" variants
            results = ModifiedImportedVariant.get_variants_for_unnormalized_variant_any_alt(chrom, position, ref)
    return results


def search_lab(search_string: str, **kwargs) -> Iterable[Lab]:
    return Lab.objects.filter(organization__active=True).filter(name__icontains=search_string)


def search_org(search_string: str, **kwargs) -> Iterable[Organization]:
    return Organization.objects.filter(active=True).filter(Q(short_name__icontains=search_string) | Q(name__icontains=search_string))
