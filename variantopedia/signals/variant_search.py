import logging
import re
from collections import defaultdict
from collections.abc import Callable, Iterable
from typing import Optional, Union

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.urls import reverse

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.cosmic import CosmicAPI
from annotation.manual_variant_entry import CreateManualVariantForbidden, check_can_create_variants
from annotation.models import AnnotationVersion
from classification.models import Classification, CreateNoClassificationForbidden
from genes.gene_copy_number import (
    COPY_NUMBER_STRING_PATTERN,
    find_gene_copy_number_events_for_string,
)
from genes.gene_fusions import FUSION_STRING_SEPARATOR, find_gene_fusions_for_string
from genes.gene_splice import SPLICE_STRING_PATTERN, find_splice_events_for_string
from genes.hgvs import (
    HGVSException,
    HGVSMatcher,
)
from library.enums.log_level import LogLevel
from library.genomics import format_chrom
from library.preview_request import PreviewData
from snpdb.clingen_allele import get_clingen_allele
from snpdb.models import (
    DBSNP_PATTERN,
    LOCUS_NO_REF_PATTERN,
    LOCUS_PATTERN,
    VARIANT_PATTERN,
    VARIANT_SYMBOLIC_PATTERN,
    Allele,
    ClinGenAllele,
    DbSNP,
    GenomeBuild,
    Sequence,
    Variant,
    VariantCoordinate,
)
from snpdb.search import (
    SearchExample,
    SearchInputInstance,
    SearchMessage,
    SearchMessageOverall,
    SearchResult,
    search_receiver,
)
from snpdb.variant_coordinate_search import get_results_from_variant_coordinate

COSMIC_PATTERN = re.compile(r"^(COS[VM])[0-9]{3,}$", re.IGNORECASE)


def get_visible_variants(search_input: SearchInputInstance, genome_build: GenomeBuild) -> QuerySet[Variant]:
    """ The build's variants, restricted on Shariant to those with a classification published to a non-admin user """
    # Search receivers filter Variants by locus/position/pk/alt and don't join to annotation
    # partitions, so a plain queryset returns identical results without the partition transformer.
    if annotation_version := AnnotationVersion.latest_or_none(genome_build, context="Search"):
        variant_qs = get_variant_queryset_for_annotation_version(annotation_version)
    else:
        variant_qs = Variant.objects.all()
    variant_qs = variant_qs.filter(Variant.get_contigs_q(genome_build))  # restrict to build
    if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not search_input.user.is_superuser:
        variant_qs = variant_qs.filter(Classification.objects.published_variant_q(search_input.user, genome_build))

    return variant_qs


class VariantExtra:

    @staticmethod
    def create_manual_variant(for_user: User, genome_build: GenomeBuild, variant_string: str) -> Optional['PreviewData']:
        try:
            check_can_create_variants(for_user)
        except CreateManualVariantForbidden:
            return None

        kwargs = {"genome_build_name": genome_build.pk, "variants_text": variant_string}
        internal_url = reverse('create_manual_variant_entry_from_text', kwargs=kwargs)
        return PreviewData(
            category="Variant",
            identifier=variant_string,
            icon="fa-solid fa-circle-plus",
            title="Click to create and annotate this variant",
            internal_url=internal_url,
            genome_builds={genome_build},
            is_operation=True
        )

    @staticmethod
    def classify_variant(variant: Variant, genome_build) -> Optional[PreviewData]:
        kwargs = {"variant_id": variant.pk, "genome_build_name": genome_build.name}
        name = "create_classification_for_variant"
        internal_url = reverse(name, kwargs=kwargs)
        parts = [str(variant)]
        return PreviewData(
            category="Variant",
            identifier=" ".join(parts),
            icon="fa-solid fa-circle-plus",
            title="Click to classify variant",
            internal_url=internal_url,
            genome_builds=variant.genome_builds,
            obj=variant,
            is_operation=True
        )

    @staticmethod
    def classify_no_variant_hgvs(for_user: User, genome_build: GenomeBuild, hgvs_string: str) -> Optional[PreviewData]:
        try:
            Classification.check_can_create_no_classification_via_web_form(for_user)
        except CreateNoClassificationForbidden:
            return None
        kwargs = {"genome_build_name": genome_build.pk, "hgvs_string": hgvs_string}
        internal_url = reverse('create_classification_from_hgvs', kwargs=kwargs)
        return PreviewData(
            category="Variant",
            identifier=hgvs_string,
            icon="fa-solid fa-circle-plus",
            title=f"Click to classify from unvalidated HGVS: '{hgvs_string}'",
            internal_url=internal_url,
            genome_builds={genome_build},
            is_operation=True
        )


@search_receiver(
    search_type=Variant,
    pattern=COSMIC_PATTERN,
    sub_name="COSMIC",
    example=SearchExample(
        note="Provide the COSMIC ID",
        examples=["COSV53567516", "COSM7659094"]
    ),
    enabled=settings.SEARCH_COSMIC_ENABLED,
)
def variant_cosmic_search(search_input: SearchInputInstance):
    # Do via API as a full table scan takes way too long with big data
    unsupported_builds = []
    supported_builds = []
    for genome_build in search_input.genome_builds:
        if CosmicAPI.supports_genome_build(genome_build):
            supported_builds.append(genome_build)
        else:
            unsupported_builds.append(genome_build)

    if unsupported_builds and not supported_builds:
        build_names = ", ".join(sorted(str(gb) for gb in unsupported_builds))
        msg = f"COSMIC search is only available for GRCh37/GRCh38, not {build_names}"
        yield SearchMessageOverall(msg, severity=LogLevel.WARNING, genome_builds=unsupported_builds)

    for genome_build in supported_builds:
        matcher = HGVSMatcher.instance(genome_build)
        cosmic = CosmicAPI(search_input.search_string, genome_build)
        results_by_variant_identifier: dict[str, list[SearchResult]] = defaultdict(list)
        hgvs_by_variant_identifier: dict[str, list[str]] = defaultdict(list)

        for hgvs_string in cosmic.get_hgvs_list():
            for result in _yield_results_from_hgvs(search_input, genome_build, matcher, hgvs_string):
                if result.preview.category == "Variant":
                    variant_identifier = result.preview.identifier
                    results_by_variant_identifier[variant_identifier].append(result)
                    hgvs_by_variant_identifier[variant_identifier].append(hgvs_string)

        for variant_identifier, results_for_record in results_by_variant_identifier.items():
            if settings.SEARCH_COSMIC_TRANSCRIPT_MESSAGES:
                transcript_hgvs = ', '.join(hgvs_by_variant_identifier[variant_identifier])
                messages = [
                    SearchMessage(f"COSMIC {search_input.search_string} resolved to: {transcript_hgvs}",
                                  severity=LogLevel.INFO)
                ]
            else:
                messages = []
            preview = results_for_record[0].preview  # These will all be the same
            yield SearchResult(preview=preview, messages=messages)


@search_receiver(
    search_type=Variant,
    pattern=LOCUS_NO_REF_PATTERN,
    sub_name="Locus w/out Ref",
    example=SearchExample(
        note="Provide chromosome and position",
        examples=["chr1:169519049"]
    )
)
def search_variant_locus_no_ref(search_input: SearchInputInstance):
    for genome_build in search_input.genome_builds:
        chrom, position = search_input.match.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        yield get_visible_variants(search_input, genome_build).filter(
            locus__contig__name=chrom,
            locus__position=position
        )


@search_receiver(
    search_type=Variant,
    pattern=LOCUS_PATTERN,
    sub_name="Locus with Ref",
    example=SearchExample(
        note="Provide chromosome and position and a reference",
        examples=["chr1:169519049 A"]
    )
)
def search_variant_locus_with_ref(search_input: SearchInputInstance):
    for genome_build in search_input.genome_builds:
        chrom, position, ref = search_input.match.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        variants_qs = get_visible_variants(search_input, genome_build)
        position_qs = variants_qs.filter(
            locus__contig__name=chrom,
            locus__position=position,
        )
        ref_qs = position_qs.filter(locus__ref__seq=ref)
        if ref_qs.exists():
            yield from ref_qs
        else:
            for v in position_qs:
                ref_messages = [
                    SearchMessage(f'No results for ref "{ref}", but found this using ref "{v.locus.ref}"',
                                  severity=LogLevel.WARNING, substituted=True)
                ]
                yield SearchResult(v.preview, messages=ref_messages)


@search_receiver(
    search_type=Variant,
    pattern=ClinGenAllele.CLINGEN_ALLELE_CODE_PATTERN,
    sub_name="ClinGen Allele ID",
    example=SearchExample(note="ClinGen Allele ID", examples=["CA285410130"])
)
def allele_search(search_input: SearchInputInstance):
    search_string = search_input.search_string
    if ClinGenAllele.looks_like_id(search_string):
        clingen_allele = get_clingen_allele(search_string)
        allele = clingen_allele.allele
        for genome_build in search_input.genome_builds:
            # This ensures we retrieve using correct permissions
            visible_variants_qs = get_visible_variants(search_input, genome_build)
            if variant := allele.variant_for_build_optional(genome_build):
                yield from visible_variants_qs.filter(pk=variant.pk)
                continue
            # if here then No variant found

            # If we are non-admin we don't want to see anything else (eg ClinGen not in our build)
            if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not search_input.user.is_superuser:
                continue

            try:
                variant_coordinate = clingen_allele.get_variant_coordinate(genome_build)
                yield from yield_no_results_for_variant_coordinate(search_input.user, genome_build, visible_variants_qs,
                                                                   variant_coordinate, [])
            except ValueError as e:
                yield SearchMessageOverall(str(e), severity=LogLevel.ERROR, genome_builds=[genome_build])


def _sv_alt_description(v: Union[Variant, VariantCoordinate]) -> str:
    return f"{v.alt}/SVLEN={v.svlen}"


def _alt_description(v: Union[Variant, VariantCoordinate]) -> str:
    """Describe the alt of a Variant or VariantCoordinate. Handles both Variant.alt (a Sequence) and
       VariantCoordinate.alt (a str) - str() abbreviates a Sequence, while Sequence.abbreviate handles the str."""
    if v.is_symbolic:
        return _sv_alt_description(v)
    return Sequence.abbreviate(str(v.alt))


def _alt_mismatch_message(variant_coordinate: VariantCoordinate, alternative_variant: Variant,
                          hgvs_matcher: Optional[HGVSMatcher]) -> str:
    """ Alts alone ("AC" vs "AT") are hard to compare for indels, so an HGVS search names both as g.HGVS """
    if hgvs_matcher:
        try:
            searched_g_hgvs = hgvs_matcher.variant_coordinate_to_g_hgvs(variant_coordinate)
            found_g_hgvs = hgvs_matcher.variant_to_g_hgvs(alternative_variant)
            return f'No results for "{searched_g_hgvs}", but found "{found_g_hgvs}" at the same position'
        except HGVSException:
            pass
    original_alt_desc = _alt_description(variant_coordinate)
    alt_alt_desc = _alt_description(alternative_variant)
    return f'No results for alt "{original_alt_desc}", but found this using alt "{alt_alt_desc}"'


def yield_no_results_for_variant_coordinate(user, genome_build: GenomeBuild, variant_qs,
                                            variant_coordinate: VariantCoordinate,
                                            search_messages: list[SearchMessage],
                                            hgvs_matcher: Optional[HGVSMatcher] = None) -> Iterable[SearchResult]:
    # manual variants
    variant_string = variant_coordinate.format()
    if cmv := VariantExtra.create_manual_variant(for_user=user, genome_build=genome_build,
                                                 variant_string=variant_string):
        yield SearchResult(cmv, messages=search_messages)

    # search for alt alts
    alts = get_results_from_variant_coordinate(genome_build, variant_qs, variant_coordinate, any_alt=True)
    for alternative_variant in alts:
        message = _alt_mismatch_message(variant_coordinate, alternative_variant, hgvs_matcher)
        alt_messages = [*search_messages, SearchMessage(message, severity=LogLevel.ERROR, substituted=True)]
        yield SearchResult(alternative_variant.preview, messages=alt_messages)


def yield_search_variant_match(search_input: SearchInputInstance, get_variant_coordinate: Callable):
    for genome_build in search_input.genome_builds:
        variant_coordinate = get_variant_coordinate(search_input.match, genome_build)
        # Coordinate is explicit/symbolic according to length
        if errors := Variant.validate(genome_build, variant_coordinate.chrom, variant_coordinate.position):
            yield SearchMessageOverall(", ".join(errors), genome_builds=[genome_build])
            continue
        if message := variant_coordinate.unsupported_symbolic_alt_message:
            yield SearchMessageOverall(message, genome_builds=[genome_build])
            continue

        search_messages = []
        if not variant_coordinate.is_symbolic:
            # Only check refs if explicit, non-symbolic refs
            calculated_ref = variant_coordinate.calculated_reference(genome_build)
            logging.info("calc ref: %s, provided ref: %s", calculated_ref, variant_coordinate.ref)
            if calculated_ref != variant_coordinate.ref:
                provided_ref = variant_coordinate.ref
                sm = SearchMessage(f'Using genomic reference "{calculated_ref}" from our build, in place of provided reference "{provided_ref}"',
                                   LogLevel.WARNING, substituted=True)
                search_messages.append(sm)
                variant_coordinate.ref = calculated_ref

        visible_variants_qs = get_visible_variants(search_input, genome_build)
        results = get_results_from_variant_coordinate(genome_build, visible_variants_qs, variant_coordinate)
        if results.exists():
            search_coordinate = variant_coordinate.as_contig_accession(genome_build)
            for v in results:
                messages = search_messages
                if v.coordinate.as_contig_accession(genome_build) != search_coordinate:  # found via ModifiedImportedVariant
                    msg = f'"{variant_coordinate.format_short()}" was normalised to this variant during VCF import'
                    messages = [*search_messages, SearchMessage(msg, LogLevel.WARNING, substituted=True)]
                yield SearchResult(v.preview, messages=messages)
        else:
            yield from yield_no_results_for_variant_coordinate(search_input.user, genome_build, visible_variants_qs,
                                                               variant_coordinate, search_messages)


VARIANT_VCF_PATTERN = re.compile(r"((?:chr)?\S*)\s+(\d+)\s+\.?\s*([GATC]+)\s+([GATC]+)", re.IGNORECASE)


@search_receiver(
    search_type=Variant,
    pattern=VARIANT_VCF_PATTERN,
    sub_name="VCF Format",
    example=SearchExample(
        note="Variant coordinate as found in a VCF (across multiple columns)",
        examples=["1 169519049 T C"]
    )
)
def variant_search_vcf(search_input: SearchInputInstance):
    return yield_search_variant_match(search_input, VariantCoordinate.from_variant_match)


VARIANT_GNOMAD_PATTERN = re.compile(r"(?:chr)?(\S*)\s*-\s*(\d+)\s*-\s*([GATC]+)\s*-\s*([GATC]+)", re.IGNORECASE)


@search_receiver(
    search_type=Variant,
    pattern=VARIANT_GNOMAD_PATTERN,
    sub_name="gnomAD Format",
    example=SearchExample(
        note="A variant coordinate as formatted by gnomAD",
        examples=["1-169519049-T-C"]
    )
)
def search_variant_gnomad(search_input: SearchInputInstance):
    return yield_search_variant_match(search_input, VariantCoordinate.from_variant_match)


@search_receiver(
    search_type=Variant,
    pattern=VARIANT_PATTERN,
    sub_name="IGV Format",
    example=SearchExample(
        note="Variant coordinate as seen in IGV",
        examples=["1:169519049 T>C"]
    )
)
def search_variant_variant(search_input: SearchInputInstance):
    return yield_search_variant_match(search_input, VariantCoordinate.from_variant_match)


@search_receiver(
    search_type=Variant,
    pattern=VARIANT_SYMBOLIC_PATTERN,
    sub_name="Variant format (symbolic)",
    example=SearchExample(
        examples=["1:169519049-169520049 <DEL>", "1:169519049-169520049 <DUP>"]
    )
)
def search_variant_symbolic(search_input: SearchInputInstance):
    return yield_search_variant_match(search_input, VariantCoordinate.from_symbolic_match)


@search_receiver(
    search_type=Variant,
    pattern=DBSNP_PATTERN,
    sub_name="dbSNP ID",
    example=SearchExample(
        note="Provide the variant's dbSNP ID",
        examples=["rs6025"]
    )
)
def search_variant_db_snp(search_input: SearchInputInstance):
    # Do via API as a full table scan takes way too long with big data
    dbsnp = DbSNP.get(search_input.search_string)

    for genome_build in search_input.genome_builds:
        matcher = HGVSMatcher.instance(genome_build)
        for data in dbsnp.get_alleles_for_genome_build(genome_build):
            if hgvs_string := data.get("hgvs"):
                search_message = SearchMessage(f'dbSNP "{search_input.search_string}" resolved to "{hgvs_string}"',
                                               severity=LogLevel.INFO)
                for r in _yield_results_from_hgvs(search_input, genome_build, matcher, hgvs_string):
                    yield r, search_message


def _yield_results_from_hgvs(search_input, genome_build, matcher, hgvs_string) -> Iterable[SearchResult]:
    variant_coordinate = matcher.get_variant_coordinate(hgvs_string)
    results = get_results_from_variant_coordinate(genome_build, get_visible_variants(search_input, genome_build),
                                                  variant_coordinate)
    if results.exists():
        yield from results
    else:
        variant_string = Variant.format_tuple(*variant_coordinate)
        if create_manual := VariantExtra.create_manual_variant(
                for_user=search_input.user,
                variant_string=variant_string,
                genome_build=genome_build
        ):
            yield SearchResult(create_manual)


DB_PREFIX_PATTERN = re.compile(fr"^(v|{settings.VARIANT_VCF_DB_PREFIX})(\d+)$")


@search_receiver(
    search_type=Variant,
    pattern=DB_PREFIX_PATTERN,
    sub_name="VCF ID",
    example=SearchExample(
        note="Provide the variant ID as it appears in VCF exports",
        examples=["v42", f"{settings.VARIANT_VCF_DB_PREFIX}4638674"]
    )
)
def search_variant_id(search_input: SearchInputInstance):
    variant_id = search_input.match.group(2)
    if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not search_input.user.is_superuser:
        for genome_build in search_input.genome_builds:
            yield get_visible_variants(search_input, genome_build).filter(pk=variant_id)
    else:
        yield Variant.objects.filter(pk=variant_id)


ALLELE_ID_SEARCH_PATTERN = re.compile(r"^a(\d+)$")

@search_receiver(
    search_type=Allele,
    pattern=ALLELE_ID_SEARCH_PATTERN,
    sub_name="Internal Allele ID",
    example=SearchExample(
        note="Provide the Allele ID as it appears in URLs prefixed with \"a\"",
        examples=["a42", "a256"]
    )
)
def search_allele_id(search_input: SearchInputInstance):
    allele_qs = Allele.objects.filter(pk=search_input.match.group(1))
    if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not search_input.user.is_superuser:
        published_qs = Classification.objects.published_for_user(search_input.user)
        allele_qs = allele_qs.filter(pk__in=published_qs.values("allele"))
    yield allele_qs


# The separators classification import accepts - never a single hyphen, which gene names contain (HLA-A)
GENE_FUSION_PATTERN = re.compile(rf"^([A-Za-z0-9.\-]+?)\s*(?:{FUSION_STRING_SEPARATOR.pattern})\s*([A-Za-z0-9.\-]+)$")


@search_receiver(
    search_type=Variant,
    pattern=GENE_FUSION_PATTERN,
    sub_name="Gene Fusion",
    example=SearchExample(
        note="A gene fusion, named by its gene pair",
        examples=["BCR::ABL1", "CD74--ROS1"]
    ),
    enabled=settings.VARIANT_GENE_LEVEL_ENABLED,
)
def search_variant_gene_fusion(search_input: SearchInputInstance):
    """ Lookup only - searching must never mint a fusion identity.

        A fusion Variant sits on the contig every build shares, so it is found once rather than per
        build, and get_visible_variants still applies each build's permissions. """
    events = find_gene_fusions_for_string(search_input.search_string)
    yield from _yield_gene_level_results(search_input, events, "Gene fusion")


@search_receiver(
    search_type=Variant,
    pattern=COPY_NUMBER_STRING_PATTERN,
    sub_name="Gene Copy Number",
    example=SearchExample(
        note="A whole-gene copy number event, named by its gene",
        examples=["EGFR amplification", "PTEN loss"]
    ),
    enabled=settings.VARIANT_GENE_LEVEL_ENABLED,
)
def search_variant_gene_copy_number(search_input: SearchInputInstance):
    """ Lookup only - searching must never mint a copy number identity. 'amp', 'gain', 'deletion'
        and 'del' all find what is written out as 'amplification' / 'loss' """
    events = find_gene_copy_number_events_for_string(search_input.search_string)
    yield from _yield_gene_level_results(search_input, events, "Gene copy number")


@search_receiver(
    search_type=Variant,
    pattern=SPLICE_STRING_PATTERN,
    sub_name="Splice Event",
    example=SearchExample(
        note="A splice event, named by its gene and the junction's label",
        examples=["AR V7", "MET exon 14 skipping"]
    ),
    enabled=settings.VARIANT_GENE_LEVEL_ENABLED,
)
def search_variant_splice_event(search_input: SearchInputInstance):
    """ Lookup only - searching must never mint a splice identity. A label shape not in
        SPLICE_STRING_PATTERN never reaches this, so widen the pattern when a new one is accepted.
        A junction named by its breakpoints is a different junction in each build, so each build the
        search covers is asked for its own label """
    events = []
    for genome_build in search_input.genome_builds:
        events.extend(find_splice_events_for_string(search_input.search_string, genome_build=genome_build))
    yield from _yield_gene_level_results(search_input, events, "Splice event")


def _yield_gene_level_results(search_input: SearchInputInstance, events, label: str):
    """ A gene-level Variant sits on the contig every build shares, so it is found once rather than
        per build, and get_visible_variants still applies each build's permissions. """
    seen = set()
    for event in events:
        if event.variant_id in seen:
            continue
        seen.add(event.variant_id)
        for genome_build in search_input.genome_builds:
            visible_variants_qs = get_visible_variants(search_input, genome_build)
            if variant := visible_variants_qs.filter(pk=event.variant_id).first():
                yield SearchResult(variant.preview,
                                   messages=[SearchMessage(f"{label} {event.canonical_str}",
                                                           severity=LogLevel.INFO)])
                break
