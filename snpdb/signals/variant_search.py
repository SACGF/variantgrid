import itertools
import logging
import re
from collections import defaultdict
from itertools import zip_longest
from typing import Optional, Iterable, Union, Callable

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.urls import reverse

from annotation.cosmic import CosmicAPI
from annotation.manual_variant_entry import check_can_create_variants, CreateManualVariantForbidden
from classification.models import Classification, CreateNoClassificationForbidden
from genes.hgvs import HGVSMatcher, HGVSException, VariantResolvingError, HGVSImplementationException, \
    HGVSNomenclatureException, HgvsOriginallyNormalized
from genes.hgvs.hgvs_converter import HgvsMatchRefAllele
from genes.models import MissingTranscript, MANE, TranscriptVersion, BadTranscript
from genes.models_enums import AnnotationConsortium
from library.enums.log_level import LogLevel
from library.genomics import format_chrom
from library.preview_request import PreviewData
from snpdb.clingen_allele import get_clingen_allele
from snpdb.models import Variant, LOCUS_PATTERN, LOCUS_NO_REF_PATTERN, DbSNP, DBSNP_PATTERN, VariantCoordinate, \
    ClinGenAllele, GenomeBuild, Contig, HGVS_UNCLEANED_PATTERN, VARIANT_PATTERN, VARIANT_SYMBOLIC_PATTERN, Allele, \
    Sequence
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, SearchResult, SearchMessageOverall, \
    SearchMessage, INVALID_INPUT
from upload.models import ModifiedImportedVariant

COSMIC_PATTERN = re.compile(r"^(COS[VM])[0-9]{3,}$", re.IGNORECASE)


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
    for genome_build in search_input.genome_builds:
        matcher = HGVSMatcher(genome_build)
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
        yield search_input.get_visible_variants(genome_build).filter(
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
        variants_qs = search_input.get_visible_variants(genome_build)
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
            visible_variants_qs = search_input.get_visible_variants(genome_build)
            if variant := allele.variant_for_build_optional(genome_build):
                yield from visible_variants_qs.filter(pk=variant.pk)
                continue
            # if here then No variant found

            # If we are non-admin we don't want to see anything else (eg ClinGen not in our build)
            if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not search_input.user.is_superuser:
                continue

            try:
                variant_coordinate = clingen_allele.get_variant_coordinate(genome_build)
                yield from _yield_no_results_for_variant_coordinate(search_input.user, genome_build, visible_variants_qs,
                                                                    variant_coordinate, [])
            except ValueError as e:
                yield SearchMessageOverall(str(e), severity=LogLevel.ERROR, genome_builds=[genome_build])


def _sv_alt_description(v: Union[Variant, VariantCoordinate]) -> str:
    return f"{v.alt}/SVLEN={v.svlen}"

def _yield_no_results_for_variant_coordinate(user, genome_build: GenomeBuild, variant_qs,
                                             variant_coordinate: VariantCoordinate,
                                             search_messages: list[SearchMessage]) -> Iterable[SearchResult]:
    # manual variants
    variant_string = variant_coordinate.format()
    if cmv := VariantExtra.create_manual_variant(for_user=user, genome_build=genome_build,
                                                 variant_string=variant_string):
        yield SearchResult(cmv, messages=search_messages)

    if variant_coordinate.is_symbolic():
        original_alt_desc = _sv_alt_description(variant_coordinate)
    else:
        original_alt_desc = Sequence.abbreviate(variant_coordinate.alt)

    # search for alt alts
    alts = get_results_from_variant_coordinate(genome_build, variant_qs, variant_coordinate, any_alt=True)
    for alternative_variant in alts:
        if alternative_variant.is_symbolic:
            alt_alt_desc = _sv_alt_description(alternative_variant)
        else:
            alt_alt_desc = Sequence.abbreviate(alternative_variant.alt)

        alt_messages = search_messages + [
            SearchMessage(f'No results for alt "{original_alt_desc}", but found this using alt "{alt_alt_desc}"',
                          severity=LogLevel.ERROR, substituted=True)]
        yield SearchResult(alternative_variant.preview, messages=alt_messages)


def get_results_from_variant_coordinate(genome_build: GenomeBuild, qs: QuerySet, variant_coordinate: VariantCoordinate,
                                        any_alt: bool = False) -> QuerySet[Variant]:
    """
    :param genome_build: genome build (used for format variant_coordinate.chromosome for variant search
    :param qs: A query set that we'll be searching inside of (except for when returning ModifiedImportVariants)
    :param variant_coordinate: The variant coordinate to lookup
    :param any_alt: If true, search without using alt and return all matches
    :return: A QuerySet of variants
    """
    chrom = format_chrom(variant_coordinate.chrom, genome_build.reference_fasta_has_chr)
    results = qs.filter(Variant.get_chrom_q(chrom),
                        locus__position=variant_coordinate.position,
                        locus__ref__seq=variant_coordinate.ref)
    if not any_alt:
        results = results.filter(alt__seq=variant_coordinate.alt, svlen=variant_coordinate.svlen)

    if not results.exists():
        if not any_alt:
            return ModifiedImportedVariant.get_variants_for_unnormalized_variant(variant_coordinate)
        else:
            # should we really be searching ModifiedImportVariants with any alt? or should that just happen for
            # the filter of "real" variants
            return ModifiedImportedVariant.get_variants_for_unnormalized_variant_any_alt(variant_coordinate)
    return results


def yield_search_variant_match(search_input: SearchInputInstance, get_variant_coordinate: Callable):
    for genome_build in search_input.genome_builds:
        variant_coordinate = get_variant_coordinate(search_input.match, genome_build)
        # Coordinate is explicit/symbolic according to length
        if errors := Variant.validate(genome_build, variant_coordinate.chrom, variant_coordinate.position):
            yield SearchMessageOverall(", ".join(errors), genome_builds=[genome_build])
            continue

        search_messages = []
        if not variant_coordinate.is_symbolic():
            # Only check refs if explicit, non-symbolic refs
            calculated_ref = variant_coordinate.calculated_reference(genome_build)
            logging.info("calc ref: %s, provided ref: %s", calculated_ref, variant_coordinate.ref)
            if calculated_ref != variant_coordinate.ref:
                provided_ref = variant_coordinate.ref
                sm = SearchMessage(f'Using genomic reference "{calculated_ref}" from our build, in place of provided reference "{provided_ref}"',
                                   LogLevel.WARNING, substituted=True)
                search_messages.append(sm)
                variant_coordinate.ref = calculated_ref

        visible_variants_qs = search_input.get_visible_variants(genome_build)
        results = get_results_from_variant_coordinate(genome_build, visible_variants_qs, variant_coordinate)
        if results.exists():
            for v in results:
                yield SearchResult(v.preview, messages=search_messages)
        else:
            yield from _yield_no_results_for_variant_coordinate(search_input.user, genome_build, visible_variants_qs,
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
        matcher = HGVSMatcher(genome_build)
        for data in dbsnp.get_alleles_for_genome_build(genome_build):
            if hgvs_string := data.get("hgvs"):
                search_message = SearchMessage(f'dbSNP "{search_input.search_string}" resolved to "{hgvs_string}"',
                                               severity=LogLevel.INFO)
                for r in _yield_results_from_hgvs(search_input, genome_build, matcher, hgvs_string):
                    yield r, search_message


def _yield_results_from_hgvs(search_input, genome_build, matcher, hgvs_string) -> Iterable[SearchResult]:
    variant_coordinate = matcher.get_variant_coordinate(hgvs_string)
    results = get_results_from_variant_coordinate(genome_build, search_input.get_visible_variants(genome_build),
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


def _search_hgvs_using_gene_symbol(
        transcript_versions,
        mane_status_by_transcript,
        hgvs_matcher: HGVSMatcher,
        search_messages: list[SearchMessage],
        hgvs_string: str,
        user: User,
        genome_build: GenomeBuild,
        variant_qs: QuerySet) -> Iterable[Union[SearchResult, SearchMessageOverall]]:

    # Group results + hgvs by result.identifier
    results_by_variant_identifier: dict[str, list[SearchResult]] = defaultdict(list)
    transcript_accessions_by_variant_identifier: dict[str, list] = defaultdict(list)
    transcript_accessions_by_exception: dict[str, list] = defaultdict(list)
    hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_string)
    hgvs_variant.gene = None

    # TODO: Sort transcript_versions somehow
    for transcript_version in transcript_versions:
        hgvs_variant.transcript = transcript_version.accession
        transcript_hgvs = hgvs_variant.format()
        tv_message = str(transcript_version.accession)
        if mane_status := mane_status_by_transcript.get(transcript_version.accession):
            tv_message += f" ({mane_status})"
        try:
            for result in _search_hgvs(transcript_hgvs, user, genome_build, variant_qs):
                if isinstance(result, SearchResult):
                    result.preview.annotation_consortia = [AnnotationConsortium(transcript_version.annotation_consortium)]
                    if result.preview.category == "Variant":
                        variant_identifier = result.preview.identifier
                        results_by_variant_identifier[variant_identifier].append(result)
                        transcript_accessions_by_variant_identifier[variant_identifier].append(tv_message)
                    else:
                        # result is a ClassifyVariantHgvs or similar, yield it and only care about real variants for the rest
                        yield result
        except Exception as e:
            transcript_accessions_by_exception[str(e)].append(tv_message)

    have_results = False
    for variant_identifier, results_for_record in results_by_variant_identifier.items():
        result_message = f"Results for: {', '.join(transcript_accessions_by_variant_identifier[variant_identifier])}"
        # Use a dict for uniqueness while preserving order
        unique_messages = {
            SearchMessage(result_message): True
        }
        # Go through messages for each result together, so they stay in same order
        for messages in zip_longest(*[r.messages for r in results_for_record]):
            for m in messages:
                if m:
                    unique_messages[m] = True

        search_messages.extend(unique_messages.keys())

        unique_annotation_consortia = set()
        for r in results_for_record:
            unique_annotation_consortia.update(r.annotation_consortia)

        # All weights should be the same, just take 1st
        have_results = True
        preview = results_for_record[0].preview  # These will all be the same
        if unique_annotation_consortia:
            preview.annotation_consortia = unique_annotation_consortia
        yield SearchResult(preview=preview, messages=search_messages)

    if transcript_accessions_by_exception:
        resolution_errors = []
        for exception_msg, tv_message_list in transcript_accessions_by_exception.items():
            resolution_errors.append(f"Error resolving transcripts: {', '.join(tv_message_list)}: {exception_msg}")

        resolution_msg = "\n".join(resolution_errors)
        if settings.SEARCH_HGVS_GENE_SYMBOL_REPORT_FAILURES:
            yield SearchResult.error_result(resolution_msg, genome_build)
        elif not have_results:
            raise VariantResolvingError(resolution_msg)

    if not have_results:
        # In some special cases, add in special messages for no result
        messages_as_strs = [str(message) for message in search_messages]
        if settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE and not settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS:
            message = "\n".join(messages_as_strs + [f"Only searched MANE transcripts: {', '.join(mane_status_by_transcript)}"])
            yield SearchMessageOverall(message)

        elif not (settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE or settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS):
            yield SearchMessageOverall("\n".join(messages_as_strs))


def _get_search_hgvs_gene_symbol_transcripts(gene_symbol, genome_build):
    transcript_versions = set()
    mane_status_by_transcript = {}
    if settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE:
        # There can be multiple MANE entries (MANE Plus Clinical and MANE Select)
        for mane in MANE.objects.filter(symbol=gene_symbol):
            for ac in AnnotationConsortium:
                if tv := mane.get_transcript_version(ac):
                    transcript_versions.add(tv)
                    mane_status_by_transcript[tv.accession] = mane.get_status_display()
    if settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS:
        for gene in gene_symbol.genes:
            tv_qs = TranscriptVersion.objects.filter(genome_build=genome_build, gene_version__gene=gene)
            # Take latest version for each transcript
            for tv in tv_qs.order_by("transcript_id", "-version").distinct("transcript_id"):
                transcript_versions.add(tv)
    return transcript_versions, mane_status_by_transcript


@search_receiver(
    search_type=Variant,
    pattern=HGVS_UNCLEANED_PATTERN,
    sub_name="HGVS",
    example=SearchExample(
        note="Provide a c. or g. HGVS",
        examples=["NM_001080463.1:c.5972T>A", "NM_000492.3(CFTR):c.1438G>T", "NC_000007.13:g.117199563G>T"]
    )
)
def search_hgvs(search_input: SearchInputInstance) -> Iterable[SearchResult]:
    if search_input.search_string.lower().startswith("medgen"):
        return [INVALID_INPUT]
    for_all_genome_builds = []
    for genome_build in search_input.genome_builds:
        results = []
        error_message = None
        try:
            results =_search_hgvs(hgvs_string=search_input.search_string, user=search_input.user,
                                  genome_build=genome_build,
                                  visible_variants=search_input.get_visible_variants(genome_build),
                                  classify=search_input.classify)
            results = list(results)  # read iterator to trigger exceptions
        except (HGVSNomenclatureException, BadTranscript) as e:
            # Show errors about HGVS the user can fix
            error_message = str(e)
        except (Exception, HGVSImplementationException) as e:
            # Hide implementation details (user can't fix)
            # log_traceback()
            # logging.info("HGVS exception type(%s): %s", type(e), str(e))
            if search_input.user.is_superuser:
                error_message = str(e)
            else:
                error_message = "Error resolving HGVS"

        if error_message:
            results = [SearchResult.error_result(error_message, genome_build)]
        if results:
            for_all_genome_builds.append(results)

    return itertools.chain(*for_all_genome_builds)


def _search_hgvs(hgvs_string: str, user: User, genome_build: GenomeBuild, visible_variants: QuerySet, classify: bool = False) -> Iterable[Union[SearchResult, SearchMessageOverall]]:
    hgvs_matcher = HGVSMatcher(genome_build)
    variant_qs = visible_variants

    # TODO, add genome build to more of the SearchMessages that are genome build specific

    # can add on search_message to objects to (stop auto-jump and show message)
    original_hgvs_string = hgvs_string
    variant_coordinate = None
    used_transcript_accession: Optional[str] = None
    kind = None
    clean_hgvs_string, _ = hgvs_matcher.clean_hgvs(hgvs_string)
    if clean_hgvs_string != hgvs_string:
        # TODO add proper support for "cleaned" where a before and after value are presented
        yield SearchMessageOverall(f'Cleaned hgvs from \n"{hgvs_string}" to\n"{clean_hgvs_string}"', severity=LogLevel.INFO)
        hgvs_string = clean_hgvs_string

    search_messages: list[SearchMessage] = []  # [SearchMessage(m) for m in hgvs_search_messages]

    try:
        vc_details = hgvs_matcher.get_variant_coordinate_and_details(hgvs_string)
        variant_coordinate, used_transcript_accession, kind, _used_converter_type, method, matches_reference, originally_normalized = vc_details
        logging.info("get_variant_coordinate_and_details - variant_coordinate=%s", variant_coordinate)
        if not matches_reference:
            # reporting on the "provided" reference is slightly problematic as it's not always provided directly, it could be indirectly

            if isinstance(matches_reference, HgvsMatchRefAllele) and matches_reference.provided_ref:
                msg = matches_reference.get_message()
                search_messages.append(SearchMessage(msg, LogLevel.WARNING, substituted=True))

        if not originally_normalized:
            if isinstance(originally_normalized, HgvsOriginallyNormalized):
                if msg := originally_normalized.get_message():
                    search_messages.append(SearchMessage(msg, LogLevel.WARNING, substituted=True))

    except MissingTranscript:
        pass
    except Contig.ContigNotInBuildError as e:
        # HGVS is valid but g.HGVS contig from a different genome build than this one
        # we don't want to throw an error for 37 contig not in 38 (if both enabled)
        # but we do want to throw one if the contig is completely unrecognised
        hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_string)
        contig_accession = hgvs_variant.transcript
        if contig := Contig.objects.filter(refseq_accession=contig_accession).first():
            # Contig is known to system
            contig_builds = contig.genomebuildcontig_set.all()
            search_builds = GenomeBuild.builds_with_annotation()
            if contig_builds.filter(genome_build__in=search_builds).exists():
                pass  # skip silently - will get in other build
            else:
                build_names = ", ".join(contig_builds.order_by("genome_build__name").values_list("genome_build__name", flat=True))
                msg = f"{contig_accession=} is from builds {build_names} which are not enabled/annotated."
                yield SearchMessageOverall(msg, severity=LogLevel.WARNING)
        else:
            known_builds = GenomeBuild.get_known_builds_comma_separated_string()
            msg = f"{contig_accession=} not in known_builds: {known_builds}"
            yield SearchMessageOverall(msg, severity=LogLevel.WARNING)

    except (ValueError, NotImplementedError) as hgvs_error:
        if gene_symbol_str := hgvs_matcher.get_gene_symbol_if_no_transcript(hgvs_string):
            has_results = False

            if settings.SEARCH_HGVS_GENE_SYMBOL:
                msg_hgvs_given_symbol = "Searching without transcript is only supported for coordinates resolved against MANE transcripts."

                mane_and_aliases = MANE.get_mane_and_aliases_list_from_symbol(gene_symbol_str)
                if not mane_and_aliases:
                    raise ValueError(
                        msg_hgvs_given_symbol + f" {gene_symbol_str} (or any aliases) have no MANE transcripts")

                gene_aliases_to_mane_symbols = set()
                genes_with_non_mane_transcripts = set()
                transcript_versions = set()
                mane_status_by_transcript = {}
                for mane, alias in mane_and_aliases:
                    if alias:
                        gene_aliases_to_mane_symbols.add(f"{alias.gene_symbol} ({alias.get_source_display()})")

                    gene_tvs, gene_mane_status = _get_search_hgvs_gene_symbol_transcripts(mane.symbol, genome_build)
                    transcript_versions.update(gene_tvs)
                    mane_status_by_transcript.update(gene_mane_status)

                    tv_qs = TranscriptVersion.objects.filter(genome_build=genome_build,
                                                             gene_version__gene__in=mane.symbol.genes)
                    transcripts = [tv.transcript for tv in transcript_versions]

                    if tv_qs.exclude(transcript__in=transcripts).distinct("transcript_id").exists():
                        # We want to make the messages the same for both genome builds, so they are collected into 1
                        genes_with_non_mane_transcripts.add(str(mane.symbol))

                # Combine messages into 1
                msg_hgvs_gene_search = msg_hgvs_given_symbol
                if gene_aliases_to_mane_symbols:
                    if len(gene_aliases_to_mane_symbols) == 1:
                        alias_word = "alias"
                    else:
                        alias_word = "aliases"
                    msg_hgvs_gene_search += f" Matched to MANE symbols via {alias_word}: {', '.join(gene_aliases_to_mane_symbols)}. "
                if genes_with_non_mane_transcripts:
                    if len(genes_with_non_mane_transcripts) == 1:
                        gene_word = "Gene symbol"
                        has_word = "has"
                    else:
                        gene_word = "Gene symbols"
                        has_word = "have"

                    msg_hgvs_gene_search += f"{gene_word} {', '.join(genes_with_non_mane_transcripts)} {has_word} "
                    msg_hgvs_gene_search += "non-MANE transcripts that may resolve to different coordinates. " + \
                        "You may wish to add a transcript or search for the gene symbol to view all results"

                yield SearchMessageOverall(msg_hgvs_gene_search, severity=LogLevel.INFO)

                try:
                    for result in _search_hgvs_using_gene_symbol(transcript_versions, mane_status_by_transcript,
                                                                 hgvs_matcher, search_messages,
                                                                 hgvs_string, user, genome_build, variant_qs):
                        if isinstance(result, SearchResult):
                            # don't count SearchMessageOverall as a result
                            has_results = True
                            # FIXME - not sure what a SearchMessage of the alias is meant to accomplish by itself?
                            # result.messages.append(SearchMessage(str(alias)))
                            yield result
                except HGVSException as e:
                    logging.error(e)
                    # might just not be a HGVS name at all
            else:
                yield SearchMessageOverall("HGVS search requires a transcript")

            if has_results:
                return
            else:  # Valid - just no results
                hgvs_error = None

        if variant_coordinate is None:
            if classify:
                search_message = SearchMessage(f"Error reading HGVS \"{hgvs_error}\"")
                if classify_no_variant := VariantExtra.classify_no_variant_hgvs(for_user=user, genome_build=genome_build, hgvs_string=original_hgvs_string):
                    yield SearchResult(classify_no_variant, messages=[search_message])

        # yield SearchMessageOverall(str(hgvs_error))
        if hgvs_error:  # May have been cleared if matched gene symbol HGVS
            raise hgvs_error

    if used_transcript_accession:
        if used_transcript_accession not in hgvs_string:
            reported = False
            # are we just inserting a version because none was provided, or are we changing it

            if '.' in used_transcript_accession:
                used_transcript_parts = used_transcript_accession.split('.')
                used_transcript_versionless = used_transcript_parts[0]
                used_transcript_version = used_transcript_parts[1]
                if (used_transcript_versionless + '.') in hgvs_string:
                    transcript_index = hgvs_string.index(used_transcript_versionless)
                    transcript_version_index = transcript_index + len(used_transcript_versionless) + 1
                    old_transcript_version = ''
                    for c in hgvs_string[transcript_version_index:]:
                        if '0' <= c <= '9':
                            old_transcript_version += c
                        else:
                            break
                    search_messages.append(SearchMessage(
                        f"Using transcript \"{used_transcript_versionless}\" version \"{used_transcript_version}\" instead of provided version \"{old_transcript_version}\"",
                        severity=LogLevel.ERROR, substituted=True))
                    reported = True

            if not reported:
                search_messages.append(SearchMessage(f"Used transcript version \"{used_transcript_accession}\"", severity=LogLevel.INFO))

        hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_string)
        # If these were in wrong order they have been switched now
        if hgvs_variant.transcript and hgvs_variant.gene:
            annotation_consortium = AnnotationConsortium.get_from_transcript_accession(used_transcript_accession)
            try:
                # If we look it up using a different method (eg ClinGen Allele registry) the used transcript accession
                # May not exist on our system. That's ok
                transcript_version = TranscriptVersion.get(used_transcript_accession, genome_build,
                                                           annotation_consortium=annotation_consortium)
                alias_symbol_strs = transcript_version.gene_version.gene_symbol.alias_meta.alias_symbol_strs
                if hgvs_variant.gene.upper() not in [a.upper() for a in alias_symbol_strs]:
                    search_messages.append(SearchMessage(f"Symbol \"{hgvs_variant.gene}\" not associated with transcript "
                                           f"{used_transcript_accession} (known symbols='{', '.join(alias_symbol_strs)}')"))
            except TranscriptVersion.DoesNotExist:
                pass  # Doesn't exist - OK just can't check

    # TODO: alter initial_score based on warning messages of alt not matching?
    # also - _lrg_get_variant_coordinate should add matches_reference to search warnings list

    if variant_coordinate:
        try:
            results = get_results_from_variant_coordinate(genome_build, variant_qs, variant_coordinate)
            variant = results.get()
            if classify:
                yield VariantExtra.classify_variant(
                    variant=variant,
                    genome_build=genome_build
                ), search_messages
            else:
                # doesn't matter what the preferred genome build is (explicit contig version)
                is_genomic = kind in ('g', 'm')
                yield SearchResult(preview=variant.preview, messages=search_messages,
                                   ignore_genome_build_mismatch=is_genomic)

        except Variant.DoesNotExist:
            variant_string_abbreviated = Variant.format_tuple(*variant_coordinate, abbreviate=True)
            search_messages.append(
                SearchMessage(f'"{hgvs_string}" resolved to genomic coordinates "{variant_string_abbreviated}"',
                              LogLevel.INFO,
                              genome_build=genome_build))
            # if we're saying what we're resolving to, no need to complain about reference_message

            yield from _yield_no_results_for_variant_coordinate(user, genome_build, variant_qs,
                                                                variant_coordinate, search_messages)


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
    yield Variant.objects.filter(pk=search_input.match.group(2))


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
    yield Allele.objects.filter(pk=search_input.match.group(1))
