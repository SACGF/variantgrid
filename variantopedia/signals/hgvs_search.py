""" Variant search by HGVS (the search_hgvs receiver): resolves c./n./g. HGVS through HGVSMatcher, retrying c. on a
    non-coding transcript as n., falls back to a gene symbol's MANE transcripts when no transcript is given, and
    reports what was substituted (transcript version, reference, normalisation) as search messages """
import itertools
import logging
import re
import sys
from collections import defaultdict
from collections.abc import Iterable
from itertools import zip_longest
from typing import Optional, Union

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.template.defaultfilters import pluralize

from genes.enums import AnnotationConsortium, MANEStatus
from genes.hgvs import (
    HGVSException,
    HGVSImplementationException,
    HGVSMatcher,
    HGVSNomenclatureException,
    HgvsOriginallyNormalized,
    VariantResolvingError,
)
from genes.hgvs.hgvs_converter import HgvsMatchRefAllele, HGVSNonCodingTranscriptException
from genes.hgvs.hgvs_matcher import VariantCoordinateAndDetails
from genes.models import MANE, BadTranscript, MissingTranscript, TranscriptVersion
from library.enums.log_level import LogLevel
from library.log_utils import report_exc_info
from snpdb.models import (
    HGVS_UNCLEANED_PATTERN,
    Contig,
    GenomeBuild,
    Variant,
)
from snpdb.search import (
    INVALID_INPUT,
    SearchExample,
    SearchInputInstance,
    SearchMessage,
    SearchMessageOverall,
    SearchResult,
    search_receiver,
)
from snpdb.variant_coordinate_search import get_results_from_variant_coordinate
from variantopedia.signals.variant_search import (
    VariantExtra,
    get_visible_variants,
    yield_no_results_for_variant_coordinate,
)


def _search_hgvs_using_gene_symbol(
        candidates: dict[str, Optional[str]],
        allele: str,
        user: User,
        genome_build: GenomeBuild,
        variant_qs: QuerySet) -> Iterable[Union[SearchResult, SearchMessageOverall]]:
    """ candidates: {transcript_accession: mane_status_or_None}, best-first, sourced from cdot.
        Each accession + allele is searched via the standard _search_hgvs path (which
        resolves the transcript version), and hits are grouped/deduped by variant. """

    # Group results + hgvs by result.identifier
    results_by_variant_identifier: dict[str, list[SearchResult]] = defaultdict(list)
    transcript_accessions_by_variant_identifier: dict[str, list] = defaultdict(list)
    transcript_accessions_by_exception: dict[str, list] = defaultdict(list)

    for accession, mane_status in candidates.items():
        transcript_hgvs = f"{accession}:{allele}"
        tv_message = accession
        if mane_status:
            tv_message += f" ({mane_status})"
        try:
            for result in _search_hgvs(transcript_hgvs, user, genome_build, variant_qs):
                if isinstance(result, SearchResult):
                    result.preview.annotation_consortia = [AnnotationConsortium.get_from_transcript_accession(accession)]
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

        unique_annotation_consortia = set()
        for r in results_for_record:
            unique_annotation_consortia.update(r.annotation_consortia)

        # All weights should be the same, just take 1st
        have_results = True
        preview = results_for_record[0].preview  # These will all be the same
        if unique_annotation_consortia:
            preview.annotation_consortia = unique_annotation_consortia
        yield SearchResult(preview=preview, messages=list(unique_messages))

    if transcript_accessions_by_exception:
        resolution_errors = []
        for exception_msg, tv_message_list in transcript_accessions_by_exception.items():
            resolution_errors.append(f"Error resolving transcripts: {', '.join(tv_message_list)}: {exception_msg}")

        resolution_msg = "\n".join(resolution_errors)
        if settings.SEARCH_HGVS_GENE_SYMBOL_REPORT_FAILURES:
            yield SearchResult.error_result(resolution_msg, genome_build)
        elif not have_results:
            raise VariantResolvingError(resolution_msg)

    if not have_results and settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE and not settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS:
        yield SearchMessageOverall(f"Only searched MANE transcripts: {', '.join(candidates)}")


# cdot MANE tag spellings -> MANEStatus display strings (for per-transcript messages)
_CDOT_TAG_TO_MANE_DISPLAY = {
    "MANE_Select": MANEStatus.MANE_SELECT.label,
    "MANE_Plus_Clinical": MANEStatus.MANE_PLUS_CLINICAL.label,
}


def _get_search_hgvs_gene_symbol_transcripts(hgvs_matcher: HGVSMatcher, gene_symbol_str: str):
    """ Source candidate transcripts for a gene symbol from cdot's ranking
        (MANE/canonical tags), honouring VG's MANE / all-transcripts settings.

        Returns (candidates, has_non_mane_transcripts) where candidates is a list of
        (transcript_accession, mane_status_or_None) best-first (cdot's versionless
        accessions), and has_non_mane_transcripts is True if the gene has non-MANE
        transcripts that were not included (used to warn they may resolve differently). """
    candidates = []
    has_non_mane_transcripts = False

    ranked, _warnings = hgvs_matcher.rank_gene_symbol_transcripts(gene_symbol_str)
    for tx_ac, _tags, matched_tag in ranked:
        mane_status = _CDOT_TAG_TO_MANE_DISPLAY.get(matched_tag)
        is_mane = mane_status is not None
        include = settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS or \
            (is_mane and settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE)
        if not include:
            if not is_mane:
                has_non_mane_transcripts = True
            continue
        candidates.append((tx_ac, mane_status))
    return candidates, has_non_mane_transcripts


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
                                  visible_variants=get_visible_variants(search_input, genome_build),
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
                report_exc_info()
                # error_message = "Error resolving HGVS"
                error_message = str(e)
                # work has been done to improve the error messages
                # planning a future change where the technical errors identify themselves rather than the
                # user-friendly errors identifying themselves

        if error_message:
            results = [SearchResult.error_result(error_message, genome_build)]
            if search_input.classify:
                results.extend(_classify_unresolved_hgvs(search_input.user, genome_build,
                                                         search_input.search_string, error_message))
        if results:
            for_all_genome_builds.append(results)

    return itertools.chain(*for_all_genome_builds)


def _classify_unresolved_hgvs(user: User, genome_build: GenomeBuild, hgvs_string: str, error_message: str) -> Iterable[SearchResult]:
    if classify_no_variant := VariantExtra.classify_no_variant_hgvs(for_user=user, genome_build=genome_build, hgvs_string=hgvs_string):
        yield SearchResult(classify_no_variant, messages=[SearchMessage(f"Error reading HGVS \"{error_message}\"")])


def _get_non_coding_c_hgvs_as_n_details(hgvs_matcher: HGVSMatcher, hgvs_string: str,
                                        search_messages: list[SearchMessage]) -> VariantCoordinateAndDetails:
    """ Import refuses c.HGVS on a non-coding transcript (it may be numbered from a coding model, so shifted by
        the 5' UTR length). Search resolves it as n. with a warning, unless the provided reference doesn't match """
    hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_string)
    hgvs_variant.kind = 'n'
    n_hgvs_string = hgvs_variant.format(max_ref_length=sys.maxsize)
    vc_details = hgvs_matcher.get_variant_coordinate_and_details(n_hgvs_string)
    matches_reference = vc_details.matches_reference
    if isinstance(matches_reference, HgvsMatchRefAllele) and not matches_reference:
        raise HGVSNomenclatureException(f'c.HGVS used on non-coding transcript. Not resolved as "{n_hgvs_string}" '
                                        'as the provided reference does not match, so it is probably numbered '
                                        'from a coding transcript')
    msg = f'c.HGVS used on non-coding transcript, resolved as "{n_hgvs_string}". If the c. was numbered from a ' \
          "coding transcript, the position will be shifted by the length of its 5' UTR"
    search_messages.append(SearchMessage(msg, LogLevel.WARNING, substituted=True))
    return vc_details


def _resolve_hgvs(hgvs_matcher: HGVSMatcher, hgvs_string: str) -> tuple[VariantCoordinateAndDetails, list[SearchMessage]]:
    """ Resolves hgvs_string to a coordinate, with messages for what the matcher had to change to get there.
        Raises what the matcher raises: MissingTranscript, Contig.ContigNotInBuildError, HGVSNomenclatureException,
        and ValueError / NotImplementedError for anything else (including a gene symbol in place of a transcript) """
    search_messages: list[SearchMessage] = []
    try:
        vc_details = hgvs_matcher.get_variant_coordinate_and_details(hgvs_string)
    except HGVSNonCodingTranscriptException:
        vc_details = _get_non_coding_c_hgvs_as_n_details(hgvs_matcher, hgvs_string, search_messages)
    logging.info("get_variant_coordinate_and_details - variant_coordinate=%s", vc_details.variant_coordinate)

    # reporting on the "provided" reference is slightly problematic as it's not always provided directly
    matches_reference = vc_details.matches_reference
    if not matches_reference and isinstance(matches_reference, HgvsMatchRefAllele) and matches_reference.provided_ref:
        search_messages.append(SearchMessage(matches_reference.get_message(), LogLevel.WARNING, substituted=True))

    originally_normalized = vc_details.originally_normalized
    if not originally_normalized and isinstance(originally_normalized, HgvsOriginallyNormalized):
        if msg := originally_normalized.get_message():
            search_messages.append(SearchMessage(msg, LogLevel.WARNING, substituted=True))
    return vc_details, search_messages


def _contig_not_in_build_message(contig_accession: str) -> Optional[SearchMessageOverall]:
    """ g.HGVS on a contig from another build. Silent if that build is searched too (a 37 contig while searching 38,
        with both enabled), a warning if the contig's builds aren't annotated or the contig is unknown """
    if contig := Contig.objects.filter(refseq_accession=contig_accession).first():
        contig_builds = contig.genomebuildcontig_set.all()
        if contig_builds.filter(genome_build__in=GenomeBuild.builds_with_annotation()).exists():
            return None
        build_names = ", ".join(contig_builds.order_by("genome_build__name").values_list("genome_build__name", flat=True))
        msg = f"{contig_accession=} is from builds {build_names} which are not enabled/annotated."
    else:
        enabled_builds = GenomeBuild.get_enabled_builds_comma_separated_string()
        msg = f"{contig_accession=} not in enabled builds: {enabled_builds}"
    return SearchMessageOverall(msg, severity=LogLevel.WARNING)


_GENE_SYMBOL_HGVS_NEEDS_MANE = "Searching without transcript is only supported for coordinates resolved against MANE transcripts."


def _gene_symbol_search_message(gene_aliases_to_mane_symbols: set[str], genes_with_non_mane_transcripts: set[str]) -> str:
    """ One message for both genome builds, so they collapse into one on the search page """
    gene_search_sentences = [_GENE_SYMBOL_HGVS_NEEDS_MANE]
    if gene_aliases_to_mane_symbols:
        gene_search_sentences.append(f"Matched to MANE symbols via alias{pluralize(gene_aliases_to_mane_symbols, 'es')}: "
                                     f"{', '.join(gene_aliases_to_mane_symbols)}.")
    if genes_with_non_mane_transcripts:
        gene_search_sentences.append(
            f"Gene symbol{pluralize(genes_with_non_mane_transcripts)} {', '.join(genes_with_non_mane_transcripts)} "
            f"{pluralize(genes_with_non_mane_transcripts, 'has,have')} "
            "non-MANE transcripts that may resolve to different coordinates. "
            "You may wish to add a transcript or search for the gene symbol to view all results")
    return " ".join(gene_search_sentences)


def _search_hgvs_without_transcript(hgvs_matcher: HGVSMatcher, hgvs_string: str, gene_symbol_str: str,
                                    hgvs_error: Exception, user: User, genome_build: GenomeBuild,
                                    variant_qs: QuerySet) -> Iterable[Union[SearchResult, SearchMessageOverall]]:
    """ "SYMBOL:allele" HGVS - searched against the MANE transcripts of the symbol (and of symbols it's an alias of) """
    if not settings.SEARCH_HGVS_GENE_SYMBOL:
        yield SearchMessageOverall("HGVS search requires a transcript")
        return

    mane_and_aliases = MANE.get_mane_and_aliases_list_from_symbol(gene_symbol_str)
    if not mane_and_aliases:
        raise ValueError(
            _GENE_SYMBOL_HGVS_NEEDS_MANE + f" {gene_symbol_str} (or any aliases) have no MANE transcripts") from hgvs_error

    gene_aliases_to_mane_symbols = set()
    genes_with_non_mane_transcripts = set()
    # {accession: mane_status} best-first, deduped across aliases
    candidates: dict[str, Optional[str]] = {}
    for mane, alias in mane_and_aliases:
        if alias:
            gene_aliases_to_mane_symbols.add(f"{alias.gene_symbol} ({alias.get_source_display()})")

        gene_candidates, has_non_mane = _get_search_hgvs_gene_symbol_transcripts(hgvs_matcher, str(mane.symbol))
        for accession, mane_status in gene_candidates:
            candidates.setdefault(accession, mane_status)
        if has_non_mane:
            genes_with_non_mane_transcripts.add(str(mane.symbol))

    yield SearchMessageOverall(_gene_symbol_search_message(gene_aliases_to_mane_symbols, genes_with_non_mane_transcripts),
                               severity=LogLevel.INFO)

    try:
        # gene-only HGVS is "SYMBOL:allele" - swap each cdot transcript in for the symbol
        allele = hgvs_string.split(":", 1)[1]
        for result in _search_hgvs_using_gene_symbol(candidates, allele, user, genome_build, variant_qs):
            if isinstance(result, SearchResult):
                yield result
    except HGVSException as e:
        logging.error(e)
        # might just not be a HGVS name at all


def _hgvs_nomenclature_error_message(hgvs_string: str, raw_message: str) -> str:
    if "char" not in raw_message and "EOF" not in raw_message:
        # this is likely not a technical message, send it directly to the user
        return raw_message
    if match := HGVSMatcher.HGVS_CLEAN_PATTERN_OPTIONAL_TRANSCRIPT.match(hgvs_string):
        return f"Invalid HGVS cDNA allele : \"{match.group('letter')}.{match.group('nomen')}\""
    return f"Invalid HGVS cDNA allele \"{hgvs_string}\""


def _transcript_substitution_messages(hgvs_string: str, used_transcript_accession: str) -> list[SearchMessage]:
    """ Says which transcript version was used, when it's not the one written in hgvs_string """
    if used_transcript_accession in hgvs_string:
        return []
    # are we just inserting a version because none was provided, or are we changing it
    if '.' in used_transcript_accession:
        used_transcript_versionless, used_transcript_version = used_transcript_accession.split('.')[:2]
        if match := re.search(rf"{re.escape(used_transcript_versionless)}\.(\d*)", hgvs_string):
            old_transcript_version = match.group(1)
            return [SearchMessage(
                f"Using transcript \"{used_transcript_versionless}\" version \"{used_transcript_version}\" instead of provided version \"{old_transcript_version}\"",
                severity=LogLevel.ERROR, substituted=True)]
    return [SearchMessage(f"Used transcript version \"{used_transcript_accession}\"", severity=LogLevel.INFO)]


def _transcript_gene_symbol_messages(hgvs_matcher: HGVSMatcher, genome_build: GenomeBuild, hgvs_string: str,
                                     used_transcript_accession: str) -> list[SearchMessage]:
    """ Warns when the gene symbol given with the transcript isn't one of the transcript's gene's symbols """
    hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_string)
    # If these were in wrong order they have been switched now
    if not (hgvs_variant.transcript and hgvs_variant.gene):
        return []
    annotation_consortium = AnnotationConsortium.get_from_transcript_accession(used_transcript_accession)
    try:
        # If we look it up using a different method (eg ClinGen Allele registry) the used transcript accession
        # May not exist on our system. That's ok
        transcript_version = TranscriptVersion.get(used_transcript_accession, genome_build,
                                                   annotation_consortium=annotation_consortium)
    except TranscriptVersion.DoesNotExist:
        return []  # Doesn't exist - OK just can't check
    alias_symbol_strs = transcript_version.gene_version.gene_symbol.alias_meta.alias_symbol_strs
    if hgvs_variant.gene.upper() in [a.upper() for a in alias_symbol_strs]:
        return []
    return [SearchMessage(f"Symbol \"{hgvs_variant.gene}\" not associated with transcript "
                          f"{used_transcript_accession} (known symbols='{', '.join(alias_symbol_strs)}')")]


def _results_for_coordinate(hgvs_matcher: HGVSMatcher, hgvs_string: str, vc_details: VariantCoordinateAndDetails,
                            search_messages: list[SearchMessage], user: User, genome_build: GenomeBuild,
                            variant_qs: QuerySet, classify: bool) -> Iterable[Union[SearchResult, tuple]]:
    variant_coordinate = vc_details.variant_coordinate
    try:
        variant = get_results_from_variant_coordinate(genome_build, variant_qs, variant_coordinate).get()
    except Variant.DoesNotExist:
        variant_string_abbreviated = Variant.format_tuple(*variant_coordinate, abbreviate=True)
        # if we're saying what we're resolving to, no need to complain about reference_message
        search_messages.append(
            SearchMessage(f'"{hgvs_string}" resolved to genomic coordinates "{variant_string_abbreviated}"',
                          LogLevel.INFO,
                          genome_build=genome_build))
        yield from yield_no_results_for_variant_coordinate(user, genome_build, variant_qs,
                                                           variant_coordinate, search_messages,
                                                           hgvs_matcher=hgvs_matcher)
        return

    if classify:
        yield VariantExtra.classify_variant(variant=variant, genome_build=genome_build), search_messages
    else:
        # doesn't matter what the preferred genome build is (explicit contig version)
        is_genomic = vc_details.kind in ('g', 'm')
        yield SearchResult(preview=variant.preview, messages=search_messages,
                           ignore_genome_build_mismatch=is_genomic)


def _search_hgvs(hgvs_string: str, user: User, genome_build: GenomeBuild, visible_variants: QuerySet, classify: bool = False) -> Iterable[Union[SearchResult, SearchMessageOverall]]:
    hgvs_matcher = HGVSMatcher.instance(genome_build)

    # TODO, add genome build to more of the SearchMessages that are genome build specific
    original_hgvs_string = hgvs_string
    clean_hgvs_string, _ = hgvs_matcher.clean_hgvs(hgvs_string)
    if clean_hgvs_string != hgvs_string:
        # TODO add proper support for "cleaned" where a before and after value are presented
        yield SearchMessageOverall(f'Cleaned hgvs from \n"{hgvs_string}" to\n"{clean_hgvs_string}"', severity=LogLevel.INFO)
        hgvs_string = clean_hgvs_string

    try:
        vc_details, search_messages = _resolve_hgvs(hgvs_matcher, hgvs_string)
    except MissingTranscript:
        return
    except Contig.ContigNotInBuildError:
        contig_accession = hgvs_matcher.create_hgvs_variant(hgvs_string).transcript
        if message := _contig_not_in_build_message(contig_accession):
            yield message
        return
    except (ValueError, NotImplementedError) as hgvs_error:
        gene_symbol_str = hgvs_matcher.get_gene_symbol_if_no_transcript(hgvs_string)
        if not gene_symbol_str:
            raise
        yield from _search_hgvs_without_transcript(hgvs_matcher, hgvs_string, gene_symbol_str, hgvs_error,
                                                   user, genome_build, visible_variants)
        return
    except HGVSNomenclatureException as hgvs_ex:
        raw_message = str(hgvs_ex)
        if classify:
            yield from _classify_unresolved_hgvs(user, genome_build, original_hgvs_string, raw_message)
        yield SearchMessageOverall(_hgvs_nomenclature_error_message(hgvs_string, raw_message))
        return

    if used_transcript_accession := vc_details.transcript_accession:
        search_messages.extend(_transcript_substitution_messages(hgvs_string, used_transcript_accession))
        search_messages.extend(_transcript_gene_symbol_messages(hgvs_matcher, genome_build, hgvs_string,
                                                                used_transcript_accession))

    # TODO: alter initial_score based on warning messages of alt not matching?
    if vc_details.variant_coordinate:
        yield from _results_for_coordinate(hgvs_matcher, hgvs_string, vc_details, search_messages,
                                           user, genome_build, visible_variants, classify)
