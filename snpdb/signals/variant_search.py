import itertools
import logging
import re
from collections import defaultdict
from itertools import zip_longest
from typing import Optional, List, Iterable, Union, Dict

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.urls import reverse

from annotation.manual_variant_entry import check_can_create_variants, CreateManualVariantForbidden
from classification.models import Classification, CreateNoClassificationForbidden
from genes.hgvs import HGVSMatcher, HGVSException
from genes.hgvs.hgvs_converter import HgvsMatchRefAllele
from genes.models import MissingTranscript, MANE, TranscriptVersion, GeneSymbol
from genes.models_enums import AnnotationConsortium
from library.enums.log_level import LogLevel
from library.genomics import format_chrom
from library.preview_request import PreviewData
from snpdb.clingen_allele import get_clingen_allele
from snpdb.models import Variant, LOCUS_PATTERN, LOCUS_NO_REF_PATTERN, DbSNP, DBSNP_PATTERN, VariantCoordinate, \
    ClinGenAllele, GenomeBuild, Contig, HGVS_UNCLEANED_PATTERN
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, SearchResult, SearchMessageOverall, \
    SearchMessage
from upload.models import ModifiedImportedVariant

COSMIC_PATTERN = re.compile(r"^(COS[VM])[0-9]+$", re.IGNORECASE)


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
        examples=["COSV53567516"]
    )
)
def variant_cosmic_search(search_input: SearchInputInstance):
    search_string = search_input.search_string.upper()
    for genome_build in search_input.genome_builds:
        variant_qs = search_input.get_visible_variants(genome_build)
        if search_input.match.group(1).upper() == "COSV":
            yield variant_qs.filter(variantannotation__cosmic_id=search_string)
        elif search_input.match.group(1).upper() == "COSM":
            yield variant_qs.filter(variantannotation__cosmic_legacy_id=search_string)


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
        yield search_input.get_visible_variants(genome_build).filter(
            locus__contig__name=chrom,
            locus__position=position,
            locus__ref__seq=ref
        )


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
            try:
                if variant := allele.variant_for_build(genome_build):
                    yield variant
                    continue
            except ValueError:
                pass
            # if there was no variant for that allele
            variant_string = clingen_allele.get_variant_string(genome_build)
            if create_manual := VariantExtra.create_manual_variant(search_input.user, genome_build=genome_build, variant_string=variant_string):
                variant_string_abbreviated = clingen_allele.get_variant_string(genome_build, abbreviate=True)
                yield create_manual, SearchMessage(f'"{clingen_allele}" resolved to "{variant_string_abbreviated}"', severity=LogLevel.INFO)


def get_results_from_variant_tuples(qs: QuerySet, data: VariantCoordinate, any_alt: bool = False) -> QuerySet[Variant]:
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
            return ModifiedImportedVariant.get_variants_for_unnormalized_variant(chrom, position, ref, alt)
        else:
            # should we really be searching ModifiedImportVariants with any alt? or should that just happen for
            # the filter of "real" variants
            return ModifiedImportedVariant.get_variants_for_unnormalized_variant_any_alt(chrom, position, ref)
    return results


def yield_search_variant_match(search_input: SearchInputInstance):
    for genome_build in search_input.genome_builds:
        chrom, position, ref, alt = search_input.match.groups()
        chrom, position, ref, alt = Variant.clean_variant_fields(chrom, position, ref, alt,
                                                                 want_chr=genome_build.reference_fasta_has_chr)
        has_results = False
        results = get_results_from_variant_tuples(search_input.get_visible_variants(genome_build), VariantCoordinate(chrom, position, ref, alt))
        if results.exists():
            has_results = True
            yield results

        if errors := Variant.validate(genome_build, chrom, position):
            yield SearchMessageOverall(", ".join(errors), genome_builds=[genome_build])
        else:
            if not has_results:
                variant_string = Variant.format_tuple(chrom, position, ref, alt)
                if create_manual := VariantExtra.create_manual_variant(search_input.user, genome_build=genome_build, variant_string=variant_string):
                    search_message = f"The variant {variant_string} does not exist in our database"
                    yield create_manual, search_message


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
    return yield_search_variant_match(search_input)


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
    return yield_search_variant_match(search_input)


VARIANT_PATTERN = re.compile(r"^(MT|(?:chr)?(?:[XYM]|\d+))\s*:\s*(\d+)[,\s]*([GATC]+)>(=|[GATC]+)$", re.IGNORECASE)


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
    return yield_search_variant_match(search_input)


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
                dbsnp_message = SearchMessage(f'dbSNP "{search_input.search_string}" resolved to "{hgvs_string}"', severity=LogLevel.INFO)
                variant_tuple = matcher.get_variant_tuple(hgvs_string)
                results = get_results_from_variant_tuples(search_input.get_visible_variants(genome_build), variant_tuple)
                if results.exists():
                    for r in results:
                        yield r, dbsnp_message
                else:
                    variant_string = Variant.format_tuple(*variant_tuple)
                    if create_manual := VariantExtra.create_manual_variant(
                        for_user=search_input.user,
                        variant_string=variant_string,
                        genome_build=genome_build
                    ):
                        yield create_manual, dbsnp_message


def _search_hgvs_using_gene_symbol(
        hgvs_matcher: HGVSMatcher,
        gene_symbol: GeneSymbol,
        search_messages: List[SearchMessage],
        hgvs_string: str,
        user: User,
        genome_build: GenomeBuild,
        variant_qs: QuerySet) -> Iterable[Union[SearchResult, SearchMessageOverall]]:

    search_messages.append(SearchMessage(f'HGVS requires transcript, given symbol "{gene_symbol}"'))
    # Group results + hgvs by result.identifier
    results_by_variant_identifier: Dict[str, List[SearchResult]] = defaultdict(list)
    transcript_accessions_by_variant_identifier: Dict[str, List] = defaultdict(list)
    hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_string)
    hgvs_variant.gene = None

    transcript_versions = set()
    mane_transcripts = set()
    other_transcripts_message = None  # Want this to be after transcripts used message

    if settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE:
        if mane := MANE.objects.filter(symbol=gene_symbol).first():
            for ac in AnnotationConsortium:
                if tv := mane.get_transcript_version(ac):
                    transcript_versions.add(tv)
                    mane_transcripts.add(tv.accession)

    if settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS:
        for gene in gene_symbol.genes:
            tv_qs = TranscriptVersion.objects.filter(genome_build=genome_build, gene_version__gene=gene)
            # Take latest version for each transcript
            for tv in tv_qs.order_by("transcript_id", "-version").distinct("transcript_id"):
                transcript_versions.add(tv)
    else:
        tv_qs = TranscriptVersion.objects.filter(genome_build=genome_build, gene_version__gene__in=gene_symbol.genes)
        transcripts = [tv.transcript for tv in transcript_versions]
        if num_other := tv_qs.exclude(transcript__in=transcripts).distinct("transcript_id").count():
            other_transcripts_message = f"Warning: {gene_symbol} has {num_other} other transcripts in " \
                                        f"{genome_build.name} that may resolve to different coordinates. " \
                                        f"You may wish to search for the gene symbol to view all results"

    # TODO: Sort transcript_versions somehow
    for transcript_version in transcript_versions:
        hgvs_variant.transcript = transcript_version.accession
        transcript_hgvs = hgvs_variant.format()
        try:
            for result in _search_hgvs(transcript_hgvs, user, genome_build, variant_qs):
                if isinstance(result, SearchResult):
                    result.preview.annotation_consortia = [AnnotationConsortium(transcript_version.annotation_consortium)]
                    if result.preview.category == "Variant":
                        variant_identifier = result.preview.identifier
                        results_by_variant_identifier[variant_identifier].append(result)
                        tv_message = str(transcript_version.accession)
                        if transcript_version.accession in mane_transcripts:
                            tv_message += " (MANE)"
                        transcript_accessions_by_variant_identifier[variant_identifier].append(tv_message)
                    else:
                        # result is a ClassifyVariantHgvs or similar, yield it and only care about real variants for the rest
                        yield result
        except Exception as e:
            # Just swallow all these errors
            logging.warning(e)

    have_results = False
    for variant_identifier, results_for_record in results_by_variant_identifier.items():
        result_message = f"Results for: {', '.join(transcript_accessions_by_variant_identifier[variant_identifier])}"
        # Use a dict for uniqueness while preserving order
        unique_messages = {
            SearchMessage(result_message): True
        }
        if other_transcripts_message:
            unique_messages[SearchMessage(other_transcripts_message)] = True
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

    if not have_results:
        # In some special cases, add in special messages for no result
        messages_as_strs = [str(message) for message in search_messages]
        if settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE and not settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS:
            message = "\n".join(messages_as_strs + [f"Only searched MANE transcripts: {', '.join(mane_transcripts)}"])
            yield SearchMessageOverall(message)

        elif not (settings.SEARCH_HGVS_GENE_SYMBOL_USE_MANE or settings.SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS):
            yield SearchMessageOverall("\n".join(messages_as_strs))


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
    for_all_genome_builds = []
    for genome_build in search_input.genome_builds:
        for_all_genome_builds.append(_search_hgvs(hgvs_string=search_input.search_string, user=search_input.user, genome_build=genome_build, visible_variants=search_input.get_visible_variants(genome_build), classify=search_input.classify))
    return itertools.chain(*for_all_genome_builds)


def _search_hgvs(hgvs_string: str, user: User, genome_build: GenomeBuild, visible_variants: QuerySet, classify: bool = False) -> Iterable[Union[SearchResult, SearchMessageOverall]]:
    hgvs_matcher = HGVSMatcher(genome_build)
    variant_qs = visible_variants

    # TODO, add genome build to more of the SearchMessages that are genome build specific

    # can add on search_message to objects to (stop auto-jump and show message)
    original_hgvs_string = hgvs_string
    variant_tuple = None
    used_transcript_accession: Optional[str] = None
    kind = None
    clean_hgvs_string, _ = hgvs_matcher.clean_hgvs(hgvs_string)
    if clean_hgvs_string != hgvs_string:
        # TODO add proper support for "cleaned" where a before and after value are presented
        yield SearchMessageOverall(f'Cleaned hgvs from \n"{hgvs_string}" to\n"{clean_hgvs_string}"', severity=LogLevel.INFO)
        hgvs_string = clean_hgvs_string

    search_messages: List[SearchMessage] = [] # [SearchMessage(m) for m in hgvs_search_messages]
    reference_message: List[SearchMessage] = []

    try:
        variant_tuple, used_transcript_accession, kind, method, matches_reference = hgvs_matcher.get_variant_tuple_used_transcript_kind_method_and_matches_reference(hgvs_string)
        if not matches_reference:
            ref_base = variant_tuple[2]

            # reporting on the "provided" reference is slightly promblematic as it's not always provided directly, it could be indirectly

            if isinstance(matches_reference, HgvsMatchRefAllele) and matches_reference.provided_ref:
                reference_message.append(SearchMessage(f'Using reference "{ref_base}" from our build, in place of provided reference "{matches_reference.provided_ref}"', LogLevel.ERROR, substituted=True))
            else:
                # if no reference was provided, do we even need to provide a message?
                # e.g. this is providing a ref for when we have a delins, e.g. delinsGT => delCCinsGT
                # reference_message.append(SearchMessage(f'Using reference "{ref_base}" from our build', LogLevel.INFO))
                pass

    except (MissingTranscript, Contig.ContigNotInBuildError):
        # contig triggered from g.HGVS from another genome build - can't do anything just return no results
        pass
    except (ValueError, NotImplementedError) as hgvs_error:
        try:
            if gene_symbol := hgvs_matcher.get_gene_symbol_if_no_transcript(hgvs_string):
                has_results = False
                for result in _search_hgvs_using_gene_symbol(hgvs_matcher, gene_symbol, search_messages,
                                                             hgvs_string, user, genome_build, variant_qs):
                    if isinstance(result, SearchResult):
                        # don't count SearchMessageOverall as a result
                        has_results = True
                        yield result
                if has_results:
                    return

        except HGVSException:
            # might just not be a HGVS name at all
            pass

        if variant_tuple is None:
            if classify:
                search_message = SearchMessage(f"Error reading HGVS \"{hgvs_error}\"")
                if classify_no_variant := VariantExtra.classify_no_variant_hgvs(for_user=user, genome_build=genome_build, hgvs_string=original_hgvs_string):
                    yield SearchResult(classify_no_variant, messages=[search_message])

        # yield SearchMessageOverall(str(hgvs_error))
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
            transcript_version = TranscriptVersion.get(used_transcript_accession, genome_build,
                                                       annotation_consortium=annotation_consortium)
            alias_symbol_strs = transcript_version.gene_version.gene_symbol.alias_meta.alias_symbol_strs
            if hgvs_variant.gene.upper() not in [a.upper() for a in alias_symbol_strs]:
                search_messages.append(SearchMessage(f"Symbol \"{hgvs_variant.gene}\" not associated with transcript "
                                       f"{used_transcript_accession} (known symbols='{', '.join(alias_symbol_strs)}')"))

    # TODO: alter initial_score based on warning messages of alt not matching?
    # also - _lrg_get_variant_tuple should add matches_reference to search warnings list

    if variant_tuple:
        try:
            results = get_results_from_variant_tuples(variant_qs, variant_tuple)
            variant = results.get()
            if classify:
                yield VariantExtra.classify_variant(
                    variant=variant,
                    genome_build=genome_build
                ), search_messages + reference_message
            else:
                is_genomic = kind == 'g'  # doesn't matter what the preferred genome build is (explicit contig version)
                yield SearchResult(preview=variant.preview, messages=search_messages + reference_message,
                                   ignore_genome_build_mismatch=is_genomic)

        except Variant.DoesNotExist:
            variant_string = Variant.format_tuple(*variant_tuple)
            variant_string_abbreviated = Variant.format_tuple(*variant_tuple, abbreviate=True)
            search_messages.append(SearchMessage(f'"{hgvs_string}" resolved to genomic coordinates "{variant_string_abbreviated}"', LogLevel.INFO, genome_build=genome_build))
            # if we're saying what we're resolving to, no need to complain about reference_message

            # manual variants
            if cmv := VariantExtra.create_manual_variant(for_user=user, genome_build=genome_build, variant_string=variant_string):
                yield SearchResult(cmv, messages=search_messages + reference_message)

            # search for alt alts
            alts = get_results_from_variant_tuples(variant_qs, variant_tuple, any_alt=True)
            for alt in alts:
                alt_messages = search_messages + [SearchMessage(f'No results for alt "{variant_tuple.alt}", but found this using alt "{alt.alt}"', severity=LogLevel.ERROR, substituted=True)]
                yield SearchResult(alt.preview, messages=alt_messages)


DB_PREFIX_PATTERN = re.compile(fr"^(v|{settings.VARIANT_VCF_DB_PREFIX})(\d+)$")


@search_receiver(
    search_type=Variant,
    pattern=DB_PREFIX_PATTERN,
    sub_name="ID",
    example=SearchExample(
        note="Provide the variant ID as it appears in VCF exports",
        examples=[f"{settings.VARIANT_VCF_DB_PREFIX}4638674"]
    )
)
def search_variant_id(search_input: SearchInputInstance):
    yield Variant.objects.filter(pk=search_input.match.group(2))
