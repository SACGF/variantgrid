"""

@see http://reg.clinicalgenome.org/doc/AlleleRegistry_1.01.xx_api_v1.pdf

A variant with coordinates ('14', 45628283, 'TTA', 'T') looks like this in
http://reg.test.genome.network/allele/CA7169043

api_response["genomicAlleles"] coordinates look like:

{"allele": "",
 "end": 45628285,
 "referenceAllele": "TA",
 "start": 45628283}

So to get out Variant coords it's easier to just use HGVS as we'd have to look up the
reference base anyway to construct it

"""
import hashlib
import itertools
import logging
import time
from typing import Dict, List, Tuple, Optional

import requests
from django.conf import settings

from genes.hgvs import HGVSMatcher
from library.django_utils import thread_safe_unique_together_get_or_create
from library.utils import iter_fixed_chunks, get_single_element
from snpdb.models import Allele, ClinGenAllele, GenomeBuild, Variant, VariantAllele, Contig, GenomeFasta
from snpdb.models.models_enums import AlleleOrigin, AlleleConversionTool, ClinGenAlleleExternalRecordType
from upload.models import ModifiedImportedVariant


class ClinGenAlleleAPIException(Exception):

    @staticmethod
    def from_api_response(api_response):
        error_type = api_response['errorType']
        description = api_response['description']
        input_line = api_response['inputLine']
        message = f"ClinGeneAllele API Error: {error_type} ({description}) for input '{input_line}'"
        return ClinGenAlleleAPIException(message)


# throw exceptions in case of errors
def clingen_allele_registry_put(url, data, login, password):
    logging.debug("Calling ClinGen API")
    # calculate token & full URL
    identity = hashlib.sha1((login + password).encode('utf-8')).hexdigest()
    gbTime = str(int(time.time()))
    token = hashlib.sha1((url + identity + gbTime).encode('utf-8')).hexdigest()
    request = url + '&gbLogin=' + login + '&gbTime=' + gbTime + '&gbToken=' + token
    # send request & parse response
    res = requests.put(request, data=data)
    if res.status_code != 200:
        raise Exception(f"Error for PUT requests: {res}")
    return res.json()


def clingen_allele_get_id_url(code):
    return settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + f"/allele/{code}"


def clingen_allele_url_get(url):
    res = requests.get(url)
    if res.status_code != 200:
        raise Exception(f"Error for GET requests: {res}")
    return res.json()


def clingen_hgvs_put_iter(hgvs_iter):
    """ Calls ClinGen in batches """
    url = settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + "/alleles?file=hgvs"
    login = settings.CLINGEN_ALLELE_REGISTRY_LOGIN
    password = settings.CLINGEN_ALLELE_REGISTRY_PASSWORD
    chunk_size = settings.CLINGEN_ALLELE_REGISTRY_MAX_RECORDS

    for hgvs_chunk in iter_fixed_chunks(hgvs_iter, chunk_size):
        data = "\n".join(hgvs_chunk)
        yield clingen_allele_registry_put(url, data, login, password)


def clingen_hgvs_put(hgvs_iter):
    return itertools.chain.from_iterable(clingen_hgvs_put_iter(hgvs_iter))


def populate_clingen_alleles_for_variants(genome_build: GenomeBuild, variants):
    """ Ensures that we have ClinGenAllele for variants
        You don't need to, but more efficient to make variants distinct """

    qs = VariantAllele.objects.filter(variant__in=variants).values_list("variant_id", flat=True)
    variant_ids_with_allele = set(qs)

    reference_variant_ids_without_alleles = []
    # These are kept in sync
    variant_ids_without_alleles = []
    variant_hgvs = []

    hgvs_matcher = HGVSMatcher(genome_build)
    num_existing_records = 0
    miv_qs = ModifiedImportedVariant.objects.filter(variant__in=variants)
    normalized_variants = set(miv_qs.values_list("variant_id", flat=True))

    for v in variants:
        variant_id = v.pk
        if variant_id not in variant_ids_with_allele:
            if v.is_reference:
                reference_variant_ids_without_alleles.append(variant_id)
            else:
                variant_ids_without_alleles.append(variant_id)
                variant_hgvs.append(hgvs_matcher.variant_to_g_hgvs(v))
        else:
            num_existing_records += 1

    num_no_record = len(variant_ids_without_alleles)
    logging.debug("ClinGeneAllele %s: %d variants have Alleles, %d without", genome_build, num_existing_records, num_no_record)

    variant_id_allele_error: List[Tuple[int, Allele, Optional[str]]] = []
    if reference_variant_ids_without_alleles:
        logging.debug("%d reference variants", len(reference_variant_ids_without_alleles))
        empty_alleles = [Allele()] * len(reference_variant_ids_without_alleles)
        reference_alleles = Allele.objects.bulk_create(empty_alleles)
        variant_id_allele_error.extend(((variant_id, allele, None) for variant_id, allele in
                                        zip(reference_variant_ids_without_alleles, reference_alleles)))

    if variant_hgvs:
        # Create Allele / ClinGenAlleles if they don't exist
        # Note: There may already exist an Allele (created for another build but not linked to this build's Variant)
        # So we have to create (ignoring errors from dupes) then pull out again.
        clingen_allele_list = []

        # bulk_create doesn't return populated objects when ignore_conflicts=True
        # so break up empty (no chance of failure) and alleles_w_clingen
        allele_no_clingen_list = []
        alleles_with_clingen_list = []
        clingen_by_variant_id = {}
        clingen_errors_by_variant_id = {}

        clingen_response = clingen_hgvs_put(variant_hgvs)
        for variant_id, api_response in zip(variant_ids_without_alleles, clingen_response):
            if "errorType" in api_response:
                clingen_errors_by_variant_id[variant_id] = api_response
                allele_no_clingen_list.append(Allele())
            else:
                # Don't store errors, retry in case we screwed up on our end
                clingen_allele_id = ClinGenAllele.get_id_from_response(api_response)
                clingen_by_variant_id[variant_id] = clingen_allele_id
                cga = ClinGenAllele(id=clingen_allele_id, api_response=api_response)
                clingen_allele_list.append(cga)
                allele = Allele(clingen_allele_id=clingen_allele_id)
                alleles_with_clingen_list.append(allele)

        if clingen_allele_list:
            ClinGenAllele.objects.bulk_create(clingen_allele_list, ignore_conflicts=True)

        allele_no_clingen_list = Allele.objects.bulk_create(allele_no_clingen_list)
        alleles_with_clingen_list = Allele.objects.bulk_create(alleles_with_clingen_list, ignore_conflicts=True)
        alleles_by_clingen = {}
        existing_allele_clingen_ids = [a.clingen_allele_id for a in alleles_with_clingen_list if a.pk is None]
        if existing_allele_clingen_ids:
            allele_qs = Allele.objects.filter(clingen_allele_id__in=existing_allele_clingen_ids)
            alleles_by_clingen = {a.clingen_allele_id: a for a in allele_qs}

        for variant_id in variant_ids_without_alleles:
            error = None
            clingen_allele_id = clingen_by_variant_id.get(variant_id)
            if clingen_allele_id:
                allele = alleles_by_clingen[clingen_allele_id]
            else:
                error = clingen_errors_by_variant_id[variant_id]
                allele = allele_no_clingen_list.pop()  # Any empty Allele will do
            variant_id_allele_error.append((variant_id, allele, error))

    if variant_id_allele_error:
        variant_allele_list = []
        for (variant_id, allele, error) in variant_id_allele_error:
            if variant_id in normalized_variants:
                origin = AlleleOrigin.IMPORTED_NORMALIZED
            else:
                origin = AlleleOrigin.IMPORTED_TO_DATABASE

            va = VariantAllele(variant_id=variant_id,
                               genome_build=genome_build,
                               allele=allele,
                               origin=origin,
                               error=error)
            variant_allele_list.append(va)

        logging.debug("Creating %d VariantAlleles", len(variant_allele_list))
        VariantAllele.objects.bulk_create(variant_allele_list, ignore_conflicts=True)


def get_variant_allele_for_variant(genome_build: GenomeBuild, variant: Variant) -> VariantAllele:
    """ Retrieves from DB or calls API then caches in DB
        Successful calls link variants in all builds (that exist)
        errors are only stored on the requesting build """

    if not variant.is_standard_variant:
        msg = f"No ClinGenAllele for non-standard variant: {variant}"
        raise ClinGenAlleleAPIException(msg)

    try:
        va = VariantAllele.objects.get(variant=variant)
        if va.needs_clinvar_call():
            va = variant_allele_clingen(genome_build, variant, existing_variant_allele=va)

    except VariantAllele.DoesNotExist:
        if settings.CLINGEN_ALLELE_REGISTRY_LOGIN:
            va = variant_allele_clingen(genome_build, variant)
        else:
            logging.debug("Not using ClinGen")
            allele = Allele.objects.create()
            va = VariantAllele.objects.create(variant_id=variant.pk,
                                              genome_build=genome_build,
                                              allele=allele,
                                              origin=AlleleOrigin.variant_origin(variant))
    return va


def variant_allele_clingen(genome_build, variant, existing_variant_allele=None) -> VariantAllele:
    """ Call ClinGen and setup VariantAllele - use existing if provided, otherwise create """
    g_hgvs = HGVSMatcher(genome_build).variant_to_g_hgvs(variant)
    ca_rep_size = len(g_hgvs)
    max_size = ClinGenAllele.CLINGEN_ALLELE_MAX_REPRESENTATION_SIZE
    if ca_rep_size > max_size:
        msg = f"Representation has size {ca_rep_size} which exceeds ClinGen max of {max_size}"
        raise ClinGenAlleleAPIException(msg)

    api_response = get_single_element(list(clingen_hgvs_put([g_hgvs])))
    va: VariantAllele
    if "errorType" in api_response:
        if existing_variant_allele:
            existing_variant_allele.error = api_response
            existing_variant_allele.save()
            va = existing_variant_allele
        else:
            allele = Allele.objects.create()
            va = VariantAllele.objects.create(variant_id=variant.pk,
                                              genome_build=genome_build,
                                              allele=allele,
                                              origin=AlleleOrigin.variant_origin(variant),
                                              error=api_response)

    else:
        clingen_allele_id = ClinGenAllele.get_id_from_response(api_response)
        kwargs = {"pk": clingen_allele_id,
                  "defaults": {"api_response": api_response}}
        clingen_allele, _ = thread_safe_unique_together_get_or_create(ClinGenAllele, **kwargs)
        if existing_variant_allele:
            # Link existing_variant_allele.allele to ClinGen (merging if already exists)
            try:
                logging.debug("ClinGen (%s) exists - merging existing %s!", clingen_allele, existing_variant_allele)
                existing_allele_with_clingen = Allele.objects.get(clingen_allele=clingen_allele)
                existing_allele_with_clingen.merge(existing_variant_allele.allele)
                existing_variant_allele.allele = existing_variant_allele
                existing_variant_allele.save()
            except Allele.DoesNotExist:
                logging.debug("Linking existing %s against clingen %s", existing_variant_allele, clingen_allele)
                existing_variant_allele.allele.clingen_allele = clingen_allele
                existing_variant_allele.allele.save()
            allele = existing_variant_allele.allele
        else:
            allele, _ = Allele.objects.get_or_create(clingen_allele=clingen_allele)

        known_variants = {variant.genome_build: variant}  # Ensure this gets linked to Allele, no matter API response
        variant_allele_by_build = link_allele_to_existing_variants(allele,
                                                                   AlleleConversionTool.CLINGEN_ALLELE_REGISTRY,
                                                                   known_variants=known_variants)
        va = variant_allele_by_build[genome_build]
    return va


def get_clingen_allele_for_variant(genome_build: GenomeBuild, variant: Variant) -> ClinGenAllele:
    """ Retrieves from DB or calls API then caches in DB   """
    va = get_variant_allele_for_variant(genome_build, variant)

    if va.allele.clingen_allele is None:
        if va.error:
            raise ClinGenAlleleAPIException.from_api_response(va.error)
        else:
            if not settings.CLINGEN_ALLELE_REGISTRY_LOGIN:
                raise ClinGenAlleleAPIException("'settings.CLINGEN_ALLELE_REGISTRY_LOGIN' not set")
            ClinGenAlleleAPIException(f"Allele {va.allele} has no clingen_allele or error")
    return va.allele.clingen_allele


def get_clingen_allele(code: str) -> ClinGenAllele:
    """ Retrieves from DB or calls API then caches in DB   """

    clingen_allele_id = ClinGenAllele.get_id_from_code(code)
    try:
        clingen_allele = ClinGenAllele.objects.get(pk=clingen_allele_id)
    except ClinGenAllele.DoesNotExist:
        url = clingen_allele_get_id_url(code)
        api_response = clingen_allele_url_get(url)
        if "errorType" in api_response:
            raise ClinGenAlleleAPIException.from_api_response(api_response)
        clingen_allele = ClinGenAllele.objects.create(pk=clingen_allele_id, api_response=api_response)

    allele, _ = Allele.objects.get_or_create(clingen_allele=clingen_allele)
    link_allele_to_existing_variants(allele, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)

    return clingen_allele


def get_clingen_alleles_from_external_code(er_type: ClinGenAlleleExternalRecordType, external_code) -> List[ClinGenAllele]:
    # TODO: We could cache this? At the moment we have to make a new API call every build

    suffix = f"/alleles?{er_type.value}={external_code}"
    url = settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + suffix

    # Returns a list for external records
    api_response_list = clingen_allele_url_get(url)
    clingen_allele_list = []
    for api_response in api_response_list:
        clingen_allele_id = ClinGenAllele.get_id_from_response(api_response)
        clingen_allele, _ = ClinGenAllele.objects.get_or_create(pk=clingen_allele_id,
                                                                defaults={"api_response": api_response})
        clingen_allele_list.append(clingen_allele)
        allele, _ = Allele.objects.get_or_create(clingen_allele=clingen_allele)
        link_allele_to_existing_variants(allele, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)

    return clingen_allele_list


def link_allele_to_existing_variants(allele: Allele, conversion_tool,
                                     known_variants=None) -> Dict[GenomeBuild, VariantAllele]:
    """ known_variants - be able to pass in variants you already know are linked to Allele. Workaround to deal with
        ClinGen Allele registry accepting a coordinate (eg NC_000001.10:g.145299792A>G GRCh37) and retrieving record
        CA1051218 but the API response won't have GRCh37 in it! """
    if known_variants is None:
        known_variants = {}

    variant_allele_by_build = {}
    # Attempt to link w/any existing Variant
    for genome_build in GenomeBuild.builds_with_annotation():
        try:
            try:
                variant_tuple = allele.clingen_allele.get_variant_tuple(genome_build)
                variant = Variant.get_from_tuple(variant_tuple, genome_build)
            except ClinGenAllele.ClinGenBuildNotInResponseError:
                variant = known_variants.get(genome_build)
                if variant is None:
                    raise

            defaults = {"genome_build": genome_build,
                        "allele": allele,
                        "origin": AlleleOrigin.variant_origin(variant),
                        "conversion_tool": conversion_tool}
            va, _ = VariantAllele.objects.get_or_create(variant=variant,
                                                        defaults=defaults)
            variant_allele_by_build[genome_build] = va
        except Variant.DoesNotExist:
            pass  # Variant may not be created for build
        except (Contig.ContigNotInBuildError, GenomeFasta.ContigNotInFastaError, ClinGenAllele.ClinGenBuildNotInResponseError) as e:
            logging.error(e)

    return variant_allele_by_build
