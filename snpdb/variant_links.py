from typing import Any

from genes.hgvs import HGVSException
from snpdb.models import Variant, GenomeBuild


def variant_link_info(variant: Variant, genome_build: GenomeBuild) -> dict[str, Any]:
    """ Needs to be passed a VariantAllele """
    link_data: dict[str, Any] = {}
    coordinate = variant.coordinate

    #FIXME this really needs to refer to SpecialEKeys but can't due to its package
    link_data['variant_coordinate'] = str(coordinate)
    link_data['variant_string'] = str(variant)
    link_data['variant_svlen'] = variant.svlen
    link_data['genome_build'] = genome_build.name

    if cta := variant.get_canonical_transcript_annotation(genome_build):
        link_data['gene_symbol'] = cta.symbol
        try:
            link_data['canonical_c_hgvs'] = cta.get_hgvs_c_with_symbol()
        except HGVSException:
            pass

    # If we can include these in link_data we'll get more links
    #
    # SpecialEKeys.REFSEQ_TRANSCRIPT_ID,
    # SpecialEKeys.P_HGVS,
    #
    # SpecialEKeys.CLINGEN_ALLELE_ID,
    # SpecialEKeys.UNIPROT_ID,
    # SpecialEKeys.CLINVAR_VARIANTION_ID,
    # SpecialEKeys.HGNC_ID,
    # SpecialEKeys.GENE_OMIM_ID,
    # SpecialEKeys.UNIPROT_ID

    return link_data
