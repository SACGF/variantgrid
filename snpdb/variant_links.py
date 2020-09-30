from typing import Dict, Any

from snpdb.models import Variant


def variant_link_info(variant: Variant) -> Dict[str, Any]:
    link_data: Dict[str, Any] = dict()
    coordinate = variant.coordinate
    #FIXME this really needs to refer to SpecialEKeys but can't due to its package
    link_data['variant_coordinate'] = f'{coordinate.chrom}:{coordinate.pos} {coordinate.ref}>{coordinate.alt}'
    link_data['variant_string'] = str(variant)
    link_data['genome_build'] = variant.genome_build.name
    link_data['canonical_c_hgvs'] = variant.canonical_c_hgvs

    # If we can include these in link_data we'll get more links
    #
    # SpecialEKeys.C_HGVS,
    # SpecialEKeys.VARIANT_COORDINATE,
    # SpecialEKeys.GENE_SYMBOL,
    # SpecialEKeys.REFSEQ_TRANSCRIPT_ID,
    # SpecialEKeys.P_HGVS,
    #
    # SpecialEKeys.CLINGEN_ALLELE_ID,
    # SpecialEKeys.GENOME_BUILD,
    # SpecialEKeys.UNIPROT_ID,
    # SpecialEKeys.CLINVAR_VARIANTION_ID,
    # SpecialEKeys.HGNC_ID,
    # SpecialEKeys.GENE_OMIM_ID,
    # SpecialEKeys.UNIPROT_ID

    return link_data
