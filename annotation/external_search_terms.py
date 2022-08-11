""" A list of (Google) search terms for a variant """
from collections import defaultdict
from typing import Tuple, List, Optional

from pyhgvs import HGVSName, InvalidHGVSName


def _get_gene_and_terms(vta, c_hgvs=True, extra_terms: Optional[List[str]] = None) -> Tuple[str, List]:
    from annotation.models import VariantTranscriptAnnotation

    if vta.gene:
        gene_symbol = str(vta.gene.get_gene_symbol(vta.version.genome_build))
    else:
        gene_symbol = None

    terms = []
    if extra_terms:
        terms.extend(extra_terms)
    hgvs_name = None
    try:
        hgvs_name = HGVSName(vta.hgvs_c)
    except (NotImplementedError, InvalidHGVSName):
        pass

    if c_hgvs and hgvs_name:
        # "7051G>A" | "7051G->A" | "7051G-->A" | "7051G/A"
        if hgvs_name.mutation_type == ">":
            cdna_coords = hgvs_name.format_cdna_coords()
            ref, alt = hgvs_name.get_ref_alt()
            for change_symbol in [">", "->", "-->", "/"]:
                terms.append(f"{cdna_coords}{ref}{change_symbol}{alt}")

    if vta.hgvs_p:
        # pyHGVS doesn't handle p.HGVS very well, so can't use HGVSName.format_protein() etc
        protein_aa3 = vta.hgvs_p.split(":p.")[1]
        terms.append(protein_aa3)
        protein_aa1 = VariantTranscriptAnnotation.amino_acid_3_to_1(protein_aa3)
        terms.append(protein_aa1)

    if hgvs_name and hgvs_name.mutation_type in ('ins', 'del', 'dup'):
        # "del ex 20" and "del exon 20"
        ex_in_terms = [(vta.intron, ["intron", "in"]),
                       (vta.exon, ["exon", "ex"])]
        for val, in_terms in ex_in_terms:
            if val:
                num = val.split("/")[0]  # looks like: "20/24"
                for t in in_terms:
                    terms.append(f"{hgvs_name.mutation_type} {t} {num}")

    return gene_symbol, terms


def _get_search_terms(variant_transcripts_list: List, formatter: str = None, **kwargs):
    gene_terms = defaultdict(set)
    for vta in variant_transcripts_list:
        gene_symbol, terms = _get_gene_and_terms(vta, **kwargs)
        if gene_symbol:
            gene_terms[gene_symbol].update(terms)

    searches = []
    for gene_symbol, terms in gene_terms.items():
        if formatter:
            gene_symbol = formatter % gene_symbol
            terms = [formatter % s for s in terms]

        and_terms = [gene_symbol]
        optional_or = " OR ".join(terms)
        if optional_or:
            and_terms.append(f"({optional_or})")
        search_terms = " AND ".join(and_terms)
        searches.append(search_terms)

    if len(searches) == 1:
        return searches[0]
    return " OR ".join(["(%s)" % s for s in searches])


def get_variant_search_terms(variant_transcripts_list: List, extra_terms: Optional[List[str]] = None):
    """
    Examples:

    BRCA1 AND ("5194-?_5277+?del" OR "del exon 20" OR "His1732_Lys1759del" OR "H1732_K1759del" OR  "del ex 20")
    BRCA1 AND ("5207T>C" OR "5326T>C" OR "Val1736Ala" OR "V1736A")
    BRCA1+[("736T>G")or("855T>G")or("L246V")or("Leu246Val")]
    BRCA2 AND ("10110" OR "10338" OR "Arg3370Arg" OR "R3370R" OR "Arg3370=" OR "Arg3370" OR "R3370")
    BRCA2 AND ("6339") BRCA2:c.6339C>T BRCA2 AND (p.Asn2113")
    "BRCA2" ("7051G>A" | "7051G->A" | "7051G-->A" | "7051G/A" | "Ala2351Thr" | "A2351T" | rs80358930)
    BRCA2 AND ("3137A" OR "3365A" OR "Glu1046Gly" OR "E1046G")
    "PTEN" ("1000A>T" | "1000A->T" | "1000A-->T" | "1000A/T" | "Asn334Tyr" | "N334Y")
    PTEN AND ("1000A" OR "Asn334Tyr" OR "N334Y")

    """

    return _get_search_terms(variant_transcripts_list, formatter='"%s"', extra_terms=extra_terms)


def get_variant_pubmed_search_terms(variant_transcripts_list: List):
    """ Examples:

            (CFTR) AND ((Arg117His) OR (R117H))

        PubMed doesn't like rsIds or c.HGVS """
    return _get_search_terms(variant_transcripts_list, formatter='(%s)', c_hgvs=False)
