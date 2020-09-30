""" A list of (Google) search terms for a variant """
from pyhgvs import HGVSName

from annotation.models import VariantTranscriptAnnotation


def get_variant_transcript_annotation_search_terms(vta: VariantTranscriptAnnotation):
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

    mandatory_terms = []
    optional_terms = []
    if vta.gene:
        gene_symbol = str(vta.gene.get_gene_symbol(vta.version.genome_build))
        mandatory_terms.append(gene_symbol)

    if vta.hgvs_c:
        try:
            # "7051G>A" | "7051G->A" | "7051G-->A" | "7051G/A"
            c_name = HGVSName(vta.hgvs_c)
            if c_name.mutation_type == ">":
                cdna_coords = c_name.format_cdna_coords()
                ref, alt = c_name.get_ref_alt()
                for change_symbol in [">", "->", "-->", "/"]:
                    optional_terms.append(f"{cdna_coords}{ref}{change_symbol}{alt}")
        except NotImplementedError:
            pass

    if vta.hgvs_p:
        # pyHGVS doesn't handle p.HGVS very well, so can't use HGVSName.format_protein() etc
        protein_aa3 = vta.hgvs_p.split(":p.")[1]
        optional_terms.append(protein_aa3)
        protein_aa1 = VariantTranscriptAnnotation.amino_acid_3_to_1(protein_aa3)
        optional_terms.append(protein_aa1)

    if vta.dbsnp_rs_id:
        optional_terms.append(vta.dbsnp_rs_id)

    if vta.variant.is_deletion:
        # "del ex 20" and "del exon 20"
        ex_in_terms = [(vta.intron, ["intron", "in"]),
                       (vta.exon, ["exon", "ex"])]
        for val, terms in ex_in_terms:
            if val:
                num = val.split("/")[0]  # looks like: "20/24"
                for t in terms:
                    optional_terms.append(f"del {t} {num}")

    # double quote strings
    optional_terms = [f'"{s}"' for s in optional_terms]
    mandatory_terms = [f'"{s}"' for s in mandatory_terms]

    optional_or = " OR ".join(optional_terms)
    mandatory_terms.append(f"({optional_or})")
    search_terms = " AND ".join(mandatory_terms)
    return search_terms
