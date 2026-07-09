import tempfile

import cyvcf2
from django.test import TestCase

from analysis.grid_export import _grid_item_to_vcf_row
from library.genomics.vcf_writer import VCFWriter
from snpdb.models import GenomeBuild
from snpdb.vcf_export_utils import get_vcf_header_from_contigs

# cyvcf2 gt_types: 0=HOM_REF, 1=HET, 3=HOM_ALT (expanded zygosity labels the grid stores)
GT_TYPE_TO_ZYGOSITY = {0: "REF", 1: "HET", 3: "HOM_ALT"}


class TestGridExportVCF(TestCase):
    """ Round-trip test for the analysis-grid VCF export (Writer 4 in the issue #1068
        consolidation). A real genotyped VCF is mapped into grid items, written out as a complete
        VCF via the grid export path, then re-read to confirm the genotypes survive. """

    def test_roundtrip(self):
        genome_build = GenomeBuild.get_name_or_alias("T2T-CHM13v2.0")
        source = list(cyvcf2.Reader("upload/test_data/vcf/t2t_brca2.vcf"))
        self.assertTrue(source)  # sanity: fixture has records

        sample_ids = [7]
        samples = ["S7"]
        info_dict = {}

        items = []
        for v in source:
            zygosity = GT_TYPE_TO_ZYGOSITY[int(v.gt_types[0])]
            items.append({
                "locus__contig__name": v.CHROM, "locus__position": v.POS,
                "variantannotation__dbsnp_rs_id": v.ID,
                "locus__ref__seq": v.REF, "alt__seq": v.ALT[0],
                "sample_7_samples_zygosity": zygosity,
                # this VCF has no depth/frequency fields - they should render as '.'
                "sample_7_samples_allele_depth": None,
                "sample_7_samples_read_depth": None,
                "sample_7_samples_allele_frequency": None,
            })

        header_lines = get_vcf_header_from_contigs(genome_build, info_dict, samples, use_accession=False)
        with tempfile.NamedTemporaryFile(mode="wt", suffix=".vcf", delete=True) as temp_file:
            with open(temp_file.name, "wt") as f:
                writer = VCFWriter(f, header_lines)
                for item in items:
                    chrom, pos, vcf_id, ref, alt, info, fmt, sample_calls = \
                        _grid_item_to_vcf_row(info_dict, item, sample_ids, samples, use_accession=False)
                    writer.write_record(chrom, pos, ref, alt, vcf_id=vcf_id, info=info,
                                        fmt=fmt, sample_calls=sample_calls)
            written = list(cyvcf2.Reader(temp_file.name))

        self.assertEqual(len(source), len(written))
        for src, out in zip(source, written):
            self.assertEqual(src.CHROM, out.CHROM)
            self.assertEqual(src.POS, out.POS)
            self.assertEqual(list(src.ALT), list(out.ALT))
            # genotype (from expanded zygosity) round-trips back to the source call
            self.assertEqual(src.genotypes[0], out.genotypes[0])
