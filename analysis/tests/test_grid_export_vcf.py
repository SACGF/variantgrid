import re
import tempfile

import cyvcf2
from django.test import TestCase

from analysis.grid_export import _grid_item_to_vcf_row
from analysis.tests.test_grid_export import GridExportTestCase
from library.genomics.vcf_writer import VCFWriter
from snpdb.models import GenomeBuild
from snpdb.vcf_export_utils import get_vcf_header_from_contigs

# cyvcf2 gt_types: 0=HOM_REF, 1=HET, 3=HOM_ALT (expanded zygosity labels the grid stores)
GT_TYPE_TO_ZYGOSITY = {0: "REF", 1: "HET", 3: "HOM_ALT"}
CONTIG_HEADER_PATTERN = re.compile(r"^##contig=<ID=([^,>]+)")


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
            with open(temp_file.name, "w") as f:
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


class TestGridExportVCFOrder(GridExportTestCase):
    """ (B) The export walks contigs in genome build order, which is the order the header declares
        them in - so records come out sorted consistently with their own header """

    def test_records_follow_header_contig_order(self):
        node = self._sample_node()
        lines = self._export_lines(node, export_type="vcf")

        header_contigs = [CONTIG_HEADER_PATTERN.match(line).group(1)
                          for line in lines if line.startswith("##contig=")]
        self.assertGreater(len(header_contigs), 1)
        contig_order = {name: i for i, name in enumerate(header_contigs)}

        records = [line.split("\t") for line in lines if not line.startswith("#")]
        self.assertEqual(len(records), node.count)
        sort_keys = [(contig_order[r[0]], int(r[1])) for r in records]
        self.assertEqual(sort_keys, sorted(sort_keys))
        # more than one contig actually appears, so the ordering is exercised
        self.assertGreater(len({k[0] for k in sort_keys}), 1)
