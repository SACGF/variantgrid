"""
Converting Ensembl's SpliceVault file to MANE RefSeq ids for RefSeq VEP runs (#2426).
"""
import gzip
import os
import tempfile

import pysam
from django.test import TestCase

from annotation.splicevault import (
    convert_splicevault_file,
    default_output_filename,
    mane_ensembl_to_refseq,
)
from genes.enums import AnnotationConsortium, MANEStatus
from genes.fake_data import (
    create_fake_gene_version,
    create_fake_transcript_version,
    insert_transcript_data,
)
from genes.models.models_gene import MANE
from snpdb.models.models_genome import GenomeBuild

HEADER = ("#chrom\tpos\tref\talt\ttranscript_id\tsplicevault_site_type\tsplicevault_site_loc\tspliceai_delta\t"
          "splicevault_out_of_frame_events\tsplicevault_top4_events\tsplicevault_site_sample_count\t"
          "splicevault_site_max_depth\n")
ROWS = [
    "chr1\t925770\tG\tGT\tENST00000342066\tDonor_loss\tchr1:925801\t0.36\tFrameshift:2/2\t"
    "Top1:CD;-35;31%;Frameshift|Top2:CD;+48;0.006%;Frameshift\t16054\t17\n",
    "chr1\t925772\tC\tG\tENST00000999999\tDonor_loss\tchr1:925801\t0.26\tFrameshift:2/2\tTop1:CD;-35;31%;Frameshift\t16054\t17\n",
    "chr2\t100\tC\tG\tENST00000342066\tAcceptor_loss\tchr2:98\t0.5\tFrameshift:0/1\tTop1:ES;2;10%;Inframe\t20\t3\n",
]


class SpliceVaultConversionTests(TestCase):

    def test_ensembl_ids_rewritten_unmatched_dropped_and_indexed(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_filename = os.path.join(tmp_dir, "SpliceVault_data_GRCh38.tsv.gz")
            with gzip.open(input_filename, "wt") as f:
                f.write(HEADER)
                f.writelines(ROWS)
            output_filename = default_output_filename(input_filename)
            self.assertEqual(os.path.basename(output_filename), "SpliceVault_data_GRCh38.refseq_mane.tsv.gz")

            conversion = convert_splicevault_file(input_filename, output_filename,
                                                  {"ENST00000342066": "NM_001385641.1"})
            self.assertEqual((conversion.rows_read, conversion.rows_written, conversion.rows_dropped), (3, 2, 1))
            self.assertEqual((conversion.transcripts_converted, conversion.transcripts_dropped), (1, 1))

            # What the plugin does - a tabix query by position - finds the rewritten row
            with pysam.TabixFile(output_filename) as tabix:
                records = [r.split("\t") for r in tabix.fetch("chr1", 925760, 925780)]
                self.assertEqual(tabix.header[0], HEADER.rstrip("\n"))
        self.assertEqual([(r[1], r[4]) for r in records], [("925770", "NM_001385641.1")])
        self.assertEqual(records[0][9], "Top1:CD;-35;31%;Frameshift|Top2:CD;+48;0.006%;Frameshift")

    def test_mane_map_keys_unversioned_enst_to_versioned_nm(self):
        """ SpliceVault's ENST has no version, a RefSeq VEP cache stable_id does """
        grch38 = GenomeBuild.grch38()
        ensembl_tv = create_fake_transcript_version(grch38)  # ENST00000300305.7
        refseq_gene_version = create_fake_gene_version(grch38, "861", "RUNX1", AnnotationConsortium.REFSEQ)
        refseq_data = {**ensembl_tv.data, "id": "NM_001754.5"}
        refseq_tv = insert_transcript_data(grch38, refseq_data, refseq_gene_version)
        MANE.objects.create(symbol_id="RUNX1", ensembl_transcript_version=ensembl_tv,
                            refseq_transcript_version=refseq_tv, status=MANEStatus.MANE_SELECT)
        self.assertEqual(mane_ensembl_to_refseq(), {"ENST00000300305": "NM_001754.5"})
