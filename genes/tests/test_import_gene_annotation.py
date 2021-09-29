from unittest import TestCase

from genes.management.commands.import_gene_annotation2 import convert_transcript_pyreference_to_pyhgvs


class TestAnnotationVCF(TestCase):
    """ These are all from GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.gz """

    def test_convert_transcript_pyreference_to_pyhgvs_positive_strand(self):
        PYREFERENCE_DATA_NM_030751_6 = {
            'biotype': ['protein_coding'],
            'cds_end': 31816192,
            'cds_start': 31608163,
            'chr': 'NC_000010.10',
            'features_by_type': {'3PUTR': [{'start': 31816192, 'stop': 31818732}],
                                 '5PUTR': [{'start': 31608144, 'stop': 31608163}],
                                 'CDS': [{'start': 31608163, 'stop': 31608221},
                                         {'start': 31749965, 'stop': 31750166},
                                         {'start': 31784707, 'stop': 31784767},
                                         {'start': 31791275, 'stop': 31791437},
                                         {'start': 31799600, 'stop': 31799803},
                                         {'start': 31803530, 'stop': 31803636},
                                         {'start': 31809053, 'stop': 31810864},
                                         {'start': 31812860, 'stop': 31813041},
                                         {'start': 31815599, 'stop': 31816192}],
                                 'exon': [{'start': 31608144, 'stop': 31608221},
                                          {'start': 31749965, 'stop': 31750166},
                                          {'start': 31784707, 'stop': 31784767},
                                          {'start': 31791275, 'stop': 31791437},
                                          {'start': 31799600, 'stop': 31799803},
                                          {'start': 31803530, 'stop': 31803636},
                                          {'start': 31809053, 'stop': 31810864},
                                          {'start': 31812860, 'stop': 31813041},
                                          {'start': 31815599, 'stop': 31818732}]},
            'is_coding': 1,
            'start': 31608144,
            'stop': 31818732,
            'strand': '+'
        }

        EXPECTED_DATA = {
            "end": 31818732,
            "chrom": "NC_000010.10",
            "exons": [[31608144, 31608221], [31749965, 31750166], [31784707, 31784767], [31791275, 31791437],
                      [31799600, 31799803], [31803530, 31803636], [31809053, 31810864], [31812860, 31813041],
                      [31815599, 31818732]],
            "start": 31608144,  # We originally had it as 31608145 but this matches exon
            "strand": "+",
            "cds_end": 31816192,
            "cds_start": 31608163,
        }

        output = convert_transcript_pyreference_to_pyhgvs(PYREFERENCE_DATA_NM_030751_6)
        self.maxDiff = None
        self.assertEquals(EXPECTED_DATA, output)

    def test_convert_transcript_pyreference_to_pyhgvs_negative_strand(self):
        PYREFERENCE_DATA_NM_004656_4 = {
            'biotype': ['protein_coding'],
            'cds_end': 52443894,
            'cds_start': 52436303,
            'chr': 'NC_000003.11',
            'features_by_type': {'3PUTR': [{'start': 52435023, 'stop': 52436303}],
                                 '5PUTR': [{'start': 52443894, 'stop': 52444024}],
                                 'CDS': [{'start': 52443857, 'stop': 52443894},
                                         {'start': 52443729, 'stop': 52443759},
                                         {'start': 52443569, 'stop': 52443624},
                                         {'start': 52442489, 'stop': 52442622},
                                         {'start': 52441973, 'stop': 52442093},
                                         {'start': 52441414, 'stop': 52441476},
                                         {'start': 52441189, 'stop': 52441332},
                                         {'start': 52440844, 'stop': 52440923},
                                         {'start': 52440268, 'stop': 52440392},
                                         {'start': 52439780, 'stop': 52439928},
                                         {'start': 52439125, 'stop': 52439310},
                                         {'start': 52438468, 'stop': 52438602},
                                         {'start': 52437431, 'stop': 52437910},
                                         {'start': 52437153, 'stop': 52437314},
                                         {'start': 52436794, 'stop': 52436887},
                                         {'start': 52436617, 'stop': 52436690},
                                         {'start': 52436303, 'stop': 52436437}],
                                 'exon': [{'start': 52443857, 'stop': 52444024},
                                          {'start': 52443729, 'stop': 52443759},
                                          {'start': 52443569, 'stop': 52443624},
                                          {'start': 52442489, 'stop': 52442622},
                                          {'start': 52441973, 'stop': 52442093},
                                          {'start': 52441414, 'stop': 52441476},
                                          {'start': 52441189, 'stop': 52441332},
                                          {'start': 52440844, 'stop': 52440923},
                                          {'start': 52440268, 'stop': 52440392},
                                          {'start': 52439780, 'stop': 52439928},
                                          {'start': 52439125, 'stop': 52439310},
                                          {'start': 52438468, 'stop': 52438602},
                                          {'start': 52437431, 'stop': 52437910},
                                          {'start': 52437153, 'stop': 52437314},
                                          {'start': 52436794, 'stop': 52436887},
                                          {'start': 52436617, 'stop': 52436690},
                                          {'start': 52435023, 'stop': 52436437}]},
            'is_coding': 1,
            'start': 52435023,
            'stop': 52444024,
            'strand': '-',
        }

        EXPECTED_DATA = {
            "end": 52444024,
            "chrom": "NC_000003.11",
            "exons": [
                [52435023, 52436437],
                [52436617, 52436690],
                [52436794, 52436887],
                [52437153, 52437314],
                [52437431, 52437910],
                [52438468, 52438602],
                [52439125, 52439310],
                [52439780, 52439928],
                [52440268, 52440392],
                [52440844, 52440923],
                [52441189, 52441332],
                [52441414, 52441476],
                [52441973, 52442093],
                [52442489, 52442622],
                [52443569, 52443624],
                [52443729, 52443759],
                [52443857, 52444024],
            ],
            "start": 52435023,
            "strand": "-",
            "cds_end": 52443894,
            "cds_start": 52436303
        }

        output = convert_transcript_pyreference_to_pyhgvs(PYREFERENCE_DATA_NM_004656_4)
        self.assertEquals(EXPECTED_DATA, output)
