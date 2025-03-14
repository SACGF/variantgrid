from snpdb.clingen_allele import ClinGenAlleleRegistryAPI, ClinGenAlleleServerException


class MockServerErrorClinGenAlleleRegistryAPI(ClinGenAlleleRegistryAPI):
    """ This always throws ClinGenAlleleServerException """

    def _put(self, url, data, chunk_size=None):
        raise ClinGenAlleleServerException(url,"PUT", 502, {"description": "Bad Gateway"})

    @classmethod
    def get(cls, url):
        raise ClinGenAlleleServerException(url,"GET", 502, {"description": "Bad Gateway"})


class MockClinGenAlleleRegistryAPI(ClinGenAlleleRegistryAPI):
    """ Cached responses from ClinGen - for unit testing """
    CANONICAL_ALLELES = {
        "CA10617208": {
            '@id': 'http://reg.genome.network/allele/CA10617208',
            'type': 'nucleotide',
            '@context': 'http://reg.genome.network/schema/allele.jsonld',
            'genomicAlleles': [{'hgvs': ['NC_000003.12:g.128480137A>T',
                                         'CM000665.2:g.128480137A>T'],
                                'chromosome': '3',
                                'coordinates': [{'end': 128480137,
                                                 'start': 128480136,
                                                 'allele': 'T',
                                                 'referenceAllele': 'A'}],
                                'referenceGenome': 'GRCh38',
                                'referenceSequence': 'http://reg.genome.network/refseq/RS000051'},
                               {'hgvs': ['NC_000003.11:g.128198980A>T', 'CM000665.1:g.128198980A>T'],
                                'chromosome': '3',
                                'coordinates': [{'end': 128198980,
                                                 'start': 128198979,
                                                 'allele': 'T',
                                                 'referenceAllele': 'A'}],
                                'referenceGenome': 'GRCh37',
                                'referenceSequence': 'http://reg.genome.network/refseq/RS000027'},
                               {'hgvs': ['NC_000003.10:g.129681670A>T'],
                                'chromosome': '3',
                                'coordinates': [{'end': 129681670,
                                                 'start': 129681669,
                                                 'allele': 'T',
                                                 'referenceAllele': 'A'}],
                                'referenceGenome': 'NCBI36',
                                'referenceSequence': 'http://reg.genome.network/refseq/RS000003'},
                               {'hgvs': ['NG_029334.1:g.18051T>A', 'LRG_295:g.18051T>A'],
                                'coordinates': [{'end': 18051,
                                                 'start': 18050,
                                                 'allele': 'A',
                                                 'referenceAllele': 'T'}],
                                'referenceSequence': 'http://reg.genome.network/refseq/RS004348'}],
            'externalRecords': {'dbSNP': [{'rs': 45463895,
                                           '@id': 'http://www.ncbi.nlm.nih.gov/snp/45463895'}],
                                'gnomAD': [{'id': '3-128198980-A-T',
                                            '@id': 'http://gnomad.broadinstitute.org/variant/3-128198980-A-T',
                                            'variant': '3:128198980 A / T'}],
                                'ClinVarAlleles': [
                                    {'@id': 'http://www.ncbi.nlm.nih.gov/clinvar/?term=292483[alleleid]',
                                     'alleleId': 292483,
                                     'preferredName': 'NM_001145661.2(GATA2):c.*882T>A'}],
                                'ClinVarVariations': [{'@id': 'http://www.ncbi.nlm.nih.gov/clinvar/variation/343116',
                                                       'RCV': ['RCV000368151'],
                                                       'variationId': 343116}],
                                'MyVariantInfo_hg19': [{'id': 'chr3:g.128198980A>T',
                                                        '@id': 'http://myvariant.info/v1/variant/chr3:g.128198980A>T?assembly=hg19'}],
                                'MyVariantInfo_hg38': [{'id': 'chr3:g.128480137A>T',
                                                        '@id': 'http://myvariant.info/v1/variant/chr3:g.128480137A>T?assembly=hg38'}]},
            'transcriptAlleles': [{'MANE': {'protein': {'RefSeq': {'hgvs': 'NP_116027.2:p.='},
                                                        'Ensembl': {'hgvs': 'ENSP00000345681.2:p.='}},
                                            'maneStatus': 'MANE Select',
                                            'nucleotide': {'RefSeq': {'hgvs': 'NM_032638.5:c.*882T>A'},
                                                           'Ensembl': {'hgvs': 'ENST00000341105.7:c.*882T>A'}},
                                            'maneVersion': '0.93'},
                                   'hgvs': ['ENST00000341105.7:c.*882T>A'],
                                   'coordinates': [{'end': 2673,
                                                    'start': 2672,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'proteinEffect': {'hgvs': 'ENSP00000345681.2:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS747961'},
                                  {'hgvs': ['ENST00000341105.6:c.*882T>A'],
                                   'coordinates': [{'end': 2657,
                                                    'start': 2656,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'proteinEffect': {'hgvs': 'ENSP00000345681.2:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS261253'},
                                  {'gene': 'http://reg.genome.network/gene/GN004171',
                                   'hgvs': ['NM_001145661.1:c.*882T>A', 'LRG_295t1:c.*882T>A'],
                                   'geneSymbol': 'GATA2',
                                   'coordinates': [{'end': 2760,
                                                    'start': 2759,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'geneNCBI_id': 2624,
                                   'proteinEffect': {'hgvs': 'NP_001139133.1:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS013724'},
                                  {'gene': 'http://reg.genome.network/gene/GN004171',
                                   'hgvs': ['NM_001145662.1:c.*882T>A'],
                                   'geneSymbol': 'GATA2',
                                   'coordinates': [{'end': 2539,
                                                    'start': 2538,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'geneNCBI_id': 2624,
                                   'proteinEffect': {'hgvs': 'NP_001139134.1:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS013725'},
                                  {'gene': 'http://reg.genome.network/gene/GN004171',
                                   'hgvs': ['NM_032638.4:c.*882T>A', 'LRG_295t2:c.*882T>A'],
                                   'geneSymbol': 'GATA2',
                                   'coordinates': [{'end': 2659,
                                                    'start': 2658,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'geneNCBI_id': 2624,
                                   'proteinEffect': {'hgvs': 'NP_116027.2:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS039235'},
                                  {'gene': 'http://reg.genome.network/gene/GN004171',
                                   'hgvs': ['NM_001145661.2:c.*882T>A'],
                                   'geneSymbol': 'GATA2',
                                   'coordinates': [{'end': 2760,
                                                    'start': 2759,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'geneNCBI_id': 2624,
                                   'proteinEffect': {'hgvs': 'NP_001139133.1:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS677924'},
                                  {'MANE': {'protein': {'RefSeq': {'hgvs': 'NP_116027.2:p.='},
                                                        'Ensembl': {'hgvs': 'ENSP00000345681.2:p.='}},
                                            'maneStatus': 'MANE Select',
                                            'nucleotide': {'RefSeq': {'hgvs': 'NM_032638.5:c.*882T>A'},
                                                           'Ensembl': {'hgvs': 'ENST00000341105.7:c.*882T>A'}},
                                            'maneVersion': '0.93'},
                                   'gene': 'http://reg.genome.network/gene/GN004171',
                                   'hgvs': ['NM_032638.5:c.*882T>A'],
                                   'geneSymbol': 'GATA2',
                                   'coordinates': [{'end': 2673,
                                                    'start': 2672,
                                                    'allele': 'A',
                                                    'referenceAllele': 'T'}],
                                   'geneNCBI_id': 2624,
                                   'proteinEffect': {'hgvs': 'NP_116027.2:p.='},
                                   'referenceSequence': 'http://reg.genome.network/refseq/RS699163'}]
        },
        "CA7019515": {"@id": "http://reg.genome.network/allele/CA7019515",
             "type": "nucleotide",
             "@context": "http://reg.genome.network/schema/allele.jsonld",
             "genomicAlleles": [{"hgvs": ["NC_000013.11:g.95186748C>T", "CM000675.2:g.95186748C>T"],
                                 "chromosome": "13",
                                 "coordinates": [{"end": 95186748,
                                                  "start": 95186747,
                                                  "allele": "T",
                                                  "referenceAllele": "C"}],
                                 "referenceGenome": "GRCh38",
                                 "referenceSequence": "http://reg.genome.network/refseq/RS000061"},
                                {"hgvs": ["NC_000013.10:g.95839002C>T", "CM000675.1:g.95839002C>T"],
                                 "chromosome": "13",
                                 "coordinates": [{"end": 95839002,
                                                  "start": 95839001,
                                                  "allele": "T",
                                                  "referenceAllele": "C"}],
                                 "referenceGenome": "GRCh37",
                                 "referenceSequence": "http://reg.genome.network/refseq/RS000037"},
                                {"hgvs": ["NC_000013.9:g.94637003C>T"],
                                 "chromosome": "13",
                                 "coordinates": [{"end": 94637003,
                                                  "start": 94637002,
                                                  "allele": "T",
                                                  "referenceAllele": "C"}],
                                 "referenceGenome": "NCBI36",
                                 "referenceSequence": "http://reg.genome.network/refseq/RS000013"},
                                {"hgvs": ["NG_050651.1:g.119699G>A"],
                                 "coordinates": [{"end": 119699,
                                                  "start": 119698,
                                                  "allele": "A",
                                                  "referenceAllele": "G"}],
                                 "referenceSequence": "http://reg.genome.network/refseq/RS616767"}],
             "externalRecords": {"ExAC": [{"id": "13-95839002-C-T",
                                           "@id": "http://exac.broadinstitute.org/variant/13-95839002-C-T",
                                           "variant": "13:95839002 C / T"}],
                                 "dbSNP": [{"rs": 145886106,
                                            "@id": "http://www.ncbi.nlm.nih.gov/snp/145886106"}],
                                 "gnomAD": [{"id": "13-95839002-C-T",
                                             "@id": "http://gnomad.broadinstitute.org/variant/13-95839002-C-T",
                                             "variant": "13:95839002 C / T"}],
                                 "MyVariantInfo_hg19": [{"id": "chr13:g.95839002C>T",
                                                         "@id": "http://myvariant.info/v1/variant/chr13:g.95839002C>T?assembly=hg19"}],
                                 "MyVariantInfo_hg38": [{"id": "chr13:g.95186748C>T",
                                                         "@id": "http://myvariant.info/v1/variant/chr13:g.95186748C>T?assembly=hg38"}]},
             "transcriptAlleles": [{"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_001105515.2:c.1498G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1630,
                                                     "start": 1629,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_001098985.1:p.Glu500Lys",
                                                      "hgvsWellDefined": "NP_001098985.1:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS011085"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_001301829.1:c.1498G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1630,
                                                     "start": 1629,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_001288758.1:p.Glu500Lys",
                                                      "hgvsWellDefined": "NP_001288758.1:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS025247"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_001301830.1:c.1273G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1405,
                                                     "start": 1404,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_001288759.1:p.Glu425Lys",
                                                      "hgvsWellDefined": "NP_001288759.1:p.Glu425Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS025248"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_005845.4:c.1498G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1630,
                                                     "start": 1629,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_005836.2:p.Glu500Lys",
                                                      "hgvsWellDefined": "NP_005836.2:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS030905"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["XM_005254025.2:c.1369G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1690,
                                                     "start": 1689,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "XP_005254082.1:p.Glu457Lys",
                                                      "hgvsWellDefined": "XP_005254082.1:p.Glu457Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS060466"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["XM_006719914.1:c.1408G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1535,
                                                     "start": 1534,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "XP_006719977.1:p.Glu470Lys",
                                                      "hgvsWellDefined": "XP_006719977.1:p.Glu470Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS073389"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["XM_011521047.1:c.949G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 988,
                                                     "start": 987,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "XP_011519349.1:p.Glu317Lys",
                                                      "hgvsWellDefined": "XP_011519349.1:p.Glu317Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS088211"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["XM_017020319.1:c.1369G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1776,
                                                     "start": 1775,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "XP_016875808.1:p.Glu457Lys",
                                                      "hgvsWellDefined": "XP_016875808.1:p.Glu457Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS571820"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["XM_017020320.2:c.1498G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1617,
                                                     "start": 1616,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "XP_016875809.1:p.Glu500Lys",
                                                      "hgvsWellDefined": "XP_016875809.1:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS571821"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["XM_017020322.1:c.1369G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1639,
                                                     "start": 1638,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "XP_016875811.1:p.Glu457Lys",
                                                      "hgvsWellDefined": "XP_016875811.1:p.Glu457Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS571823"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_001105515.3:c.1498G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1635,
                                                     "start": 1634,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_001098985.1:p.Glu500Lys",
                                                      "hgvsWellDefined": "NP_001098985.1:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS676638"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_001301829.2:c.1498G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1635,
                                                     "start": 1634,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_001288758.1:p.Glu500Lys",
                                                      "hgvsWellDefined": "NP_001288758.1:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS683434"},
                                   {"gene": "http://reg.genome.network/gene/GN000055",
                                    "hgvs": ["NM_001301830.2:c.1273G>A"],
                                    "geneSymbol": "ABCC4",
                                    "coordinates": [{"end": 1410,
                                                     "start": 1409,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "geneNCBI_id": 10257,
                                    "proteinEffect": {"hgvs": "NP_001288759.1:p.Glu425Lys",
                                                      "hgvsWellDefined": "NP_001288759.1:p.Glu425Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS683435"},
                                   {"hgvs": ["ENST00000376887.8:c.1498G>A"],
                                    "coordinates": [{"end": 1613,
                                                     "start": 1612,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "proteinEffect": {"hgvs": "ENSP00000366084.4:p.Glu500Lys",
                                                      "hgvsWellDefined": "ENSP00000366084.4:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS271356"},
                                   {"hgvs": ["ENST00000536256.3:c.1273G>A"],
                                    "coordinates": [{"end": 1392,
                                                     "start": 1391,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "proteinEffect": {"hgvs": "ENSP00000442024.1:p.Glu425Lys",
                                                      "hgvsWellDefined": "ENSP00000442024.1:p.Glu425Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS360179"},
                                   {"hgvs": ["ENST00000629385.1:c.1498G>A"],
                                    "coordinates": [{"end": 1630,
                                                     "start": 1629,
                                                     "allele": "A",
                                                     "referenceAllele": "G"}],
                                    "proteinEffect": {"hgvs": "ENSP00000487081.1:p.Glu500Lys",
                                                      "hgvsWellDefined": "ENSP00000487081.1:p.Glu500Lys"},
                                    "referenceSequence": "http://reg.genome.network/refseq/RS407216"}]},
        "CA337095804": {
            "@context": "http://reg.genome.network/schema/allele.jsonld",
            "@id": "http://reg.genome.network/allele/CA337095804",
            "communityStandardTitle": [
                "NC_012920.1:m.263A>G"
            ],
            "externalRecords": {
                "ClinVarAlleles": [
                    {
                        "@id": "http://www.ncbi.nlm.nih.gov/clinvar/?term=434777[alleleid]",
                        "alleleId": 434777,
                        "preferredName": "NC_012920.1(MT-CYB):m.263A>G"
                    }
                ],
                "ClinVarVariations": [
                    {
                        "@id": "http://www.ncbi.nlm.nih.gov/clinvar/variation/441147",
                        "RCV": [
                            "RCV000509228"
                        ],
                        "variationId": 441147
                    }
                ],
                "MyVariantInfo_hg38": [
                    {
                        "@id": "http://myvariant.info/v1/variant/chrMT:g.263A>G?assembly=hg38",
                        "id": "chrMT:g.263A>G"
                    }
                ],
                "dbSNP": [
                    {
                        "@id": "http://www.ncbi.nlm.nih.gov/snp/2853515",
                        "rs": 2853515
                    }
                ]
            },
            "genomicAlleles": [
                {
                    "chromosome": "MT",
                    "coordinates": [
                        {
                            "allele": "G",
                            "end": 263,
                            "referenceAllele": "A",
                            "start": 262
                        }
                    ],
                    "hgvs": [
                        "NC_012920.1:m.263A>G",
                        "J01415.2:m.263A>G"
                    ],
                    "referenceGenome": "GRCh38",
                    "referenceSequence": "http://reg.genome.network/refseq/RS000433"
                }
            ],
            "type": "nucleotide"
        }
    }

    CACHED_HGVS_RESPONSES = {
        "NC_000001.10:g.169519049T>C": {
            'region': '[169519048,169519049)',
            'message': 'Cannot align NC_000001.10 [169519048,169519049). ',
            'errorType': 'NoConsistentAlignment',
            'inputLine': 'NC_000001.10:g.169519049T>C',
            'description': 'Given allele cannot be mapped in consistent way to reference genome.',
            'sequenceName': 'NC_000001.10'
        },
        "NC_000003.11:g.128198980A>T": CANONICAL_ALLELES["CA10617208"],  # 37
        "NC_000003.12:g.128480137A>T": CANONICAL_ALLELES["CA10617208"],  # 38
    }

    def _put(self, url, data, chunk_size=None):
        clingen_records = []
        for line in data.split("\n"):
            if response := self.CACHED_HGVS_RESPONSES.get(line):
                clingen_records.append(response)
            else:
                raise ValueError(f"MockClinGenAlleleRegistryAPI: unknown HGVS: '{line}'")
        return clingen_records

    @classmethod
    def get_code(cls, code):
        return cls.CANONICAL_ALLELES[code]
