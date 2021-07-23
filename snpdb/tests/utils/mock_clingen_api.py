from snpdb.clingen_allele import ClinGenAlleleRegistryAPI, ClinGenAlleleServerException


class MockServerErrorClinGenAlleleRegistryAPI(ClinGenAlleleRegistryAPI):
    """ This always throws ClinGenAlleleServerException """
    def _put(self, url, data):
        raise ClinGenAlleleServerException("PUT", 502, "Bad Gateway")

    @classmethod
    def get(cls, url):
        raise ClinGenAlleleServerException("GET", 502, "Bad Gateway")


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

    def _put(self, url, data):
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
