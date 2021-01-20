from django.test.testcases import TestCase
from classification.regexes import db_ref_regexes, DbRegexes

class RegexTests(TestCase):

    def test_simple(self):
        text = 'Hello PMID: #3432 There'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'PubMed:3432')

    def test_comma_sep(self):
        text = 'Hello PMID 1111, 2222, this is something else 3333'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 2)
        self.assertEqual(str(results[0]), 'PubMed:1111')
        self.assertEqual(str(results[1]), 'PubMed:2222')

    def test_default(self):
        text = '12345'
        results = db_ref_regexes.search(text, default_regex=DbRegexes.SNP)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'SNP:12345')

    def test_multiple(self):
        text = 'PMID is great and all, but my id is ClinGen CA123 and OMIM: 555, #444 and this is a number 23432'
        results = db_ref_regexes.search(text, default_regex=DbRegexes.SNP)
        self.assertEqual(len(results), 3)
        self.assertEqual(str(results[0]), 'ClinGen:123')
        self.assertEqual(str(results[1]), 'OMIM:000444')
        self.assertEqual(str(results[2]), 'OMIM:000555')

    def test_clingen(self):
        text = 'Without the word cl-in-gen lets see if CA123456 is caught'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'ClinGen:123456')

    def test_comma_diff(self):
        text = 'MedGen:C3150878,OMIM:613616,Orphanet:ORPHA93600'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 2)
        self.assertEqual(str(results[0]), 'MedGen:C3150878')
        self.assertEqual(str(results[1]), 'OMIM:613616')
        #self.assertEqual(str(results[2]), 'Orphanet:ORPHA93600') #Orphanet not supported

    def test_comma_diff_no_rs(self):
        text = 'rs749899964, 6/249'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'SNP:749899964')

    def test_snomed(self):
        text = 'here is a SNOMEDCT:284196006 and another SNOMED-CT:433453 link'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 2)
        self.assertEqual(str(results[0]), 'SNOMED-CT:433453')
        self.assertEqual(str(results[1]), 'SNOMED-CT:284196006')

    def test_ncbi_bookshelf(self):
        text = 'NCBIBookShelf: NBK1426'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'NCBIBookShelf:NBK1426')

    def test_dividers(self):
        text = 'OMIM# 12345 PMID#23456 HPO:123456'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 3)
        self.assertEqual(str(results[0]), 'HP:0123456')
        self.assertEqual(str(results[1]), 'OMIM:012345')
        self.assertEqual(str(results[2]), 'PubMed:23456')

    def test_urls(self):
        text = "Here is a link (http://foo.com) and another one http://bar. That was the end of a sentance, lastly https://abc.com"
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 3)
        self.assertEqual(str(results[0]), 'HTTP://bar')
        self.assertEqual(str(results[1]), 'HTTP://foo.com')
        self.assertEqual(str(results[2]), 'HTTPS://abc.com')

    def test_html(self):
        text = 'Have a look at OMIM:123456 and maybe MedGen:C3150878 good hey'
        result = db_ref_regexes.link_html(text)
        expected = "Have a look at OMIM:<a href='http://www.omim.org/entry/123456'>123456</a> and maybe MedGen:<a href='https://www.ncbi.nlm.nih.gov/medgen/?term=C3150878'>C3150878</a> good hey"
        self.assertEqual(expected, result)

        text = 'OMIM:123456 and maybe MedGen:C3150878'
        result = db_ref_regexes.link_html(text)
        expected = "OMIM:<a href='http://www.omim.org/entry/123456'>123456</a> and maybe MedGen:<a href='https://www.ncbi.nlm.nih.gov/medgen/?term=C3150878'>C3150878</a>"
        self.assertEqual(expected, result)

        text = 'OMIM:123456,234567'
        result = db_ref_regexes.link_html(text)
        expected = "OMIM:<a href='http://www.omim.org/entry/123456'>123456</a>,<a href='http://www.omim.org/entry/234567'>234567</a>"
        self.assertEqual(expected, result)

        text = 'OMIM 123456'
        result = db_ref_regexes.link_html(text)
        expected = "OMIM <a href='http://www.omim.org/entry/123456'>123456</a>"
        self.assertEqual(expected, result)

        text = 'No references here:'
        result = db_ref_regexes.link_html(text)
        expected = 'No references here:'
        self.assertEqual(expected, result)

    def test_min_length(self):
        text = 'rs1 nope, rs12 nope, rs123 yup'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'SNP:123')

    def test_mondo(self):
        text = 'MONDO:001'
        results = db_ref_regexes.search(text)
        self.assertEqual(len(results), 1)
        self.assertEqual(str(results[0]), 'MONDO:0000001')
