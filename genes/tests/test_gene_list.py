from unittest import TestCase

from django.contrib.auth.models import User

from genes.models import GeneList, GeneSymbol


class TestAnnotationVCF(TestCase):
    def test_clone_gene_list(self):
        GENE_SYMBOLS = {"RUNX1", "GATA2"}

        user = User.objects.get_or_create(username='test_user')[0]
        gl = GeneList.objects.create(name="test", user=user)
        for symbol in GENE_SYMBOLS:
            gene_symbol = GeneSymbol.objects.get_or_create(symbol=symbol)[0]
            gl.genelistgenesymbol_set.create(gene_symbol=gene_symbol)

        copy = gl.clone()
        copy_symbols = set(copy.genelistgenesymbol_set.values_list("gene_symbol", flat=True))
        self.assertEqual(GENE_SYMBOLS, copy_symbols)
