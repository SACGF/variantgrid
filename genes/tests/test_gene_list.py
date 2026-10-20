from django.contrib.auth.models import User
from django.test import TestCase

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

    def test_copy_for_user(self):
        owner = User.objects.get_or_create(username='test_owner')[0]
        user = User.objects.get_or_create(username='test_user')[0]
        gl = GeneList.objects.create(name="panel", user=owner, locked=True, url="https://example.com/panel")
        gene_symbol = GeneSymbol.objects.get_or_create(symbol="RUNX1")[0]
        gl.genelistgenesymbol_set.create(gene_symbol=gene_symbol, original_name="RUNX1")

        copy = gl.copy_for_user(user)
        self.assertEqual(("Copy of panel", user, False, None), (copy.name, copy.user, copy.locked, copy.url))
        self.assertTrue(copy.can_write(user))
        self.assertFalse(copy.can_view(owner))
        self.assertEqual(["RUNX1"], list(copy.get_gene_names()))
