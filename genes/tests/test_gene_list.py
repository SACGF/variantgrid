import json

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from genes.models import CustomTextGeneList, GeneList, GeneListCategory, GeneSymbol
from snpdb.enums import ImportStatus


class TestAnnotationVCF(TestCase):
    def test_clone_gene_list(self):
        GENE_SYMBOLS = {"RUNX1", "GATA2"}

        user = User.objects.get_or_create(username='test_user')[0]
        gl = GeneList.objects.create(name="test", user=user)
        for symbol in GENE_SYMBOLS:
            gene_symbol = GeneSymbol.objects.get_or_create(symbol=symbol)[0]
            gl.genelistgenesymbol_set.create(gene_symbol=gene_symbol)

        original_pk = gl.pk
        copy = gl.clone()
        copy_symbols = set(copy.genelistgenesymbol_set.values_list("gene_symbol", flat=True))
        self.assertEqual(GENE_SYMBOLS, copy_symbols)
        self.assertEqual(original_pk, gl.pk)
        self.assertNotEqual(original_pk, copy.pk)

    def test_clone_custom_text_gene_list(self):
        user = User.objects.get_or_create(username='test_user')[0]
        gl = GeneList.objects.create(name="custom", user=user)
        custom_text_gene_list = CustomTextGeneList.objects.create(name="custom", text="RUNX1", gene_list=gl)

        copy = custom_text_gene_list.clone()
        self.assertEqual(gl, custom_text_gene_list.gene_list)
        copy.refresh_from_db()
        self.assertNotIn(copy.gene_list_id, {None, gl.pk})

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


class TestGeneListAPI(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='test_user')[0]
        cls.gene_list = GeneList.objects.create(name="api", user=cls.user, import_status=ImportStatus.SUCCESS)
        GeneSymbol.objects.get_or_create(symbol="RUNX1")

    def setUp(self):
        self.client.force_login(self.user)

    def _modify(self, modifications):
        url = reverse("api_modify_gene_list", kwargs={"pk": self.gene_list.pk})
        return self.client.post(url, data={"modifications": json.dumps(modifications)})

    def test_modify_delete_only(self):
        response = self._modify({"delete": ["RUNX1"]})
        self.assertEqual(200, response.status_code)

    def test_modify_unknown_symbol_is_bad_request(self):
        response = self._modify({"add": ["RUNX1", "NOT_A_GENE_SYMBOL"]})
        self.assertEqual(400, response.status_code)
        self.assertFalse(self.gene_list.genelistgenesymbol_set.exists())

    def test_named_gene_list_with_duplicate_names(self):
        category = GeneListCategory.objects.create(name="dupes")
        first = GeneList.objects.create(name="dupe", category=category, user=self.user,
                                        import_status=ImportStatus.SUCCESS)
        GeneList.objects.create(name="dupe", category=category, user=self.user, import_status=ImportStatus.SUCCESS)
        url = reverse("api_named_gene_list", kwargs={"category__name": "dupes", "name": "dupe"})
        response = self.client.get(url)
        self.assertEqual(200, response.status_code)
        self.assertEqual(first.pk, response.json()["pk"])
