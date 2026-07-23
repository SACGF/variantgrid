from django.test import TestCase, override_settings
from django.urls import reverse
from rest_framework.test import APIClient

from beacon.models import BeaconInboundQuery
from beacon.schema import OBSERVATIONS_DATASET_ID, CLASSIFICATIONS_DATASET_ID
from classification.enums import SpecialEKeys, SubmissionSource
from classification.enums.classification_enums import ShareLevel
from classification.models.classification import Classification, ClassificationModification
from classification.tests.models.test_utils import ClassificationTestUtils
from library.guardian_utils import public_group
from patients.models_enums import Zygosity
from snpdb.models import GenomeBuild, CohortGenotypeCollection, CohortGenotype
from snpdb.models.models_zygosity_counts import VariantZygosityCountCollection, VariantZygosityCount
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant, create_mock_allele


@override_settings(BEACON_ENABLED=True)
class BeaconGVariantsTestCase(TestCase):
    """ g_variants query: presence/absence, granularity tiers, permission tiering (§11) and
        the two datasets (observations + public classifications, §5.5). """

    def setUp(self):
        public_group()  # anonymous permission checks resolve against the public group
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()

        self.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.cohort = create_fake_cohort(self.user, self.grch37)  # proband readable by user only
        self.cgc = CohortGenotypeCollection.objects.get(cohort=self.cohort)
        self.global_collection, _ = VariantZygosityCountCollection.objects.get_or_create(
            name="global", defaults={"description": "global germline counts"})

        # Variant with a private observation only (proband HET; mother/father missing).
        self.obs_variant = slowly_create_test_variant("3", 1000, "A", "T", self.grch37)
        self._make_private_het(self.obs_variant)

        # Variant carrying a public classification (allele-linked), no observation.
        self.cls_variant = slowly_create_test_variant("3", 2000, "A", "T", self.grch37)
        self.allele = create_mock_allele(self.cls_variant, self.grch37)
        self.public_classification = self._make_public_classification(self.allele, clin_sig="5")

    # ---- helpers ----------------------------------------------------------
    def _make_private_het(self, variant):
        # samples_zygosity indexed [proband, mother, father]; only proband called (HET).
        samples_zygosity = Zygosity.HET + Zygosity.MISSING + Zygosity.MISSING
        n = len(samples_zygosity)
        CohortGenotype.objects.create(
            collection=self.cgc, variant=variant, het_count=1,
            samples_zygosity=samples_zygosity,
            samples_allele_depth=[20] * n, samples_allele_frequency=[100] * n,
            samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
            samples_phred_likelihood=[0] * n)
        # Global VZC row (stage-1 gate) - reflects the global count incl. private samples.
        VariantZygosityCount.objects.create(collection=self.global_collection, variant=variant, het_count=1)

    def _make_public_classification(self, allele, clin_sig):
        vc = Classification.create(
            user=self.user, lab=self.lab,
            data={SpecialEKeys.GENE_SYMBOL: {'value': 'BRCA1'}},
            save=True, source=SubmissionSource.API)
        vc.allele = allele
        vc.save()
        ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
        ClassificationModification.objects.create(
            classification=vc, user=self.user, source="TEST",
            published_evidence={"clinical_significance": {"value": clin_sig}},
            share_level=ShareLevel.PUBLIC.value, is_last_published=True, published=True)
        return vc

    def _query(self, variant, user=None, granularity=None):
        client = APIClient()
        if user is not None:
            client.force_authenticate(user)
        vc = variant.coordinate
        params = {
            "referenceName": vc.chrom, "start": vc.position - 1,
            "referenceBases": vc.ref, "alternateBases": vc.alt, "assemblyId": "GRCh37",
        }
        if granularity:
            params["requestedGranularity"] = granularity
        return client.get(reverse("beacon_g_variants"), params)

    @staticmethod
    def _result_set(data, dataset_id):
        for rs in data.get("response", {}).get("resultSets", []):
            if rs["id"] == dataset_id:
                return rs
        return None

    # ---- absence ----------------------------------------------------------
    def test_absent_variant(self):
        absent = slowly_create_test_variant("3", 9999, "A", "T", self.grch37)
        data = self._query(absent, user=self.user, granularity="record").json()
        self.assertFalse(data["responseSummary"]["exists"])

    # ---- permission tiering (most important, §11) -------------------------
    def test_private_observation_hidden_from_anonymous(self):
        data = self._query(self.obs_variant, user=None, granularity="record").json()
        obs = self._result_set(data, OBSERVATIONS_DATASET_ID)
        self.assertFalse(obs["exists"])

    def test_private_observation_visible_to_authorized_user(self):
        data = self._query(self.obs_variant, user=self.user, granularity="record").json()
        self.assertTrue(data["responseSummary"]["exists"])
        obs = self._result_set(data, OBSERVATIONS_DATASET_ID)
        self.assertTrue(obs["exists"])

    @override_settings(BEACON_MIN_REPORTABLE_COUNT=1)
    def test_observation_count_reported_above_floor(self):
        data = self._query(self.obs_variant, user=self.user, granularity="count").json()
        self.assertTrue(data["responseSummary"]["exists"])
        self.assertEqual(data["responseSummary"]["numTotalResults"], 1)

    @override_settings(BEACON_MIN_REPORTABLE_COUNT=5)
    def test_small_count_suppressed_by_floor(self):
        # count of 1 is below the floor: observations resultSet drops to boolean presence.
        data = self._query(self.obs_variant, user=self.user, granularity="record").json()
        obs = self._result_set(data, OBSERVATIONS_DATASET_ID)
        self.assertTrue(obs["exists"])
        self.assertNotIn("resultsCount", obs)

    # ---- granularity ------------------------------------------------------
    def test_boolean_granularity_hides_counts_and_records(self):
        data = self._query(self.cls_variant, user=self.user, granularity="boolean").json()
        self.assertIn("exists", data["responseSummary"])
        self.assertNotIn("numTotalResults", data["responseSummary"])
        self.assertNotIn("response", data)

    # ---- classification dataset (§5.5) ------------------------------------
    def test_public_classification_visible_to_anonymous(self):
        data = self._query(self.cls_variant, user=None, granularity="count").json()
        self.assertTrue(data["responseSummary"]["exists"])
        self.assertEqual(data["responseSummary"]["numTotalResults"], 1)

    def test_classification_record_tier_payload(self):
        data = self._query(self.cls_variant, user=None, granularity="record").json()
        cls = self._result_set(data, CLASSIFICATIONS_DATASET_ID)
        self.assertTrue(cls["exists"])
        self.assertEqual(cls["resultsCount"], 1)
        record = cls["results"][0]
        self.assertIn("clinicalInterpretations", record)
        self.assertEqual(record["contact"]["href"], "mailto:curator@lab.org")

    def test_per_dataset_resultsets_present(self):
        data = self._query(self.cls_variant, user=self.user, granularity="record").json()
        ids = {rs["id"] for rs in data["response"]["resultSets"]}
        self.assertEqual(ids, {OBSERVATIONS_DATASET_ID, CLASSIFICATIONS_DATASET_ID})

    # ---- audit ------------------------------------------------------------
    def test_query_writes_audit_row(self):
        self._query(self.cls_variant, user=None, granularity="count")
        row = BeaconInboundQuery.objects.latest("created")
        self.assertFalse(row.authenticated)
        self.assertTrue(row.classifications_exists)

    def test_missing_params_returns_400(self):
        response = APIClient().get(reverse("beacon_g_variants"), {"assemblyId": "GRCh37"})
        self.assertEqual(response.status_code, 400)
