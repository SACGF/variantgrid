from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from classification.enums import ShareLevel, SpecialEKeys, SubmissionSource
from classification.models import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from library.guardian_utils import all_users_group, assign_permission_to_user_and_groups
from patients.models_enums import Zygosity
from snpdb.models import (
    CohortGenotype,
    CohortGenotypeCollection,
    GenomeBuild,
    Sample,
    VariantZygosityCount,
    VariantZygosityCountCollection,
    VCFFilter,
)
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant
from snpdb.variant_sample_information import VariantSampleGenotypes


class VariantSampleGenotypesTest(TestCase):
    """ The JSON payload the client side samples grid is drawn from """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='vsg_test_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

        # create_fake_cohort packs samples as [proband, mother, father]
        cls.cohort_37 = create_fake_cohort(cls.user, cls.grch37)
        cls.cohort_38 = create_fake_cohort(cls.user, cls.grch38)

        cls.other_user = User.objects.get_or_create(username='vsg_test_other_user')[0]
        cls.other_cohort_37 = create_fake_cohort(cls.other_user, cls.grch37)

    def _create_cohort_genotype(self, cohort, variant, samples_zygosity, filters=None, samples_filters=None):
        cgc = CohortGenotypeCollection.objects.get(cohort=cohort)
        n = len(samples_zygosity)
        return CohortGenotype.objects.create(collection=cgc, variant=variant, filters=filters,
                                             samples_zygosity=samples_zygosity, samples_filters=samples_filters,
                                             samples_allele_depth=[20] * n, samples_allele_frequency=[100] * n,
                                             samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
                                             samples_phred_likelihood=[0] * n)

    def _het_proband_and_mother(self, cohort, variant, **kwargs):
        samples_zygosity = Zygosity.HET + Zygosity.HET + Zygosity.MISSING
        return self._create_cohort_genotype(cohort, variant, samples_zygosity, **kwargs)

    def _rows_by_sample_name(self, data) -> dict:
        return {row["sample_name"]: row for row in data["rows"]}

    def test_rows_restricted_to_visible_samples(self):
        variant = slowly_create_test_variant("3", 3000, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)
        self._het_proband_and_mother(self.other_cohort_37, variant)

        data = VariantSampleGenotypes(self.user, variant).to_json()
        self.assertEqual({"total": 4, "visible": 2, "invisible": 2}, data["observations"])
        self.assertEqual(2, len(data["rows"]), "Only the VCF this user can read")
        vcf_ids = {Sample.objects.get(pk=row["sample"]).vcf_id for row in data["rows"]}
        self.assertEqual({self.cohort_37.vcf_id}, vcf_ids)
        self.assertEqual({Zygosity.HET: 2}, data["zygosity_counts"])

    def test_locus_counts_cover_hidden_samples(self):
        """ Counts are of everything at the locus - only rows are permission filtered """
        variant = slowly_create_test_variant("3", 3100, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)
        self._het_proband_and_mother(self.other_cohort_37, variant)

        data = VariantSampleGenotypes(self.user, variant).to_json()
        locus_counts = data["locus_counts"]
        self.assertEqual(1, len(locus_counts))
        this_variant = locus_counts[0]
        self.assertEqual("This variant", this_variant["description"])
        self.assertEqual(4, this_variant["total"])
        self.assertEqual(4, this_variant["HET"])

    def test_filters_pass(self):
        variant = slowly_create_test_variant("3", 3200, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)

        data = VariantSampleGenotypes(self.user, variant).to_json()
        row = self._rows_by_sample_name(data)["proband"]
        self.assertEqual("PASS", row["filters"])
        self.assertIsNone(row["sample_filters"], "No FT in the VCF")

    def test_filters_single_and_multiple(self):
        vcf = self.cohort_37.vcf
        VCFFilter.objects.create(vcf=vcf, filter_code="Y", filter_id="LOW_DEPTH")

        single = slowly_create_test_variant("3", 3300, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, single, filters="X")
        multiple = slowly_create_test_variant("3", 3400, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, multiple, filters="XY")

        single_data = VariantSampleGenotypes(self.user, single).to_json()
        self.assertEqual("YOUSHALLNOTPASS", self._rows_by_sample_name(single_data)["proband"]["filters"])

        multiple_data = VariantSampleGenotypes(self.user, multiple).to_json()
        self.assertEqual("YOUSHALLNOTPASS,LOW_DEPTH", self._rows_by_sample_name(multiple_data)["proband"]["filters"])

    def test_sample_filters(self):
        variant = slowly_create_test_variant("3", 3500, "A", "T", self.grch37)
        # MISSING_FT_VALUE = no FT for that sample, "" = FT of PASS
        samples_filters = ["X", "", CohortGenotype.MISSING_FT_VALUE]
        self._het_proband_and_mother(self.cohort_37, variant, samples_filters=samples_filters)

        data = VariantSampleGenotypes(self.user, variant).to_json()
        row = self._rows_by_sample_name(data)["proband"]
        self.assertEqual("YOUSHALLNOTPASS", row["sample_filters"])

    def test_filter_codes_decoded_per_vcf(self):
        """ Shared contig - one response, two VCFs, and the same code means different things in each """
        variant = slowly_create_test_variant("MT", 3000, "A", "T", self.grch37)
        self.assertEqual({self.grch37, self.grch38}, variant.genome_builds)

        # Both VCFs use code "X", but for a different filter id
        VCFFilter.objects.filter(vcf=self.cohort_38.vcf, filter_code="X").update(filter_id="OFF_TARGET")
        # Make the second VCF's proband visible too, so we get a row from each build
        assign_permission_to_user_and_groups(self.user, Sample.objects.get(vcf=self.cohort_38.vcf, name="proband"))

        self._het_proband_and_mother(self.cohort_37, variant, filters="X")
        self._het_proband_and_mother(self.cohort_38, variant, filters="X")

        data = VariantSampleGenotypes(self.user, variant).to_json()
        filters_by_build = {row["genome_build"]: row["filters"] for row in data["rows"]}
        self.assertEqual({self.grch37.name: "YOUSHALLNOTPASS", self.grch38.name: "OFF_TARGET"}, filters_by_build)

    def test_allele_frequency_unit_for_graphs(self):
        """ The displayed AF may be percent formatted text - graphs need the unit number """
        variant = slowly_create_test_variant("3", 3700, "A", "T", self.grch37)
        cgc = CohortGenotypeCollection.objects.get(cohort=self.cohort_37)
        CohortGenotype.objects.create(collection=cgc, variant=variant,
                                      samples_zygosity=Zygosity.HET + Zygosity.HET + Zygosity.MISSING,
                                      samples_allele_frequency=[0.25, CohortGenotype.MISSING_NUMBER_VALUE, -1],
                                      samples_allele_depth=[20] * 3, samples_read_depth=[30] * 3,
                                      samples_genotype_quality=[30] * 3, samples_phred_likelihood=[0] * 3)

        rows = self._rows_by_sample_name(VariantSampleGenotypes(self.user, variant).to_json())
        self.assertEqual(0.25, rows["proband"]["allele_frequency_unit"])
        self.assertIsNone(rows["mother"]["allele_frequency_unit"], "Missing AF is left off the graph")

    def test_row_limit(self):
        """ max_rows bounds what we materialise, and leaves every count exact """
        variant = slowly_create_test_variant("3", 3750, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)
        self._het_proband_and_mother(self.other_cohort_37, variant)

        data = VariantSampleGenotypes(self.user, variant, max_rows=1).to_json()
        self.assertTrue(data["truncated"])
        self.assertEqual(1, len(data["rows"]))

        for max_rows in [None, 2, 100]:
            all_data = VariantSampleGenotypes(self.user, variant, max_rows=max_rows).to_json()
            self.assertFalse(all_data["truncated"], f"max_rows={max_rows} isn't truncated")
            self.assertEqual(2, len(all_data["rows"]), f"max_rows={max_rows} materialises all visible rows")

    def test_counts_never_depend_on_row_limit(self):
        """ Counting reads the packed arrays, so capping rows must not move a single count.
            Beacon reports these numbers, and it asks for no rows at all. """
        variant = slowly_create_test_variant("3", 3800, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)
        self._het_proband_and_mother(self.other_cohort_37, variant)

        def counts(max_rows):
            data = VariantSampleGenotypes(self.user, variant, max_rows=max_rows).to_json()
            return {k: data[k] for k in ["observations", "zygosity_counts", "locus_counts", "samples"]}

        uncapped = counts(None)
        self.assertEqual({"total": 4, "visible": 2, "invisible": 2}, uncapped["observations"])
        self.assertEqual({Zygosity.HET: 2}, uncapped["zygosity_counts"])
        self.assertEqual(4, uncapped["locus_counts"][0]["HET"])
        for max_rows in [0, 1, 2]:
            self.assertEqual(uncapped, counts(max_rows), f"counts with max_rows={max_rows}")

    def test_stop_when_full_uses_global_counts(self):
        """ Stopping the scan early means we no longer know the permission-scoped counts, so the
            payload reports the precomputed global ones and says nothing about what's visible """
        variant = slowly_create_test_variant("3", 3900, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)
        self._het_proband_and_mother(self.other_cohort_37, variant)
        vzcc, _ = VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION, defaults={"description": "test"})
        VariantZygosityCount.objects.create(collection=vzcc, variant=variant, het_count=4)

        vsg = VariantSampleGenotypes(self.user, variant, max_rows=1, stop_when_full=True)
        self.assertTrue(vsg.partial)
        data = vsg.to_json()
        self.assertTrue(data["truncated"])
        self.assertEqual(1, len(data["rows"]))
        self.assertEqual(4, data["observations"]["total"], "Global count, not what we scanned")
        self.assertIsNone(data["observations"]["visible"], "We stopped before counting these")
        self.assertIsNone(data["observations"]["invisible"])
        self.assertEqual(4, data["locus_counts"][0]["HET"], "Locus counts come from the global collection")

    def test_stop_when_full_completes_when_it_fits(self):
        """ Under the cap there's nothing to stop for, so the counts are the exact scanned ones """
        variant = slowly_create_test_variant("3", 3950, "A", "T", self.grch37)
        self._het_proband_and_mother(self.cohort_37, variant)
        self._het_proband_and_mother(self.other_cohort_37, variant)

        vsg = VariantSampleGenotypes(self.user, variant, max_rows=100, stop_when_full=True)
        self.assertFalse(vsg.partial)
        data = vsg.to_json()
        self.assertEqual({"total": 4, "visible": 2, "invisible": 2}, data["observations"])

    def test_no_observations(self):
        variant = slowly_create_test_variant("3", 3600, "A", "T", self.grch37)

        data = VariantSampleGenotypes(self.user, variant).to_json()
        self.assertEqual({"total": 0, "visible": 0, "invisible": 0}, data["observations"])
        self.assertEqual([], data["rows"])
        self.assertEqual([self.grch37.name], data["genome_builds"])


class VariantSampleGenotypesClassificationsTest(TestCase):
    """ Classification pills, attached to the row for the sample that was classified """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        ClassificationTestUtils.setUp()
        cls.lab, cls.lab_user = ClassificationTestUtils.lab_and_user()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

        cls.user = User.objects.get_or_create(username='vsg_classification_user')[0]
        cls.user.groups.add(all_users_group())  # Every logged in user is in this
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.samples = {s.name: s for s in cls.cohort.get_samples()}
        for sample in cls.samples.values():
            assign_permission_to_user_and_groups(cls.user, sample)
            assign_permission_to_user_and_groups(cls.lab_user, sample)

        cls.variant = slowly_create_test_variant("3", 6000, "A", "T", cls.grch37)
        create_mock_allele(cls.variant, cls.grch37)
        cgc = CohortGenotypeCollection.objects.get(cohort=cls.cohort)
        CohortGenotype.objects.create(collection=cgc, variant=cls.variant,
                                      samples_zygosity=Zygosity.HET * 3,
                                      samples_allele_depth=[20] * 3, samples_allele_frequency=[100] * 3,
                                      samples_read_depth=[30] * 3, samples_genotype_quality=[30] * 3,
                                      samples_phred_likelihood=[0] * 3)

    def _create_classification(self, sample, clinical_significance, share_level):
        classification = Classification.create(
            user=self.lab_user,
            lab=self.lab,
            data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": clinical_significance},
                  SpecialEKeys.ALLELE_ORIGIN: {"value": "germline"}},  # Otherwise it gets a somatic pill too
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False,
            variant=self.variant,
            sample=sample)
        classification.publish_latest(user=self.lab_user, share_level=share_level)
        return classification

    def _classifications_by_sample_name(self, user) -> dict:
        data = VariantSampleGenotypes(user, self.variant).to_json()
        return {row["sample_name"]: row["classifications"] for row in data["rows"]}

    def test_classification_attached_to_its_sample(self):
        classification = self._create_classification(self.samples["proband"], "VUS", ShareLevel.ALL_USERS)

        by_sample_name = self._classifications_by_sample_name(self.user)
        self.assertEqual([], by_sample_name["mother"], "Only the classified sample gets pills")
        proband_pills = by_sample_name["proband"]
        self.assertEqual(1, len(proband_pills))
        pill = proband_pills[0]
        self.assertEqual(classification.pk, pill["id"])
        self.assertEqual("VUS (3)", pill["label"], "Pretty value from the evidence key")
        self.assertEqual("cs cs-vus", pill["css_class"])
        self.assertEqual(str(self.lab), pill["lab"])

    def test_classification_user_cannot_read_is_absent(self):
        self._create_classification(self.samples["proband"], "VUS", ShareLevel.LAB)

        by_sample_name = self._classifications_by_sample_name(self.user)
        self.assertEqual([], by_sample_name["proband"], "Not shared outside the lab's organisation")
        lab_by_sample_name = self._classifications_by_sample_name(self.lab_user)
        self.assertEqual(1, len(lab_by_sample_name["proband"]), "Visible to the lab that curated it")
