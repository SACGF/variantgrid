from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from classification.classification_stats import (
    get_classification_counts,
    get_criteria_counts,
    get_grouped_classification_counts,
)
from classification.enums import ShareLevel, SpecialEKeys, SubmissionSource
from classification.models import Classification, EvidenceKeyMap, ImportedAlleleInfo
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.enums import ClinicalSignificance
from snpdb.models import Allele, GenomeBuild, GenomeBuildPatchVersion, Lab


@override_settings(CLASSIFICATION_STATS_USE_SHARED=True, CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True)
class ClassificationStatsTest(TestCase):
    """ The classification_graphs counts are done in SQL - these pin the counting rules they replaced """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.other_lab = Lab.objects.create(name="Labby 2", organization=self.lab.organization, city="CityA",
                                            country=self.lab.country, group_name="instx/labby2")
        self.user.groups.add(self.other_lab.group)
        self.viewer = User.objects.create(username="stats_viewer")

    def _classify(self, lab: Lab, allele: Allele, clinical_significance: str, **evidence_values) -> Classification:
        evidence = {SpecialEKeys.CLINICAL_SIGNIFICANCE: clinical_significance, SpecialEKeys.GENE_SYMBOL: "BRCA1",
                    **evidence_values}
        vc = Classification.create(user=self.user, lab=lab, data={k: {"value": v} for k, v in evidence.items()},
                                   save=True, source=SubmissionSource.API)
        vc.publish_latest(user=self.user, share_level=ShareLevel.ALL_USERS)
        allele_info = ImportedAlleleInfo.objects.create(
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                GenomeBuild.grch37()),
            imported_c_hgvs=f"NM_007294.4:c.{vc.pk}G>A", allele=allele)
        Classification.objects.filter(pk=vc.pk).update(allele_info=allele_info)
        return vc

    def test_allele_classified_by_two_labs_counts_once(self):
        """ The latest record represents the allele - the pie used to count it in each lab's slice """
        brca_allele = Allele.objects.create()
        self._classify(self.lab, brca_allele, "P")
        self._classify(self.other_lab, brca_allele, "VUS")
        self._classify(self.lab, Allele.objects.create(), "LP")

        counts = get_classification_counts(self.viewer, show_unclassified=False, unique_alleles=True)
        self.assertEqual(counts["Pathogenic"], 0)
        self.assertEqual(counts["VUS"], 1)
        self.assertEqual(counts["Likely Pathogenic"], 1)

        gene_data = get_grouped_classification_counts(self.viewer, "published_evidence", evidence_key="gene_symbol",
                                                      show_unclassified=False, allele_level=True)
        gene_counts = {series["name"]: dict(zip(series["x"], series["y"])) for series in gene_data}
        self.assertEqual(gene_counts["VUS"], {"BRCA1": 1})
        self.assertEqual(gene_counts["Pathogenic"], {"BRCA1": 0})
        self.assertEqual(gene_counts["Likely Pathogenic"], {"BRCA1": 1})

    def test_criteria_met_not_met_not_set(self):
        """ Met is any strength - NM, NA and N are "not met", a missing value is "not set" """
        key_for_label = {e.pretty_label: e.key for e in EvidenceKeyMap.instance().criteria_for("acmg")}
        pvs1, pm2 = key_for_label["PVS1"], key_for_label["PM2"]
        self._classify(self.lab, Allele.objects.create(), "P", **{pvs1: "PVS", pm2: "NM"})
        self._classify(self.lab, Allele.objects.create(), "P", **{pvs1: "PM", pm2: "N"})

        pathogenic = get_criteria_counts(self.viewer, "published_evidence")[ClinicalSignificance.PATHOGENIC]
        met, not_met, not_set = (dict(zip(series["x"], series["y"])) for series in pathogenic)
        self.assertEqual((met["PVS1"], not_met["PVS1"], not_set["PVS1"]), (1, 0, 0))
        self.assertEqual((met["PM2"], not_met["PM2"], not_set["PM2"]), (0, 1, 0))
        self.assertEqual((met["PS1"], not_met["PS1"], not_set["PS1"]), (0, 0, 1))
