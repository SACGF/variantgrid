from django.test import TestCase

from classification.enums import SpecialEKeys, SubmissionSource
from classification.enums.classification_enums import ShareLevel
from classification.models.classification import Classification, ClassificationModification
from classification.tests.models.test_utils import ClassificationTestUtils
from mme.serializers.patient_profile import mme_eligible_classifications


class EligibilityTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        # Eligibility requires the owning lab to have opted into MME with a contact.
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()

    def _classification_with_mod(self, share_level: str, is_last_published: bool = True,
                                 withdrawn: bool = False) -> Classification:
        vc = Classification.create(
            user=self.user, lab=self.lab,
            data={SpecialEKeys.GENE_SYMBOL: {'value': 'BRCA1'}},
            save=True, source=SubmissionSource.API)
        if withdrawn:
            vc.withdrawn = True
            vc.save()
        # Deterministically control the "latest published" modification.
        ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
        ClassificationModification.objects.create(
            classification=vc, user=self.user, source="TEST", delta={},
            share_level=share_level, is_last_published=is_last_published, published=True)
        return vc

    def test_public_included(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        eligible = mme_eligible_classifications()
        self.assertIn(vc.pk, {cm.classification_id for cm in eligible})

    def test_all_users_excluded(self):
        vc = self._classification_with_mod(ShareLevel.ALL_USERS.value)
        self.assertNotIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})

    def test_lab_excluded(self):
        vc = self._classification_with_mod(ShareLevel.LAB.value)
        self.assertNotIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})

    def test_unpublished_public_excluded(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value, is_last_published=False)
        self.assertNotIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})

    def test_withdrawn_public_excluded(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value, withdrawn=True)
        self.assertNotIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})

    def test_lab_not_mme_enabled_excluded(self):
        # PUBLIC record, but the owning lab has not opted into MME -> excluded.
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.lab.mme_enabled = False
        self.lab.save()
        self.assertNotIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})

    def test_lab_mme_enabled_toggle_includes(self):
        # Flipping the same PUBLIC record's lab back on includes it again.
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.lab.mme_enabled = False
        self.lab.save()
        self.assertNotIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})
        self.lab.mme_enabled = True
        self.lab.save()
        self.assertIn(vc.pk, {cm.classification_id for cm in mme_eligible_classifications()})
