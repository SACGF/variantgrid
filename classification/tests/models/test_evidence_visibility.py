from unittest.mock import patch

from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

from classification.enums.classification_enums import ShareLevel
from classification.models.classification import Classification, EvidenceVisibility
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from classification.tests.models.test_utils import ClassificationTestUtils
from library.django_utils.unittest_utils import production_query_count
from snpdb.models.models import Lab, Organization

HIDDEN = {'value': "(hidden)", 'hidden': True}


class EvidenceVisibilityTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, self.lab_user = ClassificationTestUtils.lab_and_user()
        self.lab = Lab.objects.select_related('organization').get(pk=lab.pk)
        other_org = Organization.objects.create(name='OtherOrg', group_name='otherorg')
        self.other_lab = Lab.objects.create(name='Other', organization=other_org, city='CityC',
                                            country=self.lab.country, group_name='otherorg/other')
        self.outsider = User.objects.create(username='outsider')
        self.outsider.groups.add(self.other_lab.group)

        e_key_map = EvidenceKeyMap(ordered_keys=[
            EvidenceKey(key='lab_only', max_share_level=ShareLevel.LAB.key),
            EvidenceKey(key='shared', max_share_level=ShareLevel.ALL_USERS.key),
        ])
        patcher = patch.object(EvidenceKeyMap, 'instance', return_value=e_key_map)
        patcher.start()
        self.addCleanup(patcher.stop)

        self.evidence = {'lab_only': {'value': 'secret'}, 'shared': {'value': 'open'}}

    def test_hides_from_users_outside_lab(self):
        classification = Classification(lab=self.lab, user=self.lab_user)
        outsider_view = EvidenceVisibility(self.outsider).visible_evidence(classification, self.evidence)
        self.assertEqual(outsider_view, {'lab_only': HIDDEN, 'shared': {'value': 'open'}})

        lab_view = EvidenceVisibility(self.lab_user).visible_evidence(classification, self.evidence)
        self.assertEqual(lab_view, self.evidence)

    def test_owner_sees_own_classification(self):
        classification = Classification(lab=self.lab, user=self.outsider)
        self.assertEqual(EvidenceVisibility(self.outsider).visible_evidence(classification, self.evidence), self.evidence)

    def test_lab_config_overrides_share_level(self):
        self.other_lab.classification_settings.classification_config = {'lab_only': {'max_share_level': ShareLevel.ALL_USERS.key}}
        visibility = EvidenceVisibility(self.lab_user)
        other_lab_classification = Classification(lab=self.other_lab, user=self.outsider)
        self.assertEqual(visibility.visible_evidence(other_lab_classification, self.evidence), self.evidence)

        outsider_classification = Classification(lab=self.lab, user=self.lab_user)
        outsider_view = EvidenceVisibility(self.outsider).visible_evidence(outsider_classification, self.evidence)
        self.assertEqual(outsider_view['lab_only'], HIDDEN)

    def test_queries_once_per_user(self):
        visibility = EvidenceVisibility(self.outsider)
        classifications = [Classification(lab=self.lab, user=self.lab_user) for _ in range(5)]
        with CaptureQueriesContext(connection) as ctx:
            for classification in classifications:
                visibility.visible_evidence(classification, self.evidence)
        self.assertEqual(production_query_count(ctx.captured_queries), 1)
