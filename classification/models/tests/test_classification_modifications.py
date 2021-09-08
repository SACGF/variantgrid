from django.test import TestCase

from classification.enums import SubmissionSource, SpecialEKeys
from classification.models.classification import Classification, ClassificationModification, VCBlobKeys
from classification.models.classification_ref import ClassificationRef
from classification.models.classification_utils import ClassificationJsonParams
from classification.models.tests.test_utils import ClassificationTestUtils
from library.guardian_utils import admin_bot


class ClassificationTestCaseModifications(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def ensure_latest_version(self, vc: Classification):
        latest_ver = ClassificationRef.init_from_obj(vc.user, vc.last_edited_version)
        latest_json = latest_ver.as_json(ClassificationJsonParams(current_user=vc.user, include_data=True, flatten=True))['data']

        editable_json = vc.as_json(ClassificationJsonParams(current_user=vc.user, include_data=True, flatten=True))['data']

        self.assertEqual(latest_json, editable_json, 'Latest version json does not match editable json')

    def test_immutable(self):
        """
        Test uploading values from API, make sure they get marked as immutable
        """
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = None
        try:
            vc = Classification.create(
                user=user,
                lab=lab,
                lab_record_id=None,
                data={
                    SpecialEKeys.C_HGVS: {'value': 'c.301A>C'},
                    SpecialEKeys.G_HGVS: {'value': '5678'},
                    'foo': {'value': 'bar'}
                },
                save=True,
                source=SubmissionSource.API,
                make_fields_immutable=False)

            c_hgvs = vc.evidence.get(SpecialEKeys.C_HGVS)
            self.assertEqual(c_hgvs.get('immutable'), 'variantgrid')

            # Added validation message to C_HGVS
            # have't changed G_HGVS at all but still want it to become immutable
            patch_immutable = {
                SpecialEKeys.C_HGVS: None,
                SpecialEKeys.G_HGVS: {'note': 'bar'}
            }

            vc.patch_value(
                patch=patch_immutable,
                clear_all_fields=False,
                user=admin_bot(),
                source=SubmissionSource.API,
                leave_existing_values=False,
                save=True,
                make_patch_fields_immutable=True
            )
            g_hgvs = vc.evidence.get(SpecialEKeys.G_HGVS) or {}
            self.assertEqual(g_hgvs.get('value'), "5678")
            self.assertEqual(g_hgvs.get('note'), "bar")
            self.assertEqual(g_hgvs.get('immutable'), SubmissionSource.VARIANT_GRID)

            c_hgvs = vc.evidence.get(SpecialEKeys.C_HGVS) or {}
            self.assertEqual(c_hgvs.get('value'), 'c.301A>C')
            # At time cant test this as a process is adding extra validation
            # self.assertEqual(len(c_hgvs.get('validation', [])), 1)
            self.assertEqual(c_hgvs.get('immutable'), SubmissionSource.VARIANT_GRID)

            foo_val = vc.evidence.get('foo')
            self.assertEqual(foo_val.get('value'), 'bar')
            self.assertEqual(foo_val.get('immutable'), None)

        finally:
            if vc:
                vc.delete()

    def test_nullify(self):
        """
        Test that uploading null wipes out value, note, explain, xrefs but sets immutability
        """
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = None
        try:
            vc = Classification.create(
                user=user,
                lab=lab,
                lab_record_id=None,
                data={
                    SpecialEKeys.C_HGVS: {'value': 'g.301A>C'},
                    SpecialEKeys.G_HGVS: {'value': '5678'}
                },
                save=True,
                source=SubmissionSource.API,
                make_fields_immutable=False)

            vc.patch_value(
                patch={SpecialEKeys.REFERENCE_DEPTH: {'value': '12', 'note': 'hello there', 'explain': 'Explain this'}},
                user=user,
                source=SubmissionSource.API,
                save=True,
                make_patch_fields_immutable=True
            )

            vc.patch_value(
                patch={SpecialEKeys.REFERENCE_DEPTH: {'explain': 'Explain this 2'}},
                user=user,
                source=SubmissionSource.API,
                save=True,
                make_patch_fields_immutable=True
            )
            ref_depth = vc.evidence.get(SpecialEKeys.REFERENCE_DEPTH)
            # value, note, explain are all merged
            self.assertEqual(ref_depth, {'value': 12, 'note': 'hello there', 'explain': 'Explain this 2', 'immutable': SubmissionSource.API})

            # try adding PMID so we also have db_refs to deal with and make sure they get cleared out
            vc.patch_value(
                patch={SpecialEKeys.REFERENCE_DEPTH: {'value': 'PMID:123456', 'note': 'hello there', 'explain': 'Explain this'}},
                user=user,
                source=SubmissionSource.API,
                save=True,
                make_patch_fields_immutable=True
            )

            vc.patch_value(
                patch={SpecialEKeys.REFERENCE_DEPTH: None},
                user=user,
                source=SubmissionSource.API,
                save=True,
                make_patch_fields_immutable=True
            )
            ref_depth = vc.evidence.get(SpecialEKeys.REFERENCE_DEPTH)
            self.assertEqual(ref_depth, {'immutable': SubmissionSource.API})

        finally:
            if vc:
                vc.delete()

    def test_modifying(self):
        """
        Tests some general modifications
        - ensure immutability is respected
        - a nullify
        - test using overwrite
        :return:
        """
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = None
        try:
            vc = Classification.create(
                user=user,
                lab=lab,
                lab_record_id=None,
                data={"x": 1},
                save=True,
                source=SubmissionSource.API,
                make_fields_immutable=True)

            # Test initial creation
            initial_version_qs = ClassificationModification.objects.filter(classification=vc)
            self.assertEqual(initial_version_qs.count(), 1, 'Should be 1 and only 1 modification')
            x = vc.evidence.get('x')
            self.assertEqual(x.get('value'), "1")
            self.assertEqual(x.get('immutable'), SubmissionSource.API)

            self.ensure_latest_version(vc)

            # Test regular patch
            vc.patch_value(
                patch={"a": {"value":"2"}},
                clear_all_fields=False,
                user=user,
                source=SubmissionSource.API,
                leave_existing_values=False,
                save=True,
                make_patch_fields_immutable=True
            )
            x = vc.evidence.get('x')
            self.assertEqual(x.get('value'), "1")
            self.assertEqual(x.get('immutable'), SubmissionSource.API)
            a = vc.evidence.get('a')
            self.assertEqual(a.get('value'), "2")
            self.assertEqual(a.get('immutable'), SubmissionSource.API)

            self.ensure_latest_version(vc)

            # Test patch that wont override existing values
            vc.patch_value(
                patch={"x": 3, "b": 4},
                clear_all_fields=False,
                user=user,
                source=SubmissionSource.VARIANT_GRID,
                leave_existing_values=True,
                save=True,
                make_patch_fields_immutable=False
            )
            x = vc.evidence.get('x')
            self.assertEqual(x.get('value'), "1")
            self.assertEqual(x.get('immutable'), SubmissionSource.API)

            b = vc.evidence.get('b')
            self.assertEqual(b.get('value'), "4")
            self.assertEqual(b.get('immutable'), None)

            self.ensure_latest_version(vc)

            # Test nulling values
            vc.patch_value(
                patch={"a": None, "x": {"value": None}},
                clear_all_fields=False,
                user=user,
                source=SubmissionSource.VARIANT_GRID,
                leave_existing_values=False,
                save=True,
                make_patch_fields_immutable=True
            )
            self.assertEqual('a' in vc.evidence, False, f'Found a as {vc.evidence.get("a")}')

            x = vc.evidence.get('x')
            self.assertEqual(x.get('value'), None)
            self.assertEqual(x.get('immutable'), SubmissionSource.VARIANT_GRID)

            self.ensure_latest_version(vc)

            # Test resetting
            vc.patch_value(
                patch={"c": True},
                clear_all_fields=True,
                user=user,
                source=SubmissionSource.VARIANT_GRID,
                leave_existing_values=False,
                save=True,
                make_patch_fields_immutable=False
            )
            self.assertEqual('x' in vc.evidence, False)
            self.assertEqual('b' in vc.evidence, False)
            c = vc.evidence.get('c')
            self.assertEqual(c.get('value'), 'True')

            self.ensure_latest_version(vc)
        finally:
            if vc:
                vc.delete()

    def xtest_strip_html(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={},
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=True)

        vc.patch_value(
            patch={
                # single text has HTML stripped
                SpecialEKeys.CONDITION: 'Big <b>bold</b> Toe',
                # text areas keep their HTML
                SpecialEKeys.LITERATURE: 'Some <i>italic</i> literature'
            },
            user=user,
            source=SubmissionSource.API,
            save=True
        )
        self.assertEqual('Big bold Toe', vc.get(SpecialEKeys.CONDITION))
        self.assertEqual('Some <i>italic</i> literature', vc.get(SpecialEKeys.LITERATURE))

    def test_test_mode(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                SpecialEKeys.C_HGVS: 'g.301A>C',
                SpecialEKeys.ZYGOSITY: 'zany'
            },
            save=False,
            source=SubmissionSource.API,
            make_fields_immutable=True)
        self.assertEqual(vc.id, None)

    def test_regexes(self):
        """
        Basic test to make sure refs are being scanned for - there is a dedicated
        test for parsing the text.
        """
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                SpecialEKeys.C_HGVS: 'g.301A>C',
                SpecialEKeys.LITERATURE: {'value': 'PMID: 123456', 'note': 'PMID: 555444'}
            },
            save=False,
            source=SubmissionSource.API,
            make_fields_immutable=True)
        literature = vc.evidence.get(SpecialEKeys.LITERATURE)
        db_refs = literature.get(VCBlobKeys.DB_REFS.value)
        self.assertEqual(len(db_refs), 2)
        puby_1 = db_refs[0]
        puby_1.pop('internal_id')
        self.assertEqual(puby_1, {
            'db': 'PubMed',
            'id': 'PubMed: 123456',
            'idx': '123456',
            'url': 'https://www.ncbi.nlm.nih.gov/pubmed/?term=123456'
        })

        puby_2 = db_refs[1]
        puby_2.pop('internal_id')
        self.assertEqual(puby_2, {
            'db': 'PubMed',
            'id': 'PubMed: 555444',
            'idx': '555444',
            'url': 'https://www.ncbi.nlm.nih.gov/pubmed/?term=555444'
        })
