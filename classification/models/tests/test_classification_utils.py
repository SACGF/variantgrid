from django.test import TestCase

from classification.models import ValidationMerger, EvidenceMixin, Classification


class ClassificationTestCaseUtils(TestCase):

    def test_validation_merger(self):
        """
        Test the ValidationMerger utility to add new errors, maintain old ones.
        Thinking the way we store/re-calculate validation is starting to have a bit of cold smell
        """
        vm = ValidationMerger()
        vm.tested({'a', 'b', 'c', 'd', 'e'}, 'bad')
        vm.tested({'a', 'b', 'c', 'd'}, 'good')

        vm.add_message(key='a', code='bad', message='xxx', severity='error')
        vm.add_message(key='b', code='bad', message='xxx', severity='error')
        vm.add_message(key='c', code='bad', message='xxx', severity='error')

        patch = {
            'a': None,
            'b': {'value': 'Pb', 'validation': [{'code': 'z_untested'}]},
            'c': {'value': 'Pc'},
        }
        evidence = {
            'a': {'value': 'Oa'},
            'd': {'value': 'Od', 'validation': [{'code': 'good'}, {'code': 'z_untested'}]},  # should be added as existing error is being removed
            'e': {'value': 'Od', 'validation': [{'code': 'good'}]}  # should not be patched as no new errors added or existing errors removed
        }
        vm.apply(patch, evidence)

        self.assertEqual(set(patch.keys()), {'a', 'b', 'c', 'd'})

        self.assertEqual(len(patch['d']['validation']), 1)  # should have good removed
        self.assertEqual(patch['d']['validation'][0]['code'], 'z_untested')  # but z_untested left alone

        self.assertEqual(patch['a'].get('value', None), None)  # remade entry with None value
        self.assertEqual(patch['a']['validation'][0]['code'], 'bad')  # and new bad validation

        self.assertEqual(len(patch['b']['validation']), 2)  # has z_untested and new bad

    def test_tidy(self):
        """
        Test the to_patch method, ensure literals get changed to value
        and that malformed keys are normalised
        :return:
        """
        patch_input = {
            "a": None,
            "x": {"value": None},
            "weird  *  key": {'value': 3, 'explain': 'it is weird'},
            17: True
        }
        output = EvidenceMixin.to_patch(patch_input)
        expected = {
            'a': None,
            'x': {'value': None},
            'weird_key': {'value': 3, 'explain': 'it is weird'},
            '17': {'value': True, 'explain': None, 'note': None}
        }
        self.assertEqual(output, expected)

    def test_double_tidy(self):
        """
        Test that to_patch has no effect on already to_patch data
        i.e. if the data is already good it should have no effect
        """
        rando_input = {
            'a': None,
            'b': 3,
            'c': {'value': 234},
            4: False,
            5: 0
        }
        expected = {
            'a': None,
            'b': {'value': 3, 'explain': None, 'note': None},
            'c': {'value': 234},
            '4': {'value': False, 'explain': None, 'note': None},
            '5': {'value': 0, 'explain': None, 'note': None}
        }
        tidy_once = EvidenceMixin.to_patch(rando_input)
        tidy_twice = EvidenceMixin.to_patch(tidy_once)
        self.assertEqual(tidy_once, expected)
        self.assertEqual(tidy_once, tidy_twice)

    def test_is_valid_transcript(self):
        self.assertTrue(Classification.is_supported_transcript("NM_002739.5(PRKCG):c.1397T>C"))
        self.assertFalse(Classification.is_supported_transcript("NX_023343.1(RNU4ATAC):n.50G>A"))
