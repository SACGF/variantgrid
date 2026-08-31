from django import forms
from django.contrib.auth.models import User
from django.test import TestCase

from library.forms import ROFormMixin
from snpdb.models import GenomeBuild
from snpdb.models.models_cohort import Cohort


class ReadOnlyDisplayForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = Cohort
        fields = ['name', 'user']
        read_only_display = ('user', )


class ROFormMixinTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user("ro_owner")
        User.objects.create_user("ro_someone_else")
        cls.cohort = Cohort.objects.create(name="ro cohort", user=cls.user,
                                           genome_build=GenomeBuild.grch38())

    def test_read_only_display_renders_value_not_choices(self):
        form = ReadOnlyDisplayForm(instance=self.cohort)
        user_html = str(form['user'])
        self.assertIn("ro_owner", user_html)
        self.assertNotIn("ro_someone_else", user_html)
        self.assertNotIn("<select", user_html)

    def test_read_only_display_disables_only_that_field(self):
        form = ReadOnlyDisplayForm(instance=self.cohort)
        self.assertTrue(form.fields['user'].disabled)
        self.assertFalse(form.fields['name'].disabled)
