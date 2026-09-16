from django.test import SimpleTestCase

from variantgrid.deployment_validation.tool_version_checks import check_vcf_split_pipe


class ToolVersionChecksTest(SimpleTestCase):
    def test_vcf_split_pipe(self):
        """ The real import split stage: GNU split --filter running bash | bgzip, header prepended to each chunk.
            Needs the bgzip binary (CI installs the tabix package) """
        self.assertTrue(check_vcf_split_pipe())
