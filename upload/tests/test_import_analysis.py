import json
import os
import tempfile

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.analysis_import_export import analysis_export_to_dict
from analysis.models import Analysis, AnalysisEdge, MergeNode, SampleNode
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.models.models_enums import ImportSource
from upload.import_task_factories.import_task_factories import AnalysisImportTaskFactory
from upload.models import FileUpload, UploadedAnalysis, UploadedFileTypes
from upload.tasks.import_analysis_task import ImportAnalysisTask
from upload.uploaded_file_type import get_import_task_factory_from_extension


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class ImportAnalysisTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="import_analysis_user", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

    def _write_json(self, data) -> str:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            return f.name

    def _exported_analysis_file(self) -> str:
        analysis = Analysis(genome_build=self.grch37, name="exported analysis")
        analysis.set_defaults_and_save(self.user)
        sample_node = SampleNode.objects.create(analysis=analysis)
        AnalysisEdge.objects.create(parent=sample_node, child=MergeNode.objects.create(analysis=analysis))
        return self._write_json(analysis_export_to_dict(analysis))

    def test_factory_claims_exported_analysis(self):
        filename = self._exported_analysis_file()
        factory = get_import_task_factory_from_extension(self.user, filename, "json")
        self.assertIsInstance(factory, AnalysisImportTaskFactory)

    def test_factory_ignores_unrelated_json(self):
        filename = self._write_json({"some": "other file"})
        self.assertIsNone(get_import_task_factory_from_extension(self.user, filename, "json"))

    def test_process_items_creates_analysis(self):
        filename = self._exported_analysis_file()
        file_upload = FileUpload.objects.create(user=self.user, name=os.path.basename(filename),
                                                path=filename, file_type=UploadedFileTypes.ANALYSIS,
                                                import_source=ImportSource.COMMAND_LINE)

        self.assertEqual(ImportAnalysisTask.process_items(file_upload), 1)

        uploaded_analysis = UploadedAnalysis.objects.get(file_upload=file_upload)
        self.assertEqual(uploaded_analysis.analysis.analysisnode_set.count(), 2)
        self.assertEqual(uploaded_analysis.analysis.user, self.user)
        self.assertEqual(uploaded_analysis.genome_build, self.grch37)
        self.assertEqual(uploaded_analysis.get_data_url(), uploaded_analysis.analysis.get_absolute_url())
