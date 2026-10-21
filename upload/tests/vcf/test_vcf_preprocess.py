"""
Tests for the pieces of preprocess_vcf that make the unsorted-VCF sort-and-retry work (#127):
how the pipe is built for each pass, and the reset that puts the pipeline back to a clean state
in between. The full pipe needs bcftools and a real reference fasta, so TestRunPipe runs run_pipe on toy
stages instead, checking which stage a failure is blamed on (#1005).
"""
import copy
import os
import shlex
import shutil
import tempfile
from subprocess import CalledProcessError
from unittest import skipUnless
from unittest.mock import patch

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from snpdb.models import VCF, GenomeBuild, ImportStatus
from snpdb.models.models_enums import ImportSource
from upload.enums import UploadedFileTypes, UploadStepTaskType
from upload.models import (
    FileUpload,
    SimpleVCFImportInfo,
    ToolVersion,
    UploadedVCF,
    UploadPipeline,
    UploadStep,
)
from upload.vcf.vcf_preprocess import (
    SORT_VCF_SUB_STEP,
    UNSORTED_VCF_MESSAGE,
    VCF_CLEAN_AND_FILTER_SUB_STEP,
    PipeStageError,
    PreprocessFiles,
    _build_pipe_commands,
    _get_preprocess_files,
    _reset_for_retry,
    create_sub_step,
    preprocess_vcf,
    run_pipe,
)

_TEST_DIR = tempfile.mkdtemp(prefix="test_vcf_preprocess_")
_TEST_FASTA = os.path.join(_TEST_DIR, "placeholder.fna")


def _fake_annotation_settings() -> dict:
    annotation = copy.deepcopy(settings.ANNOTATION)
    annotation[settings.BUILD_GRCH37]["reference_fasta"] = _TEST_FASTA
    return annotation


@override_settings(ANNOTATION=_fake_annotation_settings())
class TestPreprocessVCFPipe(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        with open(_TEST_FASTA, "w", encoding="utf-8") as f:
            f.write(">placeholder\nN\n")  # Only its path is used - it goes to 'bcftools norm --fasta-ref'
        with open(_TEST_FASTA + ".fai", "w", encoding="utf-8") as f:
            for accession, length in cls.genome_build.standard_contigs.values_list("refseq_accession", "length"):
                f.write(f"{accession}\t{length}\t0\t60\t61\n")

        user = User.objects.create_user(username="preprocess_test_user", password="x")
        vcf_filename = os.path.join(_TEST_DIR, "test.vcf")
        with open(vcf_filename, "w", encoding="utf-8") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        file_upload = FileUpload.objects.create(user=user, name="test.vcf", path=vcf_filename,
                                                file_type=UploadedFileTypes.VCF,
                                                import_source=ImportSource.COMMAND_LINE)
        cls.upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        vcf = VCF.objects.create(name="test.vcf", genotype_samples=1, genome_build=cls.genome_build,
                                 import_status=ImportStatus.SUCCESS, user=user, date=timezone.now())
        UploadedVCF.objects.create(file_upload=file_upload, upload_pipeline=cls.upload_pipeline, vcf=vcf)
        cls.upload_step = UploadStep.objects.create(upload_pipeline=cls.upload_pipeline, name="Preprocess VCF",
                                                    sort_order=0, task_type=UploadStepTaskType.CELERY,
                                                    input_filename=vcf_filename)
        cls.vcf_filename = vcf_filename

    def setUp(self):
        super().setUp()
        # create_sub_step shells out to 'bcftools --version', which isn't a given in a test environment
        tool_version, _ = ToolVersion.objects.get_or_create(name="bcftools", version="test")
        patcher = patch("upload.vcf.vcf_preprocess.get_bcftools_tool_version", return_value=tool_version)
        patcher.start()
        self.addCleanup(patcher.stop)

        self.files = _get_preprocess_files(self.upload_pipeline, self.vcf_filename,
                                           os.path.join(_TEST_DIR, "cleaned.header.vcf"), disable_swap=False)
        self.addCleanup(self._remove_files)

    def _remove_files(self):
        for filename in self.files.written_by_pipe:
            if os.path.exists(filename):
                os.remove(filename)

    @staticmethod
    def _clean_and_filter_args(pipe_commands) -> list[str]:
        return pipe_commands[VCF_CLEAN_AND_FILTER_SUB_STEP]

    def test_first_pass_checks_sorted_and_has_no_sort_stage(self):
        pipe_commands, sub_steps = _build_pipe_commands(self.upload_step, self.files)
        args = self._clean_and_filter_args(pipe_commands)
        self.assertIn("--unsorted-marker-file", args)
        self.assertIn(self.files.unsorted_marker, args)
        self.assertNotIn("--allow-unsorted", args)
        self.assertNotIn(SORT_VCF_SUB_STEP, pipe_commands)
        self.assertNotIn(SORT_VCF_SUB_STEP, sub_steps)

    def test_sort_pass_adds_sort_stage_after_clean_and_filter(self):
        pipe_commands, sub_steps = _build_pipe_commands(self.upload_step, self.files, sort_vcf=True)
        args = self._clean_and_filter_args(pipe_commands)
        self.assertIn("--allow-unsorted", args)
        self.assertNotIn("--unsorted-marker-file", args)

        stage_names = list(pipe_commands)
        self.assertEqual(stage_names.index(SORT_VCF_SUB_STEP),
                         stage_names.index(VCF_CLEAN_AND_FILTER_SUB_STEP) + 1)
        self.assertEqual(stage_names.index(SORT_VCF_SUB_STEP) + 1,
                         stage_names.index(UploadStep.NORMALIZE_SUB_STEP))

        sort_cmd = pipe_commands[SORT_VCF_SUB_STEP]
        self.assertIn("sort", sort_cmd)
        temp_dir = sort_cmd[sort_cmd.index("--temp-dir") + 1]
        # Trailing separator, otherwise bcftools uses it as a temp filename prefix rather than a directory
        self.assertTrue(temp_dir.endswith(os.sep))
        self.assertTrue(os.path.isdir(temp_dir))

        # Its runtime/stderr shows on the upload page alongside the other tool steps
        self.assertEqual(sub_steps[SORT_VCF_SUB_STEP].name, SORT_VCF_SUB_STEP)

    def test_reset_for_retry_clears_sub_steps_split_files_and_stats(self):
        _, sub_steps = _build_pipe_commands(self.upload_step, self.files)
        self.assertEqual(self.upload_step.substep_set.count(), len(sub_steps))

        split_filename = os.path.join(self.files.split_vcf_dir, "test00.vcf.gz")
        for filename in [split_filename, self.files.unsorted_marker, self.files.skipped_records_stats]:
            with open(filename, "w", encoding="utf-8") as f:
                f.write("partial\n")

        _reset_for_retry(self.files, sub_steps)

        self.assertEqual(self.upload_step.substep_set.count(), 0)
        self.assertFalse(os.path.exists(split_filename))
        self.assertFalse(os.path.exists(self.files.unsorted_marker))
        self.assertFalse(os.path.exists(self.files.skipped_records_stats))

    def test_only_one_normalize_sub_step_after_reset_and_rebuild(self):
        """ ModifiedImportedVariants.get_for_pipeline() does a .get() on the normalize step """
        _, sub_steps = _build_pipe_commands(self.upload_step, self.files)
        _reset_for_retry(self.files, sub_steps)
        _build_pipe_commands(self.upload_step, self.files, sort_vcf=True)

        qs = self.upload_pipeline.uploadstep_set.filter(name=UploadStep.NORMALIZE_SUB_STEP)
        self.assertEqual(qs.count(), 1)

    def test_swap_skipped_file_only_when_disabling_swap(self):
        files = _get_preprocess_files(self.upload_pipeline, self.vcf_filename,
                                      self.files.cleaned_vcf_header, disable_swap=True)
        self.assertIn(files.swap_skipped, files.written_by_pipe)
        pipe_commands, _ = _build_pipe_commands(self.upload_step, files, disable_swap=True)
        clean_alts_cmd = pipe_commands["vcf_clean_alts"]
        self.assertIn("--disable-swap", clean_alts_cmd)
        self.assertIn(files.swap_skipped, clean_alts_cmd)

        self.assertIsNone(self.files.swap_skipped)
        self.assertNotIn("--disable-swap", _build_pipe_commands(self.upload_step, self.files)[0]["vcf_clean_alts"])


    # ------------------------------------------------------------------ retry orchestration
    def _run_preprocess_vcf(self, run_pipe_side_effect):
        with patch("upload.vcf.vcf_preprocess.run_pipe", side_effect=run_pipe_side_effect) as mock_run_pipe:
            preprocess_vcf(self.upload_step)
        return mock_run_pipe

    def _unsorted_then_ok(self):
        """ First pass writes the marker and dies (as vcf_clean_and_filter does), second pass succeeds """
        calls = []

        def side_effect(pipe_commands, *_args, **_kwargs):
            calls.append(pipe_commands)
            if len(calls) == 1:
                with open(self.files.unsorted_marker, "w", encoding="utf-8") as f:
                    f.write("VCF is not sorted. Line 7: 1:100 comes after 1:200.\n")
                raise CalledProcessError(1, "vcf_clean_and_filter")

        return side_effect, calls

    def test_unsorted_marker_triggers_sort_pass(self):
        side_effect, calls = self._unsorted_then_ok()
        self._run_preprocess_vcf(side_effect)

        self.assertEqual(len(calls), 2)
        self.assertNotIn(SORT_VCF_SUB_STEP, calls[0])
        self.assertIn(SORT_VCF_SUB_STEP, calls[1])

    def test_sort_pass_warns_the_user_records_were_reordered(self):
        side_effect, _ = self._unsorted_then_ok()
        self._run_preprocess_vcf(side_effect)

        qs = SimpleVCFImportInfo.objects.filter(upload_step__upload_pipeline=self.upload_pipeline)
        self.assertEqual([i.message_string for i in qs], [UNSORTED_VCF_MESSAGE])
        self.assertIn(UNSORTED_VCF_MESSAGE, [w.message for w in self.upload_pipeline.get_warnings()])

    def test_sorted_vcf_runs_the_pipe_once_with_no_warning(self):
        mock_run_pipe = self._run_preprocess_vcf(None)

        self.assertEqual(mock_run_pipe.call_count, 1)
        self.assertNotIn(SORT_VCF_SUB_STEP, mock_run_pipe.call_args.args[0])
        self.assertFalse(SimpleVCFImportInfo.objects.filter(message_string=UNSORTED_VCF_MESSAGE).exists())

    def test_failure_without_marker_is_raised(self):
        with self.assertRaises(CalledProcessError):
            self._run_preprocess_vcf(CalledProcessError(1, "bcftools norm"))
        self.assertFalse(SimpleVCFImportInfo.objects.filter(message_string=UNSORTED_VCF_MESSAGE).exists())

    def test_stale_marker_does_not_turn_another_failure_into_a_sort_pass(self):
        with open(self.files.unsorted_marker, "w", encoding="utf-8") as f:
            f.write("left behind by an earlier attempt\n")

        with self.assertRaises(CalledProcessError):
            self._run_preprocess_vcf(CalledProcessError(1, "bcftools norm"))
        self.assertFalse(SimpleVCFImportInfo.objects.filter(message_string=UNSORTED_VCF_MESSAGE).exists())


class TestPreprocessFiles(TestCase):
    def test_written_by_pipe_skips_absent_swap_file(self):
        files = PreprocessFiles(vcf_filename="in.vcf", cleaned_vcf_header="header.vcf", split_vcf_dir="split",
                                unsorted_marker="unsorted.txt", skipped_contigs_stats="contigs.tsv",
                                skipped_records_stats="records.tsv", skipped_filters_stats="filters.tsv",
                                skipped_alts_stats="alts.tsv", converted_alts_stats="converted.tsv")
        self.assertEqual(len(files.written_by_pipe), 6)


@override_settings(IMPORT_PROCESSING_DIR=os.path.join(_TEST_DIR, "import_processing"))
class TestRunPipe(TestCase):
    """ The real run_pipe on toy stages - no bcftools or fasta needed except for the #1179 repro """
    VCF_FILENAME = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf", "grch38_brca1.vcf")
    NO_FORMAT_HEADER_VCF_FILENAME = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf",
                                                 "no_format_header.grch38_brca1.vcf")

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.create_user(username="run_pipe_test_user", password="x")
        file_upload = FileUpload.objects.create(user=user, name="test.vcf", path=cls.VCF_FILENAME,
                                                file_type=UploadedFileTypes.VCF,
                                                import_source=ImportSource.COMMAND_LINE)
        cls.upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        cls.upload_step = UploadStep.objects.create(upload_pipeline=cls.upload_pipeline, name="Preprocess VCF",
                                                    sort_order=0, task_type=UploadStepTaskType.CELERY)

    def _run_pipe(self, pipe_commands: dict) -> dict:
        sub_steps = {name: create_sub_step(self.upload_step, name, command, tool_version=None)
                     for name, command in pipe_commands.items()}
        run_pipe(pipe_commands, sub_steps, dict(os.environ), self.upload_pipeline)
        return sub_steps

    def _failed_pipe(self, pipe_commands: dict) -> PipeStageError:
        with self.assertRaises(PipeStageError) as cm:
            self._run_pipe(pipe_commands)
        return cm.exception

    def test_failing_stage_is_named_with_its_own_stderr_and_code(self):
        pipe_commands = {"a": ["cat", self.VCF_FILENAME],
                         "b": ["sh", "-c", "cat > /dev/null; echo boom >&2; exit 3"],
                         "c": ["cat"]}
        e = self._failed_pipe(pipe_commands)
        self.assertEqual(e.stage.name, "b")
        self.assertEqual(e.returncode, 3)
        self.assertEqual(e.stderr, "boom")
        self.assertEqual(e.cmd, shlex.join(pipe_commands["b"]))
        self.assertIn("'b' (2 of 3)", str(e))
        self.assertIn("a=0 b=3 c=0", str(e))

        sub_steps = {s.name: s for s in self.upload_step.substep_set.all()}
        self.assertEqual(sub_steps["b"].error_message, "boom")
        self.assertEqual(sub_steps["a"].output_text, "")
        self.assertFalse(sub_steps["a"].error_message)

    def test_upstream_sigpipe_and_broken_pipe_are_not_blamed(self):
        e = self._failed_pipe({"yes": ["yes"], "fail": ["sh", "-c", "exit 3"], "c": ["cat"]})
        self.assertEqual([s.returncode for s in e.stages], [141, 3, 0])
        self.assertEqual(e.stage.name, "fail")

        # Python ignores SIGPIPE, so a manage.py stage exits 1 with a traceback instead of 141
        e = self._failed_pipe({"python": ["python3", "-c", "print('x\\n' * 10**7)"],
                               "fail": ["sh", "-c", "exit 3"], "c": ["cat"]})
        self.assertIn("BrokenPipeError", e.stages[0].stderr)
        self.assertEqual(e.stage.name, "fail")

    def test_success_writes_every_sub_steps_stderr(self):
        sub_steps = self._run_pipe({"a": ["cat", self.VCF_FILENAME],
                                    "b": ["sh", "-c", "echo note >&2; cat"],
                                    "c": ["cat"]})
        self.assertEqual(sub_steps["b"].output_text, "note")
        self.assertEqual(sub_steps["a"].output_text, "")
        self.assertFalse(sub_steps["b"].error_message)

    @skipUnless(shutil.which(settings.BCFTOOLS_COMMAND), "bcftools not installed")
    def test_bcftools_norm_dying_on_undeclared_format_is_blamed(self):
        """ #1179: bcftools norm dies mid-pipe on FORMAT tags missing from the header """
        e = self._failed_pipe({
            "cat": ["cat", self.NO_FORMAT_HEADER_VCF_FILENAME],
            UploadStep.NORMALIZE_SUB_STEP: [settings.BCFTOOLS_COMMAND, "norm", "--multiallelics=-", "-"],
            "remove_header": [settings.BCFTOOLS_COMMAND, "view", "--no-header", "-"],
            "c": ["cat"],
        })
        self.assertEqual(e.stage.name, UploadStep.NORMALIZE_SUB_STEP)
        self.assertIn("[E::vcf_format]", e.stderr)
