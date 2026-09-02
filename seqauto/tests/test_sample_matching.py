"""Tests for matching VCF samples to SampleSheet.csv rows (SACGF/variantgrid_sapath#440)."""
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from seqauto.models import (
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingSample,
)
from seqauto.models.models_enums import DataGeneration
from seqauto.models.models_seqauto import get_samples_by_sequencing_sample
from snpdb.models import VCF, Sample


class GetSamplesBySequencingSampleTests(TestCase):
    def setUp(self):
        seq_model, _ = SequencerModel.objects.get_or_create(
            model="MiSeq", data_naming_convention=DataGeneration.MISEQ)
        sequencer, _ = Sequencer.objects.get_or_create(name="MS1", sequencer_model=seq_model)
        sequencing_run = SequencingRun.objects.create(name="RUN_001", sequencer=sequencer)
        self.sample_sheet = SampleSheet.objects.create(
            sequencing_run=sequencing_run, path="/data/RUN_001/SampleSheet.csv", hash="SHA1")
        self.user = User.objects.create(username="testuser")

    def _make_sequencing_samples(self, sample_names):
        return [SequencingSample.objects.create(sample_sheet=self.sample_sheet, sample_id=name,
                                                sample_name=name, sample_number=i, barcode="ACGT")
                for i, name in enumerate(sample_names, start=1)]

    def _make_vcf(self, name, sample_names):
        vcf = VCF.objects.create(name=name, date=timezone.now(), user=self.user,
                                 genotype_samples=len(sample_names))
        samples = [Sample.objects.create(vcf=vcf, name=sname, vcf_sample_name=sname)
                   for sname in sample_names]
        return vcf, samples

    def test_exact_and_sample_number_suffix(self):
        sequencing_samples = self._make_sequencing_samples(["ABC-1", "ABC-2"])
        vcf, samples = self._make_vcf("RUN_001", ["ABC_1", "ABC-2_S2"])
        expected = dict(zip(sequencing_samples, samples))
        self.assertEqual(get_samples_by_sequencing_sample(sequencing_samples, vcf), expected)

    def test_wrong_sample_number_not_matched(self):
        sequencing_samples = self._make_sequencing_samples(["ABC-1"])
        vcf, _ = self._make_vcf("RUN_001", ["ABC-1_S7"])
        self.assertEqual(get_samples_by_sequencing_sample(sequencing_samples, vcf), {})

    def test_prefix_of_another_sample_not_matched(self):
        """ The old startswith match linked "ABC1_extra" to sequencing sample "ABC1" """
        sequencing_samples = self._make_sequencing_samples(["ABC1"])
        vcf, _ = self._make_vcf("RUN_001", ["ABC1_extra"])
        self.assertEqual(get_samples_by_sequencing_sample(sequencing_samples, vcf), {})

    def test_single_sample_falls_back_to_vcf_name(self):
        """ SpliceGirls VCFs name their only sample "SAMPLE" """
        sequencing_samples = self._make_sequencing_samples(["ABC-1"])
        vcf, samples = self._make_vcf("ABC_1", ["SAMPLE"])
        self.assertEqual(get_samples_by_sequencing_sample(sequencing_samples, vcf),
                         {sequencing_samples[0]: samples[0]})
