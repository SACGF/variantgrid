from django.db import migrations, models

BAM = "B"
CRAM = "C"


def _set_cram_file_types(apps, _schema_editor):
    """ CRAMs were sent up as 'bam_file' and copied onto the sample as a BAM """
    AlignmentFile = apps.get_model("seqauto", "AlignmentFile")
    SampleFilePath = apps.get_model("snpdb", "SampleFilePath")
    AlignmentFile.objects.filter(path__iendswith=".cram").update(file_type=CRAM)
    SampleFilePath.objects.filter(file_type=BAM, file_path__iendswith=".cram").update(file_type=CRAM)


def _link_alignment_files_to_samples(apps, _schema_editor):
    """ Only a sequencing sample with exactly 1 BAM, present when its VCF imported, was copied onto the Sample """
    AlignmentFile = apps.get_model("seqauto", "AlignmentFile")
    SampleFromSequencingSample = apps.get_model("seqauto", "SampleFromSequencingSample")
    SampleFilePath = apps.get_model("snpdb", "SampleFilePath")

    existing = set(SampleFilePath.objects.values_list("sample_id", "file_path"))
    paths_by_sequencing_sample = {}
    for sequencing_sample_id, path, file_type in AlignmentFile.objects.values_list("sequencing_sample_id", "path",
                                                                                   "file_type"):
        paths_by_sequencing_sample.setdefault(sequencing_sample_id, []).append((path, file_type))

    new_sample_file_paths = []
    for sample_id, sequencing_sample_id in SampleFromSequencingSample.objects.values_list("sample_id",
                                                                                         "sequencing_sample_id"):
        for path, file_type in paths_by_sequencing_sample.get(sequencing_sample_id, []):
            if (sample_id, path) not in existing:
                new_sample_file_paths.append(SampleFilePath(sample_id=sample_id, file_path=path, file_type=file_type))
    if new_sample_file_paths:
        print(f"Linking {len(new_sample_file_paths)} alignment files to samples")
        SampleFilePath.objects.bulk_create(new_sample_file_paths, batch_size=2000)


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0051_one_off_tso500_tumor_fraction_percent"),
        ("snpdb", "0278_remove_user_award_titles"),
    ]

    operations = [
        migrations.RenameModel("BamFile", "AlignmentFile"),
        migrations.RenameField("flagstats", "bam_file", "alignment_file"),
        migrations.RenameField("singlesamplevcf", "bam_file", "alignment_file"),
        migrations.RenameField("qc", "bam_file", "alignment_file"),
        migrations.AddField(
            model_name="alignmentfile",
            name="file_type",
            field=models.CharField(choices=[("B", "BAM"), ("C", "CRAM"), ("E", "BED"), ("V", "VCF")], default="B",
                                   max_length=1),
        ),
        migrations.RunPython(_set_cram_file_types, migrations.RunPython.noop),
        migrations.RunPython(_link_alignment_files_to_samples, migrations.RunPython.noop),
    ]
