# Generated by Django 4.2.10 on 2024-06-20 03:05
import logging
import os

from django.db import migrations

from library.utils import sha256sum_str, file_sha256sum


def _md5_to_sha256(apps, _schema_editor):
    Sequence = apps.get_model('snpdb', 'Sequence')
    GenomeFasta = apps.get_model('snpdb', 'GenomeFasta')

    num_seq_records = Sequence.objects.count()
    logging.info("Updating hash on %d Sequence records (may take a few mins)", num_seq_records)
    seq_records = []
    for seq in Sequence.objects.all():
        seq.seq_sha256_hash = sha256sum_str(seq.seq)
        seq_records.append(seq)

    if seq_records:
        Sequence.objects.bulk_update(seq_records, fields=['seq_sha256_hash'], batch_size=2000)

    gf_records = []
    for gf in GenomeFasta.objects.all():
        if os.path.exists(gf.index_filename):
            gf.index_sha256sum = file_sha256sum(gf.index_filename)
        else:
            logging.warning("Genome fasta index: '%s' does not exist", gf.index_filename)
            gf.index_sha256sum = "old-md5sum-hash-" + gf.index_md5sum
        gf_records.append(gf)
    if gf_records:
        GenomeFasta.objects.bulk_update(gf_records, fields=['index_sha256sum'])


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0136_genomefasta_index_sha256sum_sequence_seq_sha256_hash'),
    ]

    operations = [
        migrations.RunPython(_md5_to_sha256)
    ]
