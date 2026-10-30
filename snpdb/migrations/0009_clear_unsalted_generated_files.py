"""
Graphs and pedigree charts cached before their file names were salted with SECRET_KEY keep being served from
MEDIA_ROOT under names worked out from object ids. Delete those rows and files so they regenerate under the
salted name on next view.
"""
import os

from django.conf import settings
from django.db import migrations

PEDIGREE_CHART = "PedigreeChart"


def _remove_file(filename):
    if filename and os.path.isfile(filename):
        os.unlink(filename)


def _clear_unsalted_generated_files(apps, _schema_editor):
    CachedGeneratedFile = apps.get_model("snpdb", "CachedGeneratedFile")

    # A graph's file used to be <GENERATED_DIR>/<generator>/<params_hash>.png, now it's a salted hmac of params_hash.
    # A pedigree chart's params_hash itself is now salted, so the old ones can't be told apart: clear them all
    stale_ids = []
    generators = {PEDIGREE_CHART}
    for cgf in CachedGeneratedFile.objects.exclude(filename__isnull=True).only("pk", "generator", "params_hash",
                                                                              "filename"):
        generator_dir = os.path.join(settings.GENERATED_DIR, cgf.generator)
        unsalted_graph = (os.path.dirname(cgf.filename) == generator_dir
                          and os.path.basename(cgf.filename).split(".")[0] == cgf.params_hash)
        if cgf.generator == PEDIGREE_CHART or unsalted_graph:
            stale_ids.append(cgf.pk)
            generators.add(cgf.generator)
            _remove_file(cgf.filename)
    CachedGeneratedFile.objects.filter(pk__in=stale_ids).delete()
    CachedGeneratedFile.objects.filter(generator=PEDIGREE_CHART).delete()

    # Files left behind by rows that were deleted or never saved their filename
    for generator in generators:
        generator_dir = os.path.join(settings.GENERATED_DIR, generator)
        if not os.path.isdir(generator_dir):
            continue
        kept = set(CachedGeneratedFile.objects.filter(generator=generator).exclude(filename__isnull=True)
                   .values_list("filename", flat=True))
        for name in os.listdir(generator_dir):
            filename = os.path.join(generator_dir, name)
            if filename not in kept:
                _remove_file(filename)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0008_dbsnp_bigint_pk'),
    ]

    operations = [
        migrations.RunPython(_clear_unsalted_generated_files, reverse_code=migrations.RunPython.noop),
    ]
