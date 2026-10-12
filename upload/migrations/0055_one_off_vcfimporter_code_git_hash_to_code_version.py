from django.db import migrations
from django.db.models import Min


def _code_git_hash_to_code_version(apps, _schema_editor):
    """ One CodeVersion per hash, created in the order each hash was first used """
    CodeVersion = apps.get_model("snpdb", "CodeVersion")
    VCFImporter = apps.get_model("upload", "VCFImporter")

    first_used = VCFImporter.objects.values("code_git_hash").annotate(first_pk=Min("pk")).order_by("first_pk")
    for code_git_hash in first_used.values_list("code_git_hash", flat=True):
        code_version = CodeVersion.objects.get_or_create(git_hash=code_git_hash)[0]
        VCFImporter.objects.filter(code_git_hash=code_git_hash).update(code_version=code_version)


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0054_vcfimporter_code_version'),
    ]

    operations = [
        migrations.RunPython(_code_git_hash_to_code_version, migrations.RunPython.noop, elidable=True),
    ]
