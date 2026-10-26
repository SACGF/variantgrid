import configparser
import re

from django.conf import settings
from django.db import migrations


def _populate_mutational_signature_info(apps, _schema_editor):
    """ From the signature descriptions shipped in MUTATIONAL_SIGNATURE_DATA_DIR """
    MutationalSignatureInfo = apps.get_model("annotation", "MutationalSignatureInfo")

    config = configparser.ConfigParser(interpolation=None)
    config.read(settings.MUTATIONAL_SIGNATURE_INFO_FILE)
    signatures = []
    for section in config.sections():
        if m := re.match(r"Signature (\d+)", section):
            sdata = config[section]
            comments = sdata["Comments"]
            if comments == "N/A":
                comments = ''
            signatures.append(MutationalSignatureInfo(signature_id=int(m.group(1)),
                                                      cancer_types=sdata["Cancer types"],
                                                      proposed_aetiology=sdata["Proposed aetiology"],
                                                      additional_mutational_features=sdata["Additional mutational features"],
                                                      comments=comments))
    MutationalSignatureInfo.objects.bulk_create(signatures)


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('annotation', '0002_initial'),
    ]

    operations = [
        migrations.RunPython(_populate_mutational_signature_info),
    ]
