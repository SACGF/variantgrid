# Generated by Django 4.2.10 on 2024-06-20 03:53
import os

from django.db import migrations

from library.utils import file_sha256sum


def _md5_to_sha256(apps, _schema_editor):
    DBNSFPGeneAnnotationVersion = apps.get_model('annotation', 'DBNSFPGeneAnnotationVersion')
    ClinVarVersion = apps.get_model('annotation', 'ClinVarVersion')
    HumanProteinAtlasAnnotationVersion = apps.get_model('annotation', 'HumanProteinAtlasAnnotationVersion')

    DBNSFP_GENE_MD5_TO_SHA256 = {
        # This is the only one we've used really - hasn't changed in ages
        "2313e5c5fcb5557ef31d677922cb4ae9": "b00c30e51d761e7c46744ba228ddb36b2d10b7774ccf6684ea519d007966a374",
    }
    dgav_records = []
    for dgav in DBNSFPGeneAnnotationVersion.objects.all():
        if precalc := DBNSFP_GENE_MD5_TO_SHA256.get(dgav.md5_hash):
            dgav.sha256_hash = precalc
        else:
            # Will just cause a reload and new version if someone re-imports, no big deal
            dgav.sha256_hash = "old-md5sum-hash-" + dgav.md5_hash
        dgav_records.append(dgav)

    if dgav_records:
        DBNSFPGeneAnnotationVersion.objects.bulk_update(dgav_records, ['sha256_hash'])

    CVV_MD5_TO_SHA256 = {

    }
    cvv_records = []
    for cvv in ClinVarVersion.objects.all():
        if precalc := CVV_MD5_TO_SHA256.get(cvv.md5_hash):
            cvv.sha256_hash = precalc
        elif os.path.exists(cvv.filename):
            cvv.sha256_hash = file_sha256sum(cvv.filename)
        else:
            # Will just cause a reload and new version if someone re-imports, no big deal
            cvv.sha256_hash = "old-md5sum-hash-" + cvv.md5_hash
        cvv_records.append(cvv)

    if cvv_records:
        ClinVarVersion.objects.bulk_update(cvv_records, ['sha256_hash'])

    HPA_MD5_TO_SHA256 = {
        "1ca25285e651e770baf33d122cfbc02a": "134cc4665e45677ce7c4c4df5af392230022ee33a62bade71b01d93e8687c392",
        "f61101faafd5794192ad489a678d2a94": "592ebfdddbe5bee5c7ff3d4a3a2e42d923efcf42bfa46da442355469f638513e",
    }
    hpa_records = []
    for hpa in HumanProteinAtlasAnnotationVersion.objects.all():
        if precalc := HPA_MD5_TO_SHA256.get(hpa.md5_hash):
            hpa.sha256_hash = precalc
        elif os.path.exists(hpa.filename):
            hpa.sha256_hash = file_sha256sum(hpa.filename)
        else:
            # Will just cause a reload and new version if someone re-imports, no big deal
            hpa.sha256_hash = "old-md5sum-hash-" + hpa.md5_hash
        hpa_records.append(hpa)

    if hpa_records:
        HumanProteinAtlasAnnotationVersion.objects.bulk_update(hpa_records, ['sha256_hash'])


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0104_clinvarversion_sha256_hash_and_more'),
    ]

    operations = [
        migrations.RunPython(_md5_to_sha256)
    ]