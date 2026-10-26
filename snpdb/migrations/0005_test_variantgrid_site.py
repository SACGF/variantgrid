from django.db import migrations

TEST_VARIANTGRID_SITE_ID = 10


def _create_test_variantgrid_site(apps, _schema_editor):
    """ test.variantgrid.com's own Site, so links built outside a request use its domain rather than localhost.
        Created only if the id is free - other deployments may have made their own rows by hand """
    Site = apps.get_model("sites", "Site")
    Site.objects.get_or_create(pk=TEST_VARIANTGRID_SITE_ID,
                               defaults={"domain": "test.variantgrid.com", "name": "VariantGrid Test"})


class Migration(migrations.Migration):

    dependencies = [
        ('sites', '0002_alter_domain_unique'),
        ('snpdb', '0004_alphagenome_splicing_column'),
    ]

    operations = [
        migrations.RunPython(_create_test_variantgrid_site, migrations.RunPython.noop),
    ]
