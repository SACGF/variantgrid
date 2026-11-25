from django.db import migrations, models
from django.db.models import F, Q

from manual.operations.manual_operations import ManualOperation


def _test_has_whole_chromosome_candidates(apps):
    """ A <DUP>/<DEL> over half its contig - the candidates one_off_whole_chromosome_variants checks """
    Variant = apps.get_model("snpdb", "Variant")
    half_contig = F("locus__contig__length") / 2
    return Variant.objects.filter(alt__seq__in=["<DUP>", "<DEL>"]).filter(
        Q(svlen__gt=half_contig) | Q(svlen__lt=-half_contig)).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0001_initial'),
        ('snpdb', '0018_one_off_split_shared_user_create_labs'),
    ]

    operations = [
        migrations.AlterField(
            model_name='allelemergelog',
            name='allele_linking_tool',
            field=models.CharField(choices=[('SC', 'Identical Contig/Version'), ('CA', 'ClinGen Allele Registry'), ('DB', 'dbSNP API'), ('NR', 'NCBI Remap'), ('PC', 'Picard LiftoverVCF'), ('CM', 'CrossMap'), ('BL', 'BCFtools/liftover'), ('WC', 'Same chromosome')], max_length=2),
        ),
        migrations.AlterField(
            model_name='liftoverrun',
            name='conversion_tool',
            field=models.CharField(choices=[('SC', 'Identical Contig/Version'), ('CA', 'ClinGen Allele Registry'), ('DB', 'dbSNP API'), ('NR', 'NCBI Remap'), ('PC', 'Picard LiftoverVCF'), ('CM', 'CrossMap'), ('BL', 'BCFtools/liftover'), ('WC', 'Same chromosome')], max_length=2),
        ),
        migrations.AlterField(
            model_name='variantallele',
            name='allele_linking_tool',
            field=models.CharField(choices=[('SC', 'Identical Contig/Version'), ('CA', 'ClinGen Allele Registry'), ('DB', 'dbSNP API'), ('NR', 'NCBI Remap'), ('PC', 'Picard LiftoverVCF'), ('CM', 'CrossMap'), ('BL', 'BCFtools/liftover'), ('WC', 'Same chromosome')], max_length=2),
        ),
        ManualOperation.operation_manage(["one_off_whole_chromosome_variants"],
                                         test=_test_has_whole_chromosome_candidates),
    ]
