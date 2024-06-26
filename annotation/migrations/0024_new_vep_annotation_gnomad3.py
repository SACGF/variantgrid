# Generated by Django 3.1 on 2021-02-16 05:03

from django.db import migrations

from library.django_utils import bulk_insert_class_data


def _new_vep_annotation_gnomad3(apps, _schema_editor):
    # Separate out gnomAD 2 vs 3
    GNOMAD_2 = 'g'
    GNOMAD_3 = 'n'

    # Make everything from gnomAD2 GRCh37 specific EXCEPT gnomadAF (still want that one)
    ColumnVEPField = apps.get_model("annotation", "ColumnVEPField")
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")

    grch37 = GenomeBuild.objects.get(pk="GRCh37")
    grch38 = GenomeBuild.objects.get(pk="GRCh38")

    # All existing gnomAD are now GRCh37 only (will insert new legacy one below)
    ColumnVEPField.objects.filter(vep_custom=GNOMAD_2).update(genome_build=grch37)

    # Might as well hide these now as we can - GRCh37 has 46, GRCh37 has 30
    ColumnVEPField.objects.filter(column__in=["phylop_46_way_mammalian", "phastcons_46_way_mammalian"]).update(genome_build=grch37)
    ColumnVEPField.objects.filter(column__in=["phylop_30_way_mammalian", "phastcons_30_way_mammalian"]).update(genome_build=grch38)

    COLUMN_VEP_FIELD = [
        # Legacy
        {'column': 'gnomad2_liftover_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_2, 'variant_grid_column_id': 'gnomad2_liftover_af', 'source_field': 'AF', 'category': 'F'},
        # gnomAD 3
        {'column': 'gnomad3_ac', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_ac', 'source_field': 'AC', 'category': 'F'},
        {'column': 'gnomad3_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_af', 'source_field': 'AF', 'category': 'F'},
        {'column': 'gnomad3_an', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_an', 'source_field': 'AN', 'category': 'F'},
        {'column': 'gnomad3_afr_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_afr_af', 'source_field': 'AF-afr', 'category': 'F'},
        {'column': 'gnomad3_amr_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_amr_af', 'source_field': 'AF-amr', 'category': 'F'},
        {'column': 'gnomad3_asj_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_asj_af', 'source_field': 'AF-asj', 'category': 'F'},
        {'column': 'gnomad3_eas_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_eas_af', 'source_field': 'AF-eas', 'category': 'F'},
        {'column': 'gnomad3_filtered', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_filtered', 'source_field': 'FILTER', 'category': 'F'},
        {'column': 'gnomad3_fin_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_fin_af', 'source_field': 'AF-fin', 'category': 'F'},
        {'column': 'gnomad3_hom_alt', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_hom_alt', 'source_field': 'nhomalt', 'category': 'F'},
        {'column': 'gnomad3_nfe_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_nfe_af', 'source_field': 'AF-nfe', 'category': 'F'},
        {'column': 'gnomad3_oth_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_oth_af', 'source_field': 'AF-oth', 'category': 'F'},
        {'column': 'gnomad3_popmax', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_popmax', 'source_field': 'popmax', 'category': 'F'},
        {'column': 'gnomad3_popmax_ac', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_popmax_ac', 'source_field': 'AC_popmax', 'category': 'F'},
        {'column': 'gnomad3_popmax_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_popmax_af', 'source_field': 'AF_popmax', 'category': 'F'},
        {'column': 'gnomad3_popmax_an', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_popmax_an', 'source_field': 'AN_popmax', 'category': 'F'},
        {'column': 'gnomad3_popmax_hom_alt', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_popmax_hom_alt', 'source_field': 'nhomalt_popmax', 'category': 'F'},
        {'column': 'gnomad3_sas_af', 'source_field_has_custom_prefix': True,
         'vep_custom': GNOMAD_3, 'variant_grid_column_id': 'gnomad_sas_af', 'source_field': 'AF-sas', 'category': 'F'},
    ]
    bulk_insert_class_data(apps, "annotation", [("ColumnVEPField", COLUMN_VEP_FIELD)])
    ColumnVEPField.objects.filter(vep_custom=GNOMAD_3).update(genome_build=grch38)


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0023_auto_20210216_1528'),
        ('snpdb', '0019_new_columns_gnomad3'),
    ]

    operations = [
        migrations.RunPython(_new_vep_annotation_gnomad3)
    ]
