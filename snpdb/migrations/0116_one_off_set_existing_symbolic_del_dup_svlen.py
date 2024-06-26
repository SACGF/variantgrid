# Generated by Django 4.2.9 on 2024-02-14 04:31

from django.db import migrations
from django.db.models import F, OuterRef, Subquery


def _update_symbolic_svlen(Variant, allele_seq: str, calc_svlen: F):
    print(f"Updating {allele_seq}")
    variant_subquery = Variant.objects.filter(pk=OuterRef("pk")).annotate(calc_svlen=calc_svlen).values("calc_svlen")[:1]
    Variant.objects.filter(alt__seq=allele_seq).update(svlen=Subquery(variant_subquery))


def _one_off_set_existing_symbolic_del_dup_svlen(apps, _schema_editor):
    # We have symbolic variants in the DB (del/dup) that have their length encoded via start/end
    # These were read into the system as end = variant.POS + abs(svlen_info)

    Variant = apps.get_model("snpdb", "Variant")
    _update_symbolic_svlen(Variant, "<DUP>", F("end") - F("locus__position"))
    _update_symbolic_svlen(Variant, "<DEL>", F("locus__position") - F("end"))


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0115_alter_variant_unique_together_variant_svlen_and_more'),
    ]

    operations = [
        migrations.RunPython(_one_off_set_existing_symbolic_del_dup_svlen)
    ]
