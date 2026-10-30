from django.db import migrations

# SACGF/variantgrid#2195 - autopopulation writes these on every record made in VariantGrid, so copying them carried
# nothing a curator chose and made a significance-only record look like it had content to apply. A new record gets
# its own: search_terms for its variant, and the PubMed gene count is looked up live (cached).
_AUTOPOPULATED_COPY_SCOPES = {
    'search_terms': 'A',
    'pubmed_gene_search_count': 'G',
}


def _copy_scope_none(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_AUTOPOPULATED_COPY_SCOPES).update(copy_scope='N')


def _restore_copy_scope(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    for key, copy_scope in _AUTOPOPULATED_COPY_SCOPES.items():
        EvidenceKey.objects.filter(key=key).update(copy_scope=copy_scope)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0007_avi_splicevault_utr_annotator_ekeys'),
    ]

    operations = [
        migrations.RunPython(_copy_scope_none, reverse_code=_restore_copy_scope),
    ]
