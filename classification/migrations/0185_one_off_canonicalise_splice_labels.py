from django.db import migrations

from genes.gene_splice import canonical_splice_label
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from manual.operations.manual_operations import ManualOperation

SPLICE_ALT_PREFIX = f"<{GeneLevelSymbolicAlt.SPLICE}:"


def _has_uncanonical_splice_variants(apps):
    """ A splice Variant loaded before #1835 holds the label as it was written (V7, EX14SKIP), which is
        no longer the label a classification for the same junction resolves to """
    Variant = apps.get_model("snpdb", "Variant")
    for alt_seq in Variant.objects.filter(alt__seq__startswith=SPLICE_ALT_PREFIX) \
                                  .values_list("alt__seq", flat=True):
        if parsed := GeneLevelSymbolicAlt.parse(alt_seq):
            label = parsed[3]
            if canonical_splice_label(label) != label:
                return True
    return False


class Migration(migrations.Migration):
    dependencies = [
        ("classification", "0184_one_off_revalidate_gene_level_allele_infos"),
        ("genes", "0095_canonical_splice_labels"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["splice_labels_canonicalise"]),
                        note="Re-point splice variants at the canonical label for their junction, so a "
                             "classification naming it reaches the variant we have (#1835). Writes "
                             "snpdb_variant - run 'manage.py splice_labels_canonicalise --dry-run' first",
                        test=_has_uncanonical_splice_variants),
    ]
