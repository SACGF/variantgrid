from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_evidence_keys(apps):
    # Only surface on a real deployment - a fresh/empty DB has nothing to forward-port.
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    return EvidenceKey.objects.exists()


class Migration(migrations.Migration):
    # legacy_somatic -direction forwardport loads the authoritative VG4 somatic EvidenceKey
    # config (amp:level_*, horak:*, assertion_method, clinical_significance,
    # somatic:clinical_significance, sample_type, somatic:tumor_cellularity) with their
    # point-based criteria properties. It writes namespace_overrides (added in 0125) and must
    # run after the migrations that install these keys (0118/0119/0120/0132) so its values are
    # authoritative - hence anchoring to the classification leaf. It is a ManualOperation because
    # the command is run once, out of band, during the VG4 upgrade (surfaced by upgrade.sh).
    dependencies = [
        ("classification", "0165_pathogenicity_ekeys_raw_scores_and_cnv"),
    ]

    operations = [
        ManualOperation.operation_manage(
            ["legacy_somatic", "-direction", "forwardport"],
            note="Run once when upgrading to VG4: loads the somatic EvidenceKey configuration.",
            test=_has_evidence_keys,
        ),
    ]
