from classification.models import ClassificationModification
from classification.views.classification_export_report import ClassificationReport


def check_classification_reports() -> dict:
    classification_reports = {}

    if cm := ClassificationModification.objects.last():
        classification_report = ClassificationReport(cm, cm.user)
        unknown_evidence = classification_report.get_unknown_evidence()
        valid = not unknown_evidence
        if unknown_evidence:
            unknown_evidence = ", ".join(sorted(unknown_evidence))
            fix = f"Modify ClassificationReportTemplate and fix {unknown_evidence=}"
        else:
            fix = ""

        data = {
            "valid": valid,
            "fix": fix,
        }
        classification_reports["evidence_keys_valid"] = data
    return classification_reports
