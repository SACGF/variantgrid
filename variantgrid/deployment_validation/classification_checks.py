from classification.models import ClassificationModification
from classification.views.classification_export_report import ClassificationReport


def check_classification_reports() -> dict:
    classification_reports = {}

    if cm := ClassificationModification.objects.last():
        classification_report = ClassificationReport(cm, cm.user)
        unknown_evidence = classification_report.get_unknown_evidence()
        valid = not unknown_evidence
        if unknown_evidence:
            unknown = []
            colons = []
            for k in unknown_evidence:
                if ":" in k:
                    colons.append(k)
                else:
                    unknown.append(k)
            unknown = ",".join(unknown)
            colons = ",".join(colons)
            fixes = []
            if unknown:
                fixes.append(f"modify {unknown=}")
            if colons:
                fixes.append(f"change colon (':') to underscore ('_') for: {colons=}")
            fix = f"Modify ClassificationReportTemplate: " + " and ".join(fixes)
        else:
            fix = ""

        data = {
            "valid": valid,
            "fix": fix,
        }
        classification_reports["evidence_keys_valid"] = data
    return classification_reports
