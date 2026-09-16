"""
The case report - a lab's whole report design as one Django template, rendered server side.

Deliberately empty: classification/models/classification_report_models.py imports
template_validation from here, so anything this package imports from classification.models would
be a cycle. Import the module you want.

- case_report_context.py - ReportContext / ReportVariant, the ordering rules and amp_tier
- case_report_builder.py - building, rebuilding and finalising a CaseReport
- renderers.py - HTML, PDF, DOCX and JSON from one ReportContext
- template_validation.py - the fixture render a template has to survive to be saved
"""
