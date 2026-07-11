from django.template import Library

from snpdb.models import VCF
from upload.models import UploadedVCF, SimpleVCFImportInfo

register = Library()


@register.inclusion_tag("upload/tags/vcf_import_info.html")
def vcf_import_info(vcf: VCF, examine_error_instructions):
    errors = []
    warnings = []
    try:
        upload_pipeline = vcf.uploadedvcf.uploaded_file.uploadpipeline
        for vii in upload_pipeline.get_errors():
            errors.append(vii)

        for vii in upload_pipeline.get_warnings(include_vcf=False):
            # Skipped annotation is surfaced via the dedicated banner + tab (issue #1409)
            if getattr(vii, "type", None) == SimpleVCFImportInfo.ANNOTATION_SKIPPED:
                continue
            warnings.append(vii)
    except UploadedVCF.DoesNotExist:
        msg = 'Import Error: This file has no associated import information. Please re-upload.'
        errors.append(msg)

    return {
        "errors": errors,
        "warnings": warnings,
        "examine_error_instructions": examine_error_instructions
    }
