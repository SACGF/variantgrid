import logging

from library.guardian_utils import assign_permission_to_user_and_groups


# This method is called to release samples when annotated (either if imported with all known variants, or after annotation)
def set_vcf_and_samples_import_status(vcf, import_status):
    logging.info("setting vcf %s (%d) and sample to %s", vcf, vcf.pk, import_status)
    vcf.import_status = import_status
    vcf.save()

    # Update related vcfs with same status
    vcf.sample_set.update(import_status=import_status)

    try:
        cohort = vcf.cohort
        cohort.import_status = import_status
        cohort.save()
        assign_permission_to_user_and_groups(vcf.user, cohort)
        logging.info("Set cohort %s to %s", cohort, import_status)
    except:
        pass
