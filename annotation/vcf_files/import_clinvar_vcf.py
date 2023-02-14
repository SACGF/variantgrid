import logging
import os

from django.conf import settings
from django.core.exceptions import PermissionDenied

from annotation.models import ClinVarReviewStatus, Variant
from annotation.models.models import ClinVar, ClinVarVersion
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


class BulkClinVarInserter:
    CLINVAR_INFO_MAPPINGS = {
        'ALLELEID': 'clinvar_allele_id',
        'CLNDN': 'clinvar_preferred_disease_name',
        'CLNDISDB': 'clinvar_disease_database_name',
        'CLNREVSTAT': 'clinvar_review_status',
        'CLNSIG': 'clinical_significance',
        'CLNSIGCONF': 'conflicting_clinical_significance',
        'CLNVI': 'clinvar_clinical_sources',
        'ORIGIN': 'clinvar_origin',
        'SSR': 'clinvar_suspect_reason_code',
    }

    CLINVAR_MANDATORY_FIELDS = ['clinvar_variation_id', 'clinvar_allele_id']

    CLINVAR_CSV_FIELDS = ['version_id',
                          'variant_id',
                          'clinvar_variation_id',
                          'clinvar_allele_id',
                          'clinvar_preferred_disease_name',
                          'clinvar_disease_database_name',
                          'clinvar_review_status',
                          'clinical_significance',
                          'conflicting_clinical_significance',
                          'highest_pathogenicity',
                          'clinvar_clinical_sources',
                          'clinvar_origin',
                          'clinvar_suspect_reason_code',
                          'drug_response']

    CLINSIG_TO_PATHOGENICITY = {
        "Benign": 1,
        "Likely_benign": 2,
        "Uncertain_significance": 3,
        "Likely_pathogenic": 4,
        "Pathogenic": 5,
    }

    MAX_CONFLICTING_RECORDS_MISSING_CLINSIGCONF = 1000

    def __init__(self, clinvar_version, upload_step):
        self.clinvar_version = clinvar_version
        self.variant_by_variant_hash = {}
        self.upload_step = upload_step
        self.items_processed = 0
        self.batch_id = 0
        self.conflicting_missing_clinsigconf = 0
        self.variant_pk_lookup = VariantPKLookup.factory(clinvar_version.genome_build)
        review_status_vcf_mappings_dict = dict(ClinVarReviewStatus.VCF_MAPPINGS)
        self.field_formatters = {
            "clinvar_review_status": lambda x: review_status_vcf_mappings_dict[x]
        }

    def process_variant(self, v):
        alt = self.variant_single_alt(v)
        variant_hash = self.variant_pk_lookup.add(v.CHROM, v.POS, v.REF, alt)
        self.variant_by_variant_hash[variant_hash] = v
        if len(self.variant_by_variant_hash) >= settings.SQL_BATCH_INSERT_SIZE:
            self.bulk_insert()

    def finish(self):
        self.bulk_insert()

    def bulk_insert(self):
        print("BulkClinVarInserter bulk_insert")
        clinvar_list = []
        self.variant_pk_lookup.batch_check()
        for variant_hash, variant_pk in self.variant_pk_lookup.variant_pk_by_hash.items():
            v = self.variant_by_variant_hash[variant_hash]
            try:
                cv = self.create_clinvar_for_variant_id(variant_pk, v)
            except:
                logging.error("Could not process ClinVar variant: '%s'", v)
                raise
            clinvar_list.append(cv)

        # I messed up the alt for reference being '=' vs None and found it with this check. Could possibly remove...
        if num_unknown := len(self.variant_pk_lookup.unknown_variant_coordinates):
            unknown_example = self.variant_pk_lookup.unknown_variant_coordinates[0]
            msg = f"BulkClinVarInserter, has {num_unknown} unknowns! 1st is: {unknown_example}"
            raise ValueError(msg)

        self.variant_by_variant_hash = {}
        self.variant_pk_lookup.clear()

        if clinvar_list:
            def create_clinvar_tuple(cv):
                values_list = []
                for f in BulkClinVarInserter.CLINVAR_CSV_FIELDS:
                    values_list.append(getattr(cv, f))
                return tuple(values_list)

            DELIMITER = '\t'
            PREFIX = "clinvar"
            processing_dir = self.upload_step.upload_pipeline.get_pipeline_processing_subdir(PREFIX)
            basename = f"{PREFIX}_step_{self.upload_step.pk}_batch_{self.batch_id}.tsv"
            csv_filename = os.path.join(processing_dir, basename)
            logging.info("Writing ClinVar copy CSV")
            clinvar_tuples = map(create_clinvar_tuple, clinvar_list)
            write_sql_copy_csv(clinvar_tuples, csv_filename, delimiter=DELIMITER)
            partition_table = self.clinvar_version.get_partition_table()
            logging.info("Inserting file '%s' into partition %s", csv_filename, partition_table)
            sql_copy_csv(csv_filename, partition_table, BulkClinVarInserter.CLINVAR_CSV_FIELDS, delimiter=DELIMITER)

            self.batch_id += 1
            self.items_processed += len(clinvar_list)

    def create_clinvar_for_variant_id(self, variant_id, v):
        kwargs = {"variant_id": variant_id,
                  "version": self.clinvar_version,
                  "clinvar_variation_id": v.ID}

        for info, field in BulkClinVarInserter.CLINVAR_INFO_MAPPINGS.items():
            value = v.INFO.get(info)
            if value is not None:
                formatter = self.field_formatters.get(field)
                if formatter:
                    value = formatter(value)
                kwargs[field] = value

        for field in BulkClinVarInserter.CLINVAR_MANDATORY_FIELDS:
            if kwargs.get(field) is None:
                raise ValueError(f"Mandatory field '{field}' missing: {v}")

        # clinical_significance is now a '|' separated string
        if clinical_significance := kwargs.get("clinical_significance"):
            drug_response = "drug_response" in clinical_significance
            highest_pathogenicity = 0
            if clinical_significance.startswith("Conflicting_interpretations_of_pathogenicity"):
                multiple_clinical_significances = kwargs.get("conflicting_clinical_significance")
            else:
                multiple_clinical_significances = clinical_significance

            if multiple_clinical_significances:
                for clnsig, pathogenicity in BulkClinVarInserter.CLINSIG_TO_PATHOGENICITY.items():  # low->high
                    if clnsig in multiple_clinical_significances:
                        highest_pathogenicity = pathogenicity
            else:
                # 3 out of 50504 records in 20210828 with Conflicting_interpretations_of_pathogenicity are missing
                # CLNSIGCONF (conflicting_clinical_significance), eg see
                # https://www.ncbi.nlm.nih.gov/clinvar/variation/161486/
                # As this is really low we'll just skip setting highest pathogenicity in those ones, but will die
                # If there are a few as that could mean ClinVar has changed the INFO fields
                self.conflicting_missing_clinsigconf += 1
                if self.conflicting_missing_clinsigconf > self.MAX_CONFLICTING_RECORDS_MISSING_CLINSIGCONF:
                    message = f"{self.MAX_CONFLICTING_RECORDS_MISSING_CLINSIGCONF} records had " \
                              f"CLINSIG=Conflicting_interpretations_of_pathogenicity but no CLNSIGCONF. " \
                              f"A few missing are expected but this many is likely due to INFO field changes"
                    raise ValueError(message)

            kwargs["highest_pathogenicity"] = highest_pathogenicity
            kwargs["drug_response"] = drug_response

        return ClinVar(**kwargs)

    @staticmethod
    def variant_single_alt(vcf_record):
        """ For reading decomposed VCFs """

        # ClinVar has reference entries with alt='.' - comes back as nothing
        num_alt = len(vcf_record.ALT)
        if num_alt:
            if num_alt == 1:
                return vcf_record.ALT[0]
            else:
                raise ValueError(f"Variant {vcf_record} has more than 1 alt: '{vcf_record.ALT}'!")
        return Variant.REFERENCE_ALT  # ALT[] is "." ie reference


def check_can_import_clinvar(user):
    if not user.is_staff:
        msg = "Only Authorised users can upload a ClinVar VCF"
        raise PermissionDenied(msg)


def import_clinvar_vcf(upload_step):
    """ This can run in parallel """
    logging.debug("import_clinvar_file start")

    clinvar_version = ClinVarVersion.objects.get(md5_hash=upload_step.uploaded_file.md5_hash)

    import cyvcf2
    vcf_reader = cyvcf2.VCF(upload_step.input_filename)
    bulk_inserter = BulkClinVarInserter(clinvar_version, upload_step)

    for v in vcf_reader:
        bulk_inserter.process_variant(v)

    bulk_inserter.finish()

    return bulk_inserter.items_processed
