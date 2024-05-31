import logging
import os

from django.conf import settings
from django.core.exceptions import PermissionDenied
from django_messages.admin import User

from annotation.models import ClinVarReviewStatus, Variant
from annotation.models.models import ClinVar, ClinVarVersion
from annotation.vcf_files.vcf_types import VCFVariant
from snpdb.models import VariantCoordinate
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv

"""
##fileformat=VCFv4.1
##fileDate=2024-05-28
##source=ClinVar
##reference=GRCh38
##ID=<Description="ClinVar Variation ID">
##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">
##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="allele frequencies from ExAC">
##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="allele frequencies from TGP">
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier submitted for germline classifications, e.g. OMIM:NNNNNN">
##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier for germline classifications, e.g. OMIM:NNNNNN">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status of germline classification for the Variation ID">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Aggregate germline classification for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting germline classification for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Germline classification for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:classification; multiple values are separated by a vertical bar">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">
##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Sequence Ontology id for variant type">
##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
##INFO=<ID=DBVARID,Number=.,Type=String,Description="nsv accessions from dbVar for the variant">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">
##INFO=<ID=ONCDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in ONCDISDB">
##INFO=<ID=ONCDNINCL,Number=.,Type=String,Description="For included variant: ClinVar's preferred disease name for the concept specified by disease identifiers in ONCDISDBINCL">
##INFO=<ID=ONCDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier submitted for oncogenicity classifications, e.g. MedGen:NNNNNN">
##INFO=<ID=ONCDISDBINCL,Number=.,Type=String,Description="For included variant: Tag-value pairs of disease database name and identifier for oncogenicity classifications, e.g. OMIM:NNNNNN">
##INFO=<ID=ONC,Number=.,Type=String,Description="Aggregate oncogenicity classification for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=ONCINCL,Number=.,Type=String,Description="Oncogenicity classification for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:classification; multiple values are separated by a vertical bar">
##INFO=<ID=ONCREVSTAT,Number=.,Type=String,Description="ClinVar review status of oncogenicity classification for the Variation ID">
##INFO=<ID=ONCCONF,Number=.,Type=String,Description="Conflicting oncogenicity classification for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=SCIDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in SCIDISDB">
##INFO=<ID=SCIDNINCL,Number=.,Type=String,Description="For included variant: ClinVar's preferred disease name for the concept specified by disease identifiers in SCIDISDBINCL">
##INFO=<ID=SCIDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier submitted for somatic clinial impact classifications, e.g. MedGen:NNNNNN">
##INFO=<ID=SCIDISDBINCL,Number=.,Type=String,Description="For included variant: Tag-value pairs of disease database name and identifier for somatic clinical impact classifications, e.g. OMIM:NNNNNN">
##INFO=<ID=SCIREVSTAT,Number=.,Type=String,Description="ClinVar review status of somatic clinical impact for the Variation ID">
##INFO=<ID=SCI,Number=.,Type=String,Description="Aggregate somatic clinical impact for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=SCIINCL,Number=.,Type=String,Description="Somatic clinical impact classification for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:classification; multiple values are separated by a vertical bar">
"""


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

        # new oncogenic fields
        'ONCREVSTAT': 'oncogenic_review_status',
        'ONCINCL': 'oncogenic_classification',
        'ONCCONF': 'oncogenic_conflicting_classification',
        'ONCDN': 'oncogenic_preferred_disease_name',
        'ONCDISDB': 'oncogenic_disease_database_name',

        # new somatic fields
        'SCIREVSTAT': 'somatic_review_status',
        'SCI': 'somatic_clinical_significance',

        'SCIDN': 'somatic_preferred_disease_name',
        'SCIDISDB': 'somatic_disease_database_name'
        # ##INFO=<ID=SCIDISDBINCL,Number=.,Type=String,Description="For included variant: Tag-value pairs of disease database name and identifier for somatic clinical impact classifications, e.g. OMIM:NNNNNN">
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
                          'drug_response',

                          'oncogenic_review_status',
                          'oncogenic_classification',
                          'oncogenic_conflicting_classification',
                          'oncogenic_preferred_disease_name',
                          'oncogenic_disease_database_name',

                          'somatic_review_status',
                          'somatic_clinical_significance',
                          'somatic_preferred_disease_name',
                          'somatic_disease_database_name'
                          ]

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
        self.variant_pk_lookup = VariantPKLookup(clinvar_version.genome_build)
        review_status_vcf_mappings_dict = dict(ClinVarReviewStatus.VCF_MAPPINGS)
        self.field_formatters = {
            "clinvar_review_status": lambda x: review_status_vcf_mappings_dict[x],
            "somatic_review_status": lambda x: review_status_vcf_mappings_dict[x],
            "oncogenic_review_status": lambda x: review_status_vcf_mappings_dict[x]
        }

    def process_variant(self, v: VCFVariant):
        alt = self.variant_single_alt(v)
        variant_coordinate = VariantCoordinate.from_explicit_no_svlen(v.CHROM, v.POS, v.REF, alt)
        variant_hash = self.variant_pk_lookup.add(variant_coordinate)
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

    def create_clinvar_for_variant_id(self, variant_id: int, v: VCFVariant):
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
            # FIXME do the same for somatic classification
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
                # CLNSIGCONF (conflicting_clinical_significance), e.g. see
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


def check_can_import_clinvar(user: User):
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
