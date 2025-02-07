# This file exists for PyCharm's Django Structure plugin
# pylint: disable=unused-import
from snpdb.models.models import Tag, CachedGeneratedFile, Company, Manufacturer, Wiki, ImportedWikiCollection, \
    ImportedWiki, Organization, ClinVarKey, ClinVarKeyExcludePattern, Country, State, Lab, LabHead, UserAward, \
    LabProject, SiteMessage
from snpdb.models.models_clingen_allele import ClinGenAllele
from snpdb.models.models_cohort import Cohort, CohortSample, CohortGenotypeTaskVersion, \
    CohortGenotypeCommonFilterVersion, CommonVariantClassified, CohortGenotypeCollection, CohortGenotype, Trio
from snpdb.models.models_columns import VariantGridColumn, ColumnVCFInfo, CustomColumnsCollection, CustomColumn
from snpdb.models.models_dbsnp import DbSNP
from snpdb.models.models_genome import GenomeBuild, GenomeBuildPatchVersion, Contig, GenomeBuildContig, GenomeFasta, \
    GenomeFastaContig
from snpdb.models.models_genomic_interval import GenomicIntervalsCategory, GenomicIntervalsCollection, GenomicInterval
from snpdb.models.models_somalier import SomalierVCFExtract, SomalierSampleExtract, SomalierAncestryRun, \
    SomalierAncestry, SomalierCohortRelate, SomalierTrioRelate, SomalierAllSamplesRelate, SomalierRelatePairs
from snpdb.models.models_user_settings import UserDataPrefix, TagColorsCollection, TagColor, UserPageAck, \
    UserGridConfig, SettingsOverride, GlobalSettings, OrganizationUserSettingsOverride, LabUserSettingsOverride, \
    UserSettingsOverride, SettingsInitialGroupPermission, NodeCountSettingsCollection, NodeCountSettings, UserContact
from snpdb.models.models_variant import Allele, AlleleMergeLog, Sequence, Locus, Variant, VariantWiki, VariantAllele, \
    VariantCollection, VariantCollectionRecord, AlleleSource, VariantAlleleSource, VariantAlleleCollectionSource, \
    VariantAlleleCollectionRecord, LiftoverRun, AlleleLiftover
from snpdb.models.models_vcf import Project, VCF, VCFFilter, VCFTag, Sample, SampleFilePath, SampleTag, VCFAlleleSource, \
    SampleStatsCodeVersion, SampleStats, SampleStatsPassingFilter, SampleLocusCount, SampleLabProject, \
    VCFSourceSettings, VCFBedIntersection
from snpdb.models.models_zygosity_counts import VariantZygosityCountCollection, VariantZygosityCount, \
    VariantZygosityCountForVCF, VariantZygosityCountForSample
# pylint: enable=unused-import