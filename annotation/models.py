# This file exists for PyCharm's Django Structure plugin
# pylint: disable=unused-import
from annotation.models.models import ClinVarVersion, ClinVar, ClinVarRecordCollection, ClinVarRecord, \
    ClinVarCitationsCollection, ClinVarCitation, DBNSFPGeneAnnotationVersion, DBNSFPGeneAnnotation, \
    GeneAnnotationVersion, GeneAnnotation, HumanProteinAtlasAnnotationVersion, HumanProteinAtlasTissueSample, \
    HumanProteinAtlasAnnotation, ColumnVEPField, VariantAnnotationVersion, VCFAnnotationStats, AnnotationRangeLock, \
    AnnotationRun, VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap, ManualVariantEntryCollection, \
    ManualVariantEntry, CreatedManualVariant, AnnotationVersion, CachedWebResource, GeneSymbolCitation, \
    GenePubMedCount, MutationalSignatureInfo
from annotation.models.models_citations import Citation
from annotation.models.models_gene_counts import VariantSource, SampleAnnotationVersionVariantSource, GeneCountType, \
    GeneValue, GeneValueCountCollection, GeneValueCount, CohortGeneCounts
from annotation.models.models_phenotype_match import DescriptionProcessingStatus, PhenotypeDescription, TextPhenotype, \
    TextPhenotypeSentence, TextPhenotypeMatch, PatientTextPhenotype
from annotation.models.models_sample_stats import SampleVariantAnnotationStats, \
    SampleVariantAnnotationStatsPassingFilter, SampleGeneAnnotationStats, SampleGeneAnnotationStatsPassingFilter, \
    SampleClinVarAnnotationStats, SampleClinVarAnnotationStatsPassingFilter
from annotation.models.models_version_diff import VersionDiff, VersionDiffFromToResult, VersionDiffChangeCountResult, \
    VariantAnnotationVersionDiff
# pylint: enable=unused-import
