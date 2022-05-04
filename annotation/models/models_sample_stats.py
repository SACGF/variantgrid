"""
    SampleStats: we want to store more info similar to snpdb.models.SampleStats however we need to store the annotation versions,
    as as it might change and become obsolete.

    As we want to keep the Annotation->Snpdb dependency one way only, we have to keep annotation sample stats here...

    @see snpdb.models.AbstractSampleStats

"""
from django.db import models
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel

from annotation.models.models import VariantAnnotationVersion, ClinVarVersion, GeneAnnotationVersion
from snpdb.models import Sample, SampleStatsCodeVersion
from snpdb.models.models_enums import BuiltInFilters


class AbstractVariantAnnotationStats(TimeStampedModel):
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)

    variant_dbsnp_count = models.IntegerField(default=0)
    insertions_dbsnp_count = models.IntegerField(default=0)
    snp_dbsnp_count = models.IntegerField(default=0)
    deletions_dbsnp_count = models.IntegerField(default=0)

    # Separated into zygosity as it's used by SampleNode counts
    ref_high_or_moderate_count = models.IntegerField(default=0)
    het_high_or_moderate_count = models.IntegerField(default=0)
    hom_high_or_moderate_count = models.IntegerField(default=0)
    unk_high_or_moderate_count = models.IntegerField(default=0)

    class Meta:
        abstract = True

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label):
        LABEL_COUNTS = {
            BuiltInFilters.IMPACT_HIGH_OR_MODERATE: {'ref': self.ref_high_or_moderate_count,
                                                     'het': self.het_high_or_moderate_count,
                                                     'hom': self.hom_high_or_moderate_count,
                                                     'unk': self.unk_high_or_moderate_count},
        }
        label_counts = LABEL_COUNTS[label]

        count = 0
        if zygosity_ref:
            count += label_counts['ref']

        if zygosity_het:
            count += label_counts['het']
        if zygosity_hom:
            count += label_counts['hom']
        if zygosity_unk:
            count += label_counts['unk']

        return count


class AbstractGeneAnnotationStats(TimeStampedModel):
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    gene_annotation_version = models.ForeignKey(GeneAnnotationVersion, on_delete=CASCADE)
    gene_count = models.IntegerField(default=0)

    # Separated into zygosity as it's used by SampleNode counts
    ref_omim_phenotype_count = models.IntegerField(default=0)
    hom_omim_phenotype_count = models.IntegerField(default=0)
    het_omim_phenotype_count = models.IntegerField(default=0)
    unk_omim_phenotype_count = models.IntegerField(default=0)

    class Meta:
        abstract = True

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        count = 0
        if zygosity_ref:
            count += self.ref_omim_phenotype_count
        if zygosity_het:
            count += self.het_omim_phenotype_count
        if zygosity_hom:
            count += self.hom_omim_phenotype_count
        if zygosity_unk:
            count += self.unk_omim_phenotype_count
        return count


class AbstractClinVarAnnotationStats(TimeStampedModel):
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    clinvar_version = models.ForeignKey(ClinVarVersion, on_delete=CASCADE)
    clinvar_count = models.IntegerField(default=0)

    # Separated into zygosity as it's used by SampleNode counts
    ref_clinvar_pathogenic_count = models.IntegerField(default=0)
    hom_clinvar_pathogenic_count = models.IntegerField(default=0)
    het_clinvar_pathogenic_count = models.IntegerField(default=0)
    unk_clinvar_pathogenic_count = models.IntegerField(default=0)

    class Meta:
        abstract = True

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        count = 0
        if zygosity_ref:
            count += self.ref_clinvar_pathogenic_count
        if zygosity_het:
            count += self.het_clinvar_pathogenic_count
        if zygosity_hom:
            count += self.hom_clinvar_pathogenic_count
        if zygosity_unk:
            count += self.unk_clinvar_pathogenic_count
        return count


class AbstractPassingFilter(models.Model):
    class Meta:
        abstract = True


class AbstractSampleVariantAnnotationStats(AbstractVariantAnnotationStats):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)

    class Meta:
        abstract = True

    @classmethod
    def load_version(cls, sample, annotation_version):
        return cls.objects.get(sample=sample, variant_annotation_version=annotation_version.variant_annotation_version)


class AbstractSampleGeneAnnotationStats(AbstractGeneAnnotationStats):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)

    class Meta:
        abstract = True

    @classmethod
    def load_version(cls, sample, annotation_version):
        return cls.objects.get(sample=sample, gene_annotation_version=annotation_version.gene_annotation_version)


class AbstractSampleClinVarAnnotationStats(AbstractClinVarAnnotationStats):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)

    class Meta:
        abstract = True

    @classmethod
    def load_version(cls, sample, annotation_version):
        return cls.objects.get(sample=sample, clinvar_version=annotation_version.clinvar_version)


class SampleVariantAnnotationStats(AbstractSampleVariantAnnotationStats):
    pass


class SampleVariantAnnotationStatsPassingFilter(AbstractSampleVariantAnnotationStats, AbstractPassingFilter):
    pass


class SampleGeneAnnotationStats(AbstractSampleGeneAnnotationStats):
    pass


class SampleGeneAnnotationStatsPassingFilter(AbstractSampleGeneAnnotationStats, AbstractPassingFilter):
    pass


class SampleClinVarAnnotationStats(AbstractSampleClinVarAnnotationStats):
    pass


class SampleClinVarAnnotationStatsPassingFilter(AbstractSampleClinVarAnnotationStats, AbstractPassingFilter):
    pass
