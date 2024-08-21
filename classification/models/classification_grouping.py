# from django.db.models import CASCADE, TextChoices
# from django_extensions.db.models import TimeStampedModel
#
# from classification.enums import AlleleOriginBucket, ShareLevel
# from django.db import models
#
# from snpdb.models import Allele, Lab
#
#
# class ClassificationQualityLevel(TextChoices):
#     STANDARD = "S", "Standard"
#     LEGACY = "L", "Legacy"
#     INCOMPLETE = "I", "Incomplete"
#
#
# class ClassificationClassificationBucket(TextChoices):
#     BENIGN = "B", "Bengign"
#     VUS = "V", "VUS"
#     PATHOGENIC = "P", "Pathogenic"
#     OTHER = "O", "Other"
#     CONFLICTING = "C", "Conflicting"
#
#
# class ClassificationGrouping(TimeStampedModel):
#     # key
#     allele = models.ForeignKey(Allele, on_delete=models.CASCADE)
#     lab = models.ForeignKey(Lab, on_delete=CASCADE)
#     allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices)
#     share_level = models.CharField(max_length=16, choices=ShareLevel.choices())
#     quality_level = models.CharField(max_length=1, choices=ClassificationQualityLevel.choices())
#
#     # # for discordances
#     # classification_bucket
#     #
#     # # search details
#     # gene_symbols
#     # conditions
#     # curated
#     #
#     # # cached details
#     # summary_curated
#     # summary_classification
#     # summary_somatic_clinical_significance
#     # summary_c_hgvses
#     # summary_criteria
#
