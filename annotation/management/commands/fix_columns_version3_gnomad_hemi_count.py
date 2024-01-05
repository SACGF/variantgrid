from django.core.management.base import BaseCommand
from django.db.models import F

from annotation.models import VariantAnnotationVersion, VariantAnnotation


class Command(BaseCommand):
    def handle(self, *args, **options):
        for vav in VariantAnnotationVersion.objects.filter(gnomad='4.0', columns_version=3):
            print(f"Updating {vav} with gnomad_hemi_count")
            qs = VariantAnnotation.objects.filter(version=vav,
                                                  gnomad_xy_ac__isnull=False,
                                                  gnomad_non_par=True)
            qs.update(gnomad_hemi_count=F("gnomad_xy_ac"))
