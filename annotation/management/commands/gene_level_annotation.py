from Bio import Entrez
from django.core.management import BaseCommand
from django.utils import timezone

from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from library.django_utils import chunked_queryset





class Command(BaseCommand):

    def handle(self, *args, **options):
        print(f"Started: {timezone.now()}")


        print(f"Finished: {timezone.now()}")
