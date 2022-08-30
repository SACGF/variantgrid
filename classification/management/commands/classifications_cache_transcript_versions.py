from django.core.management import BaseCommand

from classification.models import Classification


class Command(BaseCommand):

    def handle(self, *args, **options):
        for i, c in enumerate(Classification.objects.all()):
            if i % 100 == 0:
                print(f"Processed {i:,} classification")
            update_count = c.update_transcripts()
            if update_count != 2:
                print(f"For classification {c.pk} only set {update_count} of 2 transcript versions")
            c.save(update_fields=('transcript_version_grch37', 'transcript_version_grch38'))
        print("Done")
