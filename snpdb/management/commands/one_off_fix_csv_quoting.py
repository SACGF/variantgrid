import csv
import os
import zipfile

from django.core.management.base import BaseCommand
from snpdb.models import CachedGeneratedFile


class Command(BaseCommand):
    """

    """
    @staticmethod
    def fix_csv_zip(zip_filename):
        # Move the old one to a backup
        backup_name = zip_filename + ".backup"
        os.rename(zip_filename, backup_name)

        with zipfile.ZipFile(zip_filename, 'r') as zip_file:
            csv_filename = os.path.splitext(os.path.basename(zip_filename))[0]
            with zip_file.open(csv_filename) as csv_file:
                reader = csv.reader(csv_file.read().decode('utf-8').splitlines(), dialect='excel', escapechar='\\',
                                    quoting=csv.QUOTE_NONE)

            fixed_csv_filename = os.path.splitext(zip_filename)[0]
            with open(fixed_csv_filename, "wt") as out_csv_file:
                writer = csv.writer(out_csv_file, dialect='excel', quoting=csv.QUOTE_MINIMAL)
                for r in reader:
                    writer.writerow(r)

    def handle(self, *args, **options):
        AFFECTED_GENERATORS = ["export_sample_to_downloadable_file", "export_cohort_to_downloadable_file"]
        for cgf in CachedGeneratedFile.objects.filter(generator__in=AFFECTED_GENERATORS):
            self.fix_csv_zip(cgf.filename)

