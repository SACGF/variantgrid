import csv
import logging
import os
import zipfile

from django.core.management.base import BaseCommand

from snpdb.models import CachedGeneratedFile


class Command(BaseCommand):
    """

    """
    @staticmethod
    def fix_csv(csv_file, fixed_csv_filename):
        reader = csv.reader(csv_file.read().decode('utf-8').splitlines(), dialect='excel', escapechar='\\',
                            quoting=csv.QUOTE_NONE)

        logging.info("Fixing in file: %s", fixed_csv_filename)
        with open(fixed_csv_filename, "wt") as out_csv_file:
            writer = csv.writer(out_csv_file, dialect='excel', quoting=csv.QUOTE_MINIMAL)
            for r in reader:
                writer.writerow(r)

    @staticmethod
    def fix_csv_zip(zip_file, zip_filename):
        csv_filename = os.path.splitext(os.path.basename(zip_filename))[0]
        with zip_file.open(csv_filename) as csv_file:
            fixed_csv_filename = os.path.splitext(zip_filename)[0]
            Command.fix_csv(csv_file, fixed_csv_filename)

        zip_file_path = fixed_csv_filename + ".zip"
        logging.info("Putting into zip: %s", zip_file_path)
        with zipfile.ZipFile(zip_file_path, 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(fixed_csv_filename, arcname=os.path.basename(fixed_csv_filename))
        os.unlink(fixed_csv_filename)

    def handle(self, *args, **options):
        AFFECTED_GENERATORS = ["export_sample_to_downloadable_file", "export_cohort_to_downloadable_file"]
        for cgf in CachedGeneratedFile.objects.filter(generator__in=AFFECTED_GENERATORS, filename__isnull=False):
            if not (cgf.filename.endswith(".csv") or cgf.filename.endswith(".csv.zip")):
                continue

            if os.path.exists(cgf.filename):
                logging.info("Need to fix: %s", cgf.filename)
                # Move the old one to a backup
                backup_name = cgf.filename + ".backup"
                os.rename(cgf.filename, backup_name)
                if cgf.filename.endswith(".zip"):
                    with zipfile.ZipFile(backup_name, 'r') as zip_file:
                        self.fix_csv_zip(zip_file, cgf.filename)
                else:
                    with open(backup_name, "rb") as csv_file:
                        Command.fix_csv(csv_file, cgf.filename)
