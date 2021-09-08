import datetime

import pandas as pd
from django.core.management.base import BaseCommand

from library.log_utils import console_logger
from snpdb.models import Lab, LabProject, Organization

LEADER = 'Leader'
MEMBERS = 'Members'
NAME = 'Name'
INSTITUTION = 'Institution'
CITY = 'City'
COUNTRY = 'Country'
LAT = 'Lat'
LONG = 'Long'
FAMILIES = 'Families'
INVOLVED = 'Involved'
URL = 'URL'


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('LabNameLocations', help='csv file for Lab Details')

    def handle(self, *args, **options):
        filename = options["LabNameLocations"]
        logger = console_logger()

        df = pd.read_csv(filename, sep='\t', index_col=None)
        for col in [LEADER, MEMBERS, NAME, INSTITUTION, CITY, COUNTRY, LAT, LONG, FAMILIES, INVOLVED, URL]:
            if col not in df.columns:
                msg = f"Expected column '{col}' in tab separated file LabNameLocations"
                raise ValueError(msg)

        # Insert Lab Information and Project Data
        for _, row in df.iterrows():
            organization, _ = Organization.objects.get_or_create(name=row[INSTITUTION])

            lab = Lab.objects.create(name=row[NAME],
                                     organization=organization,
                                     city=row[CITY],
                                     country=row[COUNTRY],
                                     url=row[URL],
                                     lat=row[LAT],
                                     long=row[LONG])

            LabProject.objects.create(lab=lab,
                                      leader=row[LEADER],
                                      members=row[MEMBERS],
                                      families=row[FAMILIES],
                                      involved=row[INVOLVED] == "Y",
                                      date=datetime.date.today())

            print("saved lab '%s'" % row[NAME])
        logger.info("saved data")
