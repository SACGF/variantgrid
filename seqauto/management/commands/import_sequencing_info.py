import pandas as pd
from django.core.management.base import BaseCommand

from library.log_utils import console_logger
from seqauto.models import Sequencer, SequencingInfo, EnrichmentKit
from snpdb.models import Lab, LabProject

LAB_NAME = "lab name"
INSTITUTION = "institution"
DOI = "doi"
PAPER_NAME = "paper name"
YEAR_PUB = "year pub"
ENRICHMENT_KIT = "enrichment kit"
SEQUENCER = "sequencer"
SEQ_DETAILS = "seq details"
FILE_TYPE = "file type"
FILE_COUNT = "file count"


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('SequencingInfo', help='csv file for Sequencing Details')

    def handle(self, *args, **options):
        filename = options["SequencingInfo"]

        logger = console_logger()

        df = pd.read_csv(filename, sep='\t', index_col=None)
        for col in [LAB_NAME, INSTITUTION, DOI, PAPER_NAME, YEAR_PUB, ENRICHMENT_KIT, SEQUENCER, SEQ_DETAILS, FILE_TYPE, FILE_COUNT]:
            if col not in df.columns:
                msg = f"Expected column '{col}' in tab separated file SequencingInfo"
                raise ValueError(msg)

        logger.info("Loaded df")

        # Insert Sequencing Details
        for _, row in df.iterrows():
            # First get LabProject ID using Lab Name and Institution
            try:
                lab = Lab.objects.get(name=row[LAB_NAME])
                project = LabProject.objects.get(lab=lab.id)

                #get enrichment_kit id
                enrichmentkit = EnrichmentKit.objects.get(name=row[ENRICHMENT_KIT])

                #get sequencer id
                sequencer = Sequencer.objects.get(name=row[SEQUENCER])

                if sequencer.name is None:
                    print("no sequencer in database with name %s" % row[SEQUENCER])

                    SequencingInfo.objects.create(lab_project=project,
                                                  doi=row[DOI],
                                                  paper_name=row[PAPER_NAME],
                                                  year_published=row[YEAR_PUB],
                                                  enrichment_kit=enrichmentkit,
                                                  sequencer=sequencer,
                                                  seq_details=row[SEQ_DETAILS],
                                                  file_type=row[FILE_TYPE],
                                                  file_count=row[FILE_COUNT])
                    print("saved sequence info for lab '%s'" % row[LAB_NAME])
            except:
                print(f"lab '{row[LAB_NAME]}' does not exist in the database")

        logger.info("saved data")
