from django.core.management.base import BaseCommand
import pandas as pd
from seqauto.models import SoftwarePipeline, SoftwarePipelineNode, SoftwarePipelineEdge
from library.log_utils import console_logger


PIPELINE_NAME = "Pipeline Name"
PIPELINE_VERSION = "Pipeline Version"
DESCRIPTION = "Description"
PARENT_NODE = "Parent Node"
NAME = "Name"
SOFTWARE_VERSION = "Software Version"
PARAMETERS = "Parameters"
SOFTWARE_DESCRIPTION = "Software Description"
NODE_PK = "Node Primary Key"
PIPELINE_PK = "Pipeline Primary Key"


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('SoftwarePipeline', help='csv file for Software Pipeline Details')

    def handle(self, *args, **options):
        filename = options["SoftwarePipeline"]

        logger = console_logger()

        df = pd.read_csv(filename, sep='\t', index_col=None)
        for col in [PIPELINE_NAME, PIPELINE_VERSION, DESCRIPTION, PARENT_NODE, NAME, SOFTWARE_VERSION, PARAMETERS, SOFTWARE_DESCRIPTION]:
            if col not in df.columns:
                msg = f"Expected column '{col}' in tab separated file SoftwarePipeline"
                raise ValueError(msg)

        logger.info("Loaded df")

        # Insert Software Pipeline Details
        prev_pipeline = ""
        sw_pipeline_node_pk = pd.Series([])
        sw_pipeline_pk = pd.Series([])
        for i, (_, row) in enumerate(df.iterrows()):
            # First create software pipeline record using version and name
            # possible there is more than one pipeline in the import
            if prev_pipeline != row[PIPELINE_NAME]:
                # Can only allow one pipeline with the same name and version
                sw_pipeline, created = SoftwarePipeline.objects.get_or_create(name=row[PIPELINE_NAME],
                                                                              version=row[PIPELINE_VERSION])
                prev_pipeline = row[PIPELINE_NAME]
                if not created:
                    print(f"Warning pipeline '{row[PIPELINE_NAME]}' version '{row[PIPELINE_VERSION]}' already exists")

            sw_pipeline_nd = SoftwarePipelineNode.objects.create(name=row[NAME],  # @UndefinedVariable
                                                                 version=row[SOFTWARE_VERSION],
                                                                 parameters=row[PARAMETERS],
                                                                 description=row[SOFTWARE_DESCRIPTION],
                                                                 softwarepipeline=sw_pipeline)
            sw_pipeline_node_pk[i] = sw_pipeline_nd.pk
            sw_pipeline_pk[i] = sw_pipeline.pk
            print(f"saved details for pipeline '{row[PIPELINE_NAME]}' and software '{row[SOFTWARE_VERSION]}' ")

        df.insert(8, "Node Primary Key", sw_pipeline_node_pk)
        df.insert(9, "Pipeline Primary Key", sw_pipeline_pk)

        #drop pipelines that were not imported -
        df = df.dropna(subset=['Pipeline Primary Key'])

        # Insert Parent Child Relationships
        for (_, row) in df.iterrows():
            if row.notnull()[PARENT_NODE]:
                try:
                    parent_node_id = SoftwarePipelineNode.objects.get(name=row[PARENT_NODE],  # @UndefinedVariable
                                                                      softwarepipeline=row[PIPELINE_PK])  # @UndefinedVariable
                    child_node_id = SoftwarePipelineNode.objects.get(pk=row[NODE_PK])  # @UndefinedVariable
                    sw_pipeline_edge = SoftwarePipelineEdge(child=child_node_id,
                                                            parent=parent_node_id)
                    sw_pipeline_edge.save()
                    print(f"saved parent child relationship between '{row[PARENT_NODE]}' and '{row[NAME]}' ")
                except:
                    print(f"parent child relationship not saved for row with Parent '{row[PARENT_NODE]}' with Pipeline ID '{row[PIPELINE_PK]}'")

        logger.info("saved data")
