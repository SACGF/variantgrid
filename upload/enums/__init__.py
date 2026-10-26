"""
upload's enums as one namespace, outside the models package so code the models import can use them
without running upload/models/__init__. Callers write `from upload.enums import UploadedFileTypes, VCFPipelineStage`.
"""
from upload.enums.upload_enums import *
