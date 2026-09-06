"""
upload's models as one namespace (models.py: FileUpload, UploadPipeline, UploadStep and the
Uploaded* satellites; models_uploaded_files.py), so callers write `from upload.models import
UploadedVCF, UploadPipeline`.
"""
from .models import *
from .models_uploaded_files import *
