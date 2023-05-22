import inspect
from abc import ABC, abstractmethod
from typing import Iterable, Type, List

from django.conf import settings
from django.contrib.auth.models import User
from django.db import models

from library.utils import get_all_subclasses
from upload.models import UploadPipeline


class ImportTaskFactory(ABC):
    """ Subclass this to dispatch uploaded files to tasks """

    @abstractmethod
    def get_uploaded_file_type(self) -> str:
        pass

    @abstractmethod
    def get_data_classes(self) -> Iterable[Type[models.Model]]:
        """ e.g. return UploadedVCF, UploadedGeneList """
        pass

    @abstractmethod
    def get_possible_extensions(self) -> Iterable[str]:
        """ e.g. return ['csv', 'xls'] """
        pass

    def get_processing_ability(self, user: User, filename: str, file_extension: str) -> int:
        """ If you can't process EVERY file of type in extensions, overwrite this and check.
            > 0 means you can process it - the processor with the highest value will do it
        """
        return 1

    @abstractmethod
    def create_import_task(self, upload_pipeline: UploadPipeline):
        pass


def get_import_task_factories() -> List[ImportTaskFactory]:
    # Import all factory scripts into scope  so __subclasses__ works
    for i in settings.IMPORT_TASK_FACTORY_IMPORTS:
        exec(f"import {i}")

    factories = []
    for itf_class in get_all_subclasses(ImportTaskFactory):
        if not inspect.isabstract(itf_class):
            factories.append(itf_class())
#        else:
#            logging.debug("Warning: not looking at %s", itf_class)

    return factories
