import importlib
import inspect
import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Iterable
from operator import itemgetter

from django.conf import settings
from django.contrib.auth.models import User
from django.db import models

from library.utils import get_all_subclasses
from upload.models import UploadPipeline


class ImportTaskFactory(ABC):
    """ Subclass this to dispatch uploaded files to tasks """

    @property
    def enabled(self) -> bool:
        """ False withdraws this file type from upload and from the API capabilities endpoint """
        return True

    @abstractmethod
    def get_uploaded_file_type(self) -> str:
        pass

    @abstractmethod
    def get_data_classes(self) -> Iterable[type[models.Model]]:
        """ e.g. return UploadedVCF, UploadedGeneList """
        pass

    @abstractmethod
    def get_possible_extensions(self) -> Iterable[str]:
        """ e.g. return ['csv', 'xls'] """
        pass

    def get_metadata_keys(self) -> frozenset[str]:
        """ Upload metadata keys this file type accepts - anything else is rejected at upload time.
            @see upload.upload_metadata """
        return frozenset()

    def get_processing_ability(self, user: User, filename: str, file_extension: str) -> int:
        """ If you can't process EVERY file of type in extensions, overwrite this and check.
            > 0 means you can process it - the processor with the highest value will do it
        """
        return 1

    @abstractmethod
    def create_import_task(self, upload_pipeline: UploadPipeline):
        pass


def get_import_task_factories() -> list[ImportTaskFactory]:
    # Import all factory scripts into scope  so __subclasses__ works
    for i in settings.IMPORT_TASK_FACTORY_IMPORTS:
        importlib.import_module(i)

    factories = []
    for itf_class in get_all_subclasses(ImportTaskFactory):
        if not inspect.isabstract(itf_class):
            itf = itf_class()
            if itf.enabled:
                factories.append(itf)
#        else:
#            logging.debug("Warning: not looking at %s", itf_class)

    return factories


def get_import_tasks_by_extension():
    possible_tasks = defaultdict(list)
    for itf in get_import_task_factories():
        for ext in itf.get_possible_extensions():
            possible_tasks[ext].append(itf)
    return possible_tasks


def get_import_task_factory_from_extension(user, filename, file_extension):
    possible_tasks = get_import_tasks_by_extension()
    possible_for_extension = possible_tasks[file_extension]

    tasks = []
    for possible in possible_for_extension:
        processing_ability = possible.get_processing_ability(user, filename, file_extension)
        if processing_ability:
            tasks.append((int(processing_ability), possible))

    if tasks:
        logging.debug("tasks: %s", tasks)
        tasks = sorted(tasks, key=itemgetter(0), reverse=True)
        last_pa = None
        for pa, _ in tasks:
            if last_pa is not None:
                if pa == last_pa:
                    logging.warning("Task for extension %s had 2 processors with equal ability - can't decide!", file_extension)
                    return None
            else:
                last_pa = pa

        return tasks[0][1]

    logging.warning("No tasks found for %s", tasks)
    return None
