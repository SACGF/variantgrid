import abc
import logging
import os

from celery.result import AsyncResult
from django.conf import settings
from django.utils import timezone

from library.graphs import graph_base
from library.utils import import_class
from snpdb.models import CachedGeneratedFile
from snpdb.tasks.graph_generation_task import generate_graph


class CacheableGraph(graph_base.GraphBase, metaclass=abc.ABCMeta):

    def __hash__(self):
        return hash((self.get_name(), self.get_param_hash()))

    @abc.abstractmethod
    def get_params_hash(self):
        """ This should be unique for every graph produced (from parameters) """
        pass

    def get_filename(self):
        return os.path.join(settings.GENERATED_DIR, self.get_name(), self.get_params_hash() + ".png")

    def get_name(self):
        return self.__class__.__name__


def async_graph(graph_class_name, *args):
    """ Return the output from a cacheablegraph, or go and generate it """

    graph_klass = import_class(graph_class_name)
    cacheablegraph = graph_klass(*args)
    generator = cacheablegraph.get_name()
    params_hash = cacheablegraph.get_params_hash()

    (cached_graph, created) = CachedGeneratedFile.objects.get_or_create(generator=generator,
                                                                        params_hash=params_hash)
    if created or not cached_graph.task_id:
        logging.debug("Launching Celery Job for graph: generator=%s, params_hash=%s", generator, params_hash)
        async_result = generate_graph.delay(graph_class_name, *args)  # @UndefinedVariable
        cached_graph.task_id = async_result.id
        cached_graph.generate_start = timezone.now()
        cached_graph.save()
    else:
        async_result = AsyncResult(cached_graph.task_id)

    if async_result.result:
        cached_graph.save_from_async_result(async_result)

    return cached_graph
