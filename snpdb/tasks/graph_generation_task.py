import logging
import traceback

import celery

from eventlog.models import create_event
from library import file_utils
from library.enums.log_level import LogLevel
from library.log_utils import get_traceback
from library.utils import import_class


@celery.task(ignore_result=False)
def generate_graph(graph_class_name, *args):
    """
        cacheablegraph : snpdb.graphs.graphcache.CacheableGraph
        returns filename of generated graph
    """

    graph_klass = import_class(graph_class_name)
    cacheablegraph = graph_klass(*args)
    filename = cacheablegraph.get_filename()
    try:
        file_utils.mk_path_for_file(filename)
        cacheablegraph.save(filename)
    except Exception as e:
        logging.error("Problem generating graph %s", graph_klass)
        traceback.print_exc()

        details = get_traceback()
        create_event(None, "generate_graph", details=details, filename=str(graph_klass), severity=LogLevel.ERROR)

        raise e
    return filename
