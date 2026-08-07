from rest_framework.status import HTTP_200_OK

from library.django_utils.rollbar_middleware import RollbarIgnoreException


class CeleryTasksObsoleteException(RollbarIgnoreException):
    """ Throw this in a celery task to kill dependent jobs, but not report to rollbar """


class NonFatalNodeError(Exception):
    status = HTTP_200_OK


class NodeParentErrorsException(NonFatalNodeError):

    def __str__(self):
        return "The node's parents had errors"


class NodeParentNotReadyException(NonFatalNodeError):

    def __str__(self):
        return "The node's parents were not ready"


class NodeConfigurationException(NonFatalNodeError):

    def __str__(self):
        return "The node has bad configuration"


class NodeOutOfDateException(NonFatalNodeError):

    def __str__(self):
        return "The node was updated while loading."


class NodeNotFoundException(NonFatalNodeError):

    def __init__(self, node_id):
        self.node_id = node_id
        super().__init__(f"Node {node_id} not found")

    def __str__(self):
        return "The node was deleted while loading"


class NodeOutOfMemoryException(Exception):
    """ A node load exceeded the worker memory cap (RLIMIT_AS -> MemoryError). The node is
        permanently failed BEFORE this is raised so it is never re-run; raising (a plain
        Exception, not RollbarIgnoreException) lets the celery task_failure handler report it to
        Rollbar with the analysis/node context. """

    def __init__(self, analysis_id, node_id, version):
        self.analysis_id = analysis_id
        self.node_id = node_id
        self.version = version
        super().__init__(f"Node {node_id}/{version} (analysis {analysis_id}) ran out of memory "
                         f"during load() and was permanently failed.")
