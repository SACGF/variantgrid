from rest_framework.status import HTTP_500_INTERNAL_SERVER_ERROR, HTTP_404_NOT_FOUND


class NonFatalNodeError(Exception):
    status = HTTP_500_INTERNAL_SERVER_ERROR


class NodeParentErrorsException(NonFatalNodeError):
    pass


class NodeParentNotReadyException(NonFatalNodeError):
    pass


class NodeConfigurationException(NonFatalNodeError):
    pass


class NodeNotFoundException(NonFatalNodeError):
    status = HTTP_404_NOT_FOUND

    def __init__(self, node_id):
        self.node_id = node_id
        super().__init__(f"Node {node_id} not found")
