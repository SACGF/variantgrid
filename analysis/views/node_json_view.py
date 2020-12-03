from django.http.response import JsonResponse
from django.views.generic.base import View
from rest_framework.status import HTTP_200_OK

from analysis.exceptions import NonFatalNodeError, NodeNotFoundException
from analysis.models.nodes.analysis_node import AnalysisNode
from library.log_utils import log_traceback


class NodeJSONViewMixin(View):

    def _get_data(self, request, *args, **kwargs):
        return NotImplementedError()

    def json_response(self, request, *args, **kwargs):
        try:
            data = self._get_data(request, *args, **kwargs)
            status = HTTP_200_OK
        except NonFatalNodeError as e:
            log_traceback()
            data = {"message": str(e),
                    "non_fatal": True}
            status = e.status

            if isinstance(e, NodeNotFoundException):
                # Is the node deleted, or perhaps it was just out of version?
                if e.node_id:
                    if not AnalysisNode.objects.filter(pk=e.node_id).exists():
                        data["deleted_nodes"] = [e.node_id]
                    else:
                        # We could potentially redirect to the same page with new version
                        # but client will probably just ditch this and re-request latest anyway
                        pass

        return JsonResponse(data, status=status)


class NodeJSONPostView(NodeJSONViewMixin):
    def post(self, request, *args, **kwargs):
        return self.json_response(request, *args, **kwargs)


class NodeJSONGetView(NodeJSONViewMixin):
    def get(self, request, *args, **kwargs):
        return self.json_response(request, *args, **kwargs)
