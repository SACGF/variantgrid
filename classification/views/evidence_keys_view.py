from rest_framework.decorators import permission_classes
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from classification.models.evidence_key import EvidenceKeyMap


@permission_classes((AllowAny, ))
class EvidenceKeysView(APIView):

    def get(self, request, **kwargs) -> Response:
        key_map = EvidenceKeyMap.instance()
        data = [k.to_json() for k in key_map.all_keys]
        return Response(status=HTTP_200_OK, data=data)
