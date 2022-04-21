from typing import Optional

from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.utils.decorators import method_decorator
from rest_framework.generics import get_object_or_404
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK, HTTP_400_BAD_REQUEST, \
    HTTP_500_INTERNAL_SERVER_ERROR
from rest_framework.views import APIView
from classification.classification_stats import get_lab_gene_counts
from classification.enums import ClinicalSignificance
from classification.models import ClassificationRef, ClassificationJsonParams
from classification.models.classification_import_run import ClassificationImportRun
from classification.models.classification_inserter import BulkClassificationInserter
from library.utils import empty_to_none
from snpdb.models import Lab
from uicore.json.to_json import force_json


class ClassificationView(APIView):
    api_version = 1

    @method_decorator(login_required)
    def get(self, request, **kwargs) -> Response:
        """
        Provide an id in the URL to retrieve that specific record, otherwise get
        header information on all available records
        """
        record_id = empty_to_none(kwargs.get('record_id', None))
        if record_id is None:
            # bulk get has been deprecated on this URL for now
            return Response(status=HTTP_200_OK, data={})

        else:
            record_ref = ClassificationRef.init_from_str(request.user, record_id)
            record_ref.check_exists()
            record_ref.check_security()

            include_config = request.query_params.get('config', '').lower() == 'true'

            flatten = request.query_params.get('data', '').lower() == 'flat'
            include_data = request.query_params.get('data', '').lower() != 'false'

            response_data = record_ref.as_json(ClassificationJsonParams(
                current_user=request.user,
                include_data=include_data,
                flatten=flatten,
                include_lab_config=include_config,
                api_version=self.api_version))

            return Response(status=HTTP_200_OK, data=response_data)

    @method_decorator(login_required)
    def post(self, request, **kwargs) -> Response:
        """ Create a new record """

        if not settings.UPLOAD_ENABLED:
            raise PermissionDenied("Uploads are currently disabled (settings.UPLOAD_ENABLED=False)")

        user: User = request.user
        record_id = empty_to_none(kwargs.get('record_id', None))

        importer = BulkClassificationInserter(user=user, api_version=self.api_version)
        data = request.data
        if isinstance(data, list):
            json_data = []
            for record_data in data:
                result = importer.insert(record_data)
                json_data.append(result)

        else:
            if isinstance(data.get('records'), list):

                records = data.get('records')
                complete_identifier = None
                classification_import_run: Optional[ClassificationImportRun] = None

                if import_id := data.get('import_id'):
                    # prefix import_id with username, so users can't overwrite each other
                    import_id = f"{user.username}#{import_id}"
                    status = data.get('status')
                    completed = status == 'complete'
                    if completed:
                        complete_identifier = import_id
                    classification_import_run = ClassificationImportRun.record_classification_import(identifier=import_id)

                per_json_data = list()
                for record_data in records:
                    result = importer.insert(record_data, import_run=classification_import_run)
                    if classification_import_run:
                        classification_import_run.increment_status(result.status)
                    per_json_data.append(result)

                if classification_import_run:
                    classification_import_run.save()
                json_data = {"results": per_json_data}

                if complete_identifier:
                    # have to mark is_complete=True at the end after import has run
                    ClassificationImportRun.record_classification_import(
                        identifier=import_id,
                        is_complete=True)

            else:
                # single record
                json_data = importer.insert(data, record_id)
                if json_data.internal_error:
                    return Response(status=HTTP_500_INTERNAL_SERVER_ERROR, data=force_json(json_data))

        importer.finish()

        return Response(status=HTTP_200_OK, data=force_json(json_data))


class LabGeneClassificationCountsView(APIView):
    """ Returns a dict of {gene_symbol: {clinical_significance: classification_counts}} """

    def get(self, request, *args, **kwargs):
        lab = get_object_or_404(Lab, pk=kwargs["lab_id"])

        lab_gene_counts = get_lab_gene_counts(request.user, lab)
        classification_counts = {}
        for gene_symbol, clinical_significance_count in lab_gene_counts.items():
            total = 0
            summary = []
            for cs, cs_display in ClinicalSignificance.SHORT_LABELS.items():
                if count := clinical_significance_count.get(cs):
                    summary.append(f"{cs_display}: {count}")
                    total += count
            classification_counts[gene_symbol] = {
                "total": total,
                "summary": ", ".join(summary),
            }

        data = {
            "lab_name": lab.name,
            "gene_symbols": lab_gene_counts.keys(),
            "classification_counts": classification_counts,
        }
        return Response(data)
