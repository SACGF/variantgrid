import collections
from typing import Optional, Dict, Any

from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import transaction
from django.utils.decorators import method_decorator
from django.utils.timezone import now
from lazy import lazy
from rest_framework.generics import get_object_or_404
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK, HTTP_400_BAD_REQUEST, \
    HTTP_500_INTERNAL_SERVER_ERROR
from rest_framework.views import APIView
from threadlocals.threadlocals import get_current_user

from classification.classification_stats import get_lab_gene_counts
from classification.enums import SubmissionSource, ShareLevel, ClinicalSignificance
from classification.models import ClassificationRef, ClassificationImport, \
    ClassificationJsonParams, PatchMeta, Classification
from classification.models.classification import ClassificationProcessError, \
    ClassificationModification
from classification.models.classification_patcher import patch_merge_age_units, patch_fuzzy_age
from classification.models.evidence_mixin import EvidenceMixin, VCStore
from classification.models.flag_types import classification_flag_types
from classification.tasks.classification_import_task import process_classification_import_task
from eventlog.models import create_event
from library.log_utils import report_exc_info, report_message
from library.utils import empty_to_none, DebugTimer
from snpdb.models import Lab, GenomeBuild
from snpdb.models.models_enums import ImportSource


class BulkInserter:

    def __init__(self, user: User, api_version=1, force_publish=False):
        """
        @param force_publish if True, then validation errors will not stop records
        from being published
        """
        self.user = user
        self._import_for_genome_build: Dict[Any, ClassificationImport] = dict()
        self.single_insert = False
        self.api_version = api_version
        self.force_publish = force_publish
        self.record_count = 0
        self.new_record_count = 0
        self.start = now()
        self.debug_timer = DebugTimer()

    def import_for(self, genome_build: GenomeBuild, transcript: str) -> Optional[ClassificationImport]:
        """
        Returns the ClassificationImport record that a classification should attach to to have its variant processed
        """
        if transcript:
            if not Classification.is_supported_transcript(transcript):
                report_message(message="Unsupported transcript type imported - will not attempt variant match", extra_data={
                    "target": transcript
                })
                return None

        if existing := self._import_for_genome_build.get(genome_build):
            return existing
        new_import = ClassificationImport.objects.create(genome_build=genome_build, user=self.user)
        self._import_for_genome_build[genome_build] = new_import
        return new_import

    def all_imports(self):
        return self._import_for_genome_build.values()

    @transaction.atomic
    def insert(self, data: dict, record_id: str = None, request=None):
        debug_timer = self.debug_timer

        if self.single_insert:
            raise ClassificationProcessError('If record_id is provided in URL cannot insert more than one record')

        self.record_count += 1
        data_copy = dict()
        data_copy.update(data)
        data = data_copy

        try:
            record_ref: ClassificationRef
            source = ClassificationView.verify_source(data)
            user = self.user

            if record_id:
                self.single_insert = True
                record_ref = ClassificationRef.init_from_str(user, record_id)
                if data.get('id'):
                    raise ClassificationProcessError('If providing id part, do not provide id in URL')
            else:
                id_part = data.pop('id', None)
                if not id_part:
                    keys_label = ', '.join(data.keys())
                    raise ClassificationProcessError(f'Must provide "id" segment with submission - keys = ({keys_label})')
                if isinstance(id_part, collections.Mapping):
                    id_part['user'] = user
                    record_ref = ClassificationRef.init_from_parts(**id_part)
                else:
                    record_ref = ClassificationRef.init_from_str(user, str(id_part))

            record_ref.check_security(must_be_writable=True)

            operation = None
            operation_data = None

            # if you say test=true, we just want validation messages back
            save = not data.pop('test', False)

            share_level = ShareLevel.from_key(data.pop('publish', None)) or ShareLevel.from_key(data.pop('share', None))
            if share_level == ShareLevel.CURRENT_USER:
                share_level = None

            delete_value = data.pop('delete', None)
            requested_delete = delete_value is True
            requested_undelete = delete_value is False

            debug_timer.tick("Security Check")

            for op in ['create', 'upsert', 'overwrite', 'data', 'patch']:
                op_data = data.pop(op, None)
                if op == 'data':
                    op = 'upsert'

                if op_data is not None:
                    if operation:
                        raise ClassificationProcessError('Can only provide 1 of create, overwrite, data, patch or upsert')
                    else:
                        operation = op
                        operation_data = op_data or {}

            # if we're deleting, don't do anything else
            if requested_delete:
                operation = None
                operation_data = None
                share_level = None

            if operation_data:
                operation_data = EvidenceMixin.to_patch(operation_data)

            # converts data before we even look at it
            # good place to manipulate multi-keys, e.g.
            # if we get age and age_units separately and want to merge
            def pre_process(base_data: VCStore):
                nonlocal operation_data
                patch_meta = PatchMeta(patch=operation_data, existing=base_data)
                patch_merge_age_units(patch_meta)
                if settings.VARIANT_CLASSIFICATION_AUTOFUZZ_AGE:
                    patch_fuzzy_age(patch_meta)

            record = None
            patch_messages = None
            patched_keys = None
            if operation:
                # operation modifiers
                immutable = source == SubmissionSource.API and not data.pop('editable', False)

                if record_ref.version is not None:
                    raise ClassificationProcessError("Can't perform operation on version")

                if operation == 'create' and record_ref.exists():
                    raise ClassificationProcessError('Record already exists, cannot create')

                if operation == 'patch' and not record_ref.exists():
                    raise ClassificationProcessError('Record does not exist, cannot patch')

                if operation in ('create', 'upsert', 'overwrite') and not record_ref.exists():
                    # creating a new record

                    debug_timer.tick("Preparing to Insert")

                    pre_process({})
                    record = record_ref.create(
                        source=source,
                        data=operation_data,
                        save=save,
                        make_fields_immutable=immutable)

                    debug_timer.tick("Inserted")

                    # We only want to link the variant on initial create - as patching may cause issues
                    # as classifications change variants. So the variant is immutable
                    if save:
                        self.new_record_count += 1

                        genome_build = record.get_genome_build()
                        classification_import = self.import_for(genome_build=genome_build, transcript=record.transcript)
                        record.classification_import = classification_import

                        # mark the fact that we're searching for a variant
                        if not record.variant:
                            if not classification_import:
                                record.set_variant(message=f"Transcript {record.transcript} is not of a currently supported type", failed=True)
                            else:
                                record.set_variant()

                        record.save()

                        debug_timer.tick("Prepare for Variant Resolution")

                        record.publish_latest(user=user)

                        debug_timer.tick("Published")
                else:
                    # patching existing records
                    record = record_ref.record
                    record.check_can_write(user)

                    # THE ACTUAL PATCHING OF DATA VALUES
                    pre_process(record.evidence)
                    patch_result = record.patch_value(operation_data,
                                                      clear_all_fields=operation == 'overwrite',
                                                      user=user,
                                                      source=source,
                                                      save=save,
                                                      make_patch_fields_immutable=immutable)
                    patch_messages = patch_result['messages']
                    patched_keys = patch_result['modified']

                    debug_timer.tick("Update Existing Record")
            else:
                record = record_ref.record

            if not record:
                if requested_delete:
                    return {"deleted": True}
                else:
                    raise ClassificationProcessError('Record not found')

            data_request = data.pop('return_data', False)
            flatten = data_request == 'flat'
            include_data = data_request or flatten
            if include_data == 'changes':
                include_data = bool(patched_keys)

            if patch_messages is None:
                patch_messages = []

            if data:
                patch_messages.append({'code': 'unexpected_parameters', 'message': 'Unexpected parameters ' + str(data)})

            if save:
                if requested_delete:
                    record.check_can_write(user)

                    params = ClassificationJsonParams(
                        version=record_ref.version,
                        include_data=include_data,
                        current_user=user,
                        flatten=flatten,
                        include_lab_config=data.pop('config', False),
                        api_version=self.api_version)

                    response_data = record.as_json(params)

                    if record.share_level_enum >= ShareLevel.ALL_USERS:
                        record.set_withdrawn(user=user, withdraw=True)
                        response_data['withdrawn'] = True
                    else:
                        response_data['deleted'] = True
                        record.delete()

                    return response_data

                elif share_level:
                    record.check_can_write(user)

                    if share_level < record.share_level_enum:
                        # raise ValueError('Record was already published at ' + record.share_level_enum.key + ', can no longer publish at ' + share_level.key)
                        patch_messages.append({'code': 'shared_higher', 'message': f"The record is already shared as '{record.share_level_enum.key}', re-sharing at that level instead of '{share_level}'"})
                        record.publish_latest(share_level=record.share_level_enum, user=user)

                    elif not self.force_publish and record.has_errors():
                        patch_messages.append({'code': 'share_failure', 'message': 'Cannot share record with errors'})

                    else:
                        record.publish_latest(share_level=share_level, user=user)
                        patch_messages.append({'code': 'shared', 'message': 'Latest revision is now shared with ' + str(share_level)})

                if requested_undelete:
                    record.set_withdrawn(user=user, withdraw=False)

                if record.withdrawn:
                    # Withdrawn records are no longer automatically un-withdrawn
                    # record.set_withdrawn(user=user, withdraw=False)
                    patch_messages.append({'code': 'withdrawn', 'message': 'The record you updated is withdrawn and will not appear to users'})

            else:
                patch_messages.append({'code': 'test_mode', 'message': 'Test mode on, no changes have been saved'})

            # Check to see if the current record matches the previous one exactly
            # then if the most recent record isn't marked as published, delete it.
            # - note we only compare against the very immediate previous record
            recent_modifications = list(ClassificationModification.objects.filter(classification=record.id).order_by('-created')[:2])
            resolution = 'closed'
            if recent_modifications:
                latest_mod = recent_modifications[0]
                if not latest_mod.is_last_edited:
                    raise ValueError('Assertion error - expected most recent record to be last edited')

                if len(recent_modifications) >= 2:
                    previous_mod = recent_modifications[1]

                    # if this latest record isn't published, and it matches the previous record exactly, ew can delete it
                    if not latest_mod.published and latest_mod.evidence == previous_mod.evidence:
                        latest_mod.delete()
                        previous_mod.is_last_edited = True
                        previous_mod.save()
                        latest_mod = previous_mod

                has_modifications = not latest_mod.published

                has_errors = record.has_errors()
                resolution = 'closed'
                if has_errors:
                    resolution = 'in_error'
                elif has_modifications:
                    resolution = 'open'

            if save:
                record.flag_collection_safe.ensure_resolution(
                    flag_type=classification_flag_types.classification_outstanding_edits,
                    resolution=resolution
                )

            json_data = record.as_json(ClassificationJsonParams(
                version=record_ref.version,
                include_data=include_data or not save,
                current_user=user,
                flatten=flatten,
                include_lab_config=data.pop('config', False),
                api_version=self.api_version
            ))

            debug_timer.tick("Generate Response")

            if patch_messages:
                json_data['patch_messages'] = patch_messages
            return json_data

        except ClassificationProcessError as ve:
            # expected error
            report_exc_info(request=request)
            return {'fatal_error': str(ve)}
        except Exception as e:
            # unexpected error
            report_exc_info(request=request)
            return {'internal_error': str(e)}

    def finish(self):
        debug_timer: DebugTimer = self.debug_timer

        if settings.VARIANT_CLASSIFICATION_MATCH_VARIANTS:
            for vc_import in self.all_imports():
                task = process_classification_import_task.si(vc_import.pk, ImportSource.API)
                task.apply_async()

        debug_timer.tick("Setup Async Variant Matching")

        if count := self.record_count:
            time_taken = now() - self.start
            total_time = time_taken.total_seconds()
            time_per_record = total_time / count
            create_event(user=get_current_user(), name="classification_import", details=f"{count} records imported {self.new_record_count} new, avg record processing {time_per_record:.3f}s\n{debug_timer}")


class ClassificationView(APIView):
    api_version = 1

    @staticmethod
    def verify_source(data) -> SubmissionSource:
        source = SubmissionSource(data.pop('source', SubmissionSource.API))
        if not source.is_valid_user_source():
            raise ValueError('Illegal value for source, should be "' + SubmissionSource.API.value + '" or "' + SubmissionSource.FORM.value + '"')
        return source

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

        user = request.user
        record_id = empty_to_none(kwargs.get('record_id', None))

        importer = BulkInserter(user=user, api_version=self.api_version)
        data = request.data
        if isinstance(data, list):
            json_data = []
            for record_data in data:
                result = importer.insert(record_data)
                json_data.append(result)

        else:
            if isinstance(data.get('records'), list):
                json_data = []
                for record_data in data.get('records'):
                    result = importer.insert(record_data)
                    json_data.append(result)
                json_data = {"results": json_data}
            else:
                json_data = importer.insert(data, record_id)
                if 'fatal_error' in json_data:
                    return Response(status=HTTP_400_BAD_REQUEST, data=json_data)
                elif 'internal_error' in json_data:
                    return Response(status=HTTP_500_INTERNAL_SERVER_ERROR, data=json_data)

        importer.finish()

        return Response(status=HTTP_200_OK, data=json_data)


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
