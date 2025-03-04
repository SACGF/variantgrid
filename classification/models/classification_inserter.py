from functools import cached_property
from typing import Optional, Mapping

from django.contrib.auth.models import User
from django.db import transaction
from django.utils.timezone import now

# from classification.enums import ForceUpdate
from classification.enums import ShareLevel, SubmissionSource, SpecialEKeys
from classification.models import ClassificationProcessError, ClassificationRef, \
    EvidenceMixin, classification_flag_types, ClassificationJsonParams, ClassificationModification, \
    ClassificationPatchResponse, ClassificationImportRun
from classification.models.classification_utils import ClassificationPatchStatus
from classification.models.variant_resolver import VariantResolver
from eventlog.models import create_event
from library.log_utils import report_exc_info
from library.utils import DebugTimer, reset_timer


class BulkClassificationInserter:

    def __init__(self, user: User, api_version=1, force_publish=False):
        """
        @param force_publish if True, then validation errors will not stop records
        from being published
        """
        self.user = user
        self.variant_resolver = VariantResolver(user=user)
        self.single_insert = False
        self.api_version = api_version
        self.force_publish = force_publish
        self.record_count = 0
        self.new_record_count = 0
        self.start = now()

    @cached_property
    def debug_timer(self) -> DebugTimer:
        return reset_timer()

    @staticmethod
    def verify_source(data) -> SubmissionSource:
        source = SubmissionSource(data.pop('source', SubmissionSource.API))
        if not source.is_valid_user_source():
            raise ValueError(
                'Illegal value for source, should be "' + SubmissionSource.API.value + '" or "' + SubmissionSource.FORM.value + '"')
        return source

    @transaction.atomic
    def insert(
            self,
            data: dict,
            record_id: Optional[str] = None,
            submission_source: Optional[SubmissionSource] = None,
            import_run: Optional[ClassificationImportRun] = None) -> ClassificationPatchResponse:

        debug_timer = self.debug_timer
        patch_response = ClassificationPatchResponse()

        if self.single_insert:
            raise ClassificationProcessError('If record_id is provided in URL cannot insert more than one record')

        self.record_count += 1
        data_copy = {}
        data_copy.update(data)
        data = data_copy
        record_ref: Optional[ClassificationRef] = None

        try:
            record_ref: ClassificationRef
            source = submission_source or BulkClassificationInserter.verify_source(data)
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
                    raise ClassificationProcessError(
                        f'Must provide "id" segment with submission - keys = ({keys_label})')
                if isinstance(id_part, Mapping):
                    id_part['user'] = user
                    record_ref = ClassificationRef.init_from_parts(**id_part)
                else:
                    record_ref = ClassificationRef.init_from_str(user, str(id_part))

            record_ref.check_security(must_be_writable=True)

            operation: Optional[str] = None
            operation_data: Optional[dict] = None

            # if you say test=true, we just want validation messages back
            save = not data.pop('test', False)

            share_level = ShareLevel.from_key(data.pop('publish', None)) or ShareLevel.from_key(data.pop('share', None))
            if share_level == ShareLevel.CURRENT_USER:
                share_level = None

            delete_value = data.pop('delete', None)
            delete_reason = data.pop('delete_reason', None)
            # delete can be
            # True : delete (or withdraw if shared)
            # False : unwithdraw
            # "withdraw": (withdraw regardless of share level)
            requested_withdraw = delete_value == "withdraw"
            requested_delete = delete_value is True
            requested_undelete = delete_value is False

            debug_timer.tick("Security Check")

            # -- ENSURE REQUEST IS WELL FORMATTED --

            # note we don't actually use ForceUpdate for anything at the moment
            # at one time (never in prod) was used to stop changes to only source ID, curation date
            # force: Optional[ForceUpdate] = None
            # if force_str := data.pop('force', None):
            #     try:
            #         force = ForceUpdate(force_str)
            #     except ValueError:
            #         pass

            for op in ['create', 'upsert', 'overwrite', 'data', 'patch', 'patch-empty']:
                op_data = data.pop(op, None)
                if op == 'data':
                    op = 'upsert'

                if op_data is not None:
                    if operation:
                        raise ClassificationProcessError(
                            'Can only provide 1 of create, overwrite, data, patch or upsert')
                    else:
                        operation = op
                        operation_data = op_data or {}

            # if we're deleting, don't do anything else
            if requested_delete or requested_withdraw:
                operation = None
                operation_data = None
                share_level = None

            new_source_id: Optional[str] = None
            if operation_data:
                if source_id := operation_data.pop(SpecialEKeys.SOURCE_ID, None):
                    if isinstance(source_id, Mapping):
                        source_id = source_id.get('value')
                    new_source_id = source_id

                operation_data = EvidenceMixin.to_patch(operation_data)

            record = None
            is_new_record = False
            if operation:
                # operation modifiers
                immutable = source == SubmissionSource.API and not data.pop('editable', False)

                if record_ref.version is not None:
                    raise ClassificationProcessError("Can't perform operation on version")

                if operation == 'create' and record_ref.exists():
                    raise ClassificationProcessError('Record already exists, cannot create')

                if operation in ('patch', 'patch-empty') and not record_ref.exists():
                    raise ClassificationProcessError('Record does not exist, cannot patch')

                if operation in ('create', 'upsert', 'overwrite') and not record_ref.exists():
                    # creating a new record

                    debug_timer.tick("Preparing to Insert")

                    patch_response.status = ClassificationPatchStatus.NEW
                    record, create_response = record_ref.create_with_response(
                        source=source,
                        data=operation_data,
                        save=save,
                        make_fields_immutable=immutable)

                    patch_response += create_response
                    debug_timer.tick("Inserted")

                    # We only want to link the variant on initial create - as patching may cause issues
                    # as classifications change variants. So the variant is immutable
                    if save:
                        self.new_record_count += 1

                        # FIXME here's where we can see if we've already matched a varaint

                        if allele_info := record.attempt_set_variant_info_from_pre_existing_imported_allele_info():
                            # note will only queue to resolve if it's a newly created allele info
                            self.variant_resolver.queue_resolve(allele_info)

                        record.save()

                        # TODO remove this code if publishing just the once works
                        # debug_timer.tick("Prepare for Variant Resolution")
                        # PUBLISH NUMBER 1
                        # record.publish_latest(user=user)
                        is_new_record = True

                        # debug_timer.tick("Publish Complete")
                else:
                    ignore_if_only_patch: Optional[set[str]] = None
                    # if some updates are not-meaningful and by themselves shouldn't cause a new version of the record to be made
                    # if source == SubmissionSource.API and not force:
                    #     ignore_if_only_patch = {"curation_date", "source_id"}

                    # patching existing records
                    record = record_ref.record
                    record.check_can_write(user)

                    # THE ACTUAL PATCHING OF DATA VALUES
                    patch_response += record.patch_value(operation_data,
                                                         clear_all_fields=operation == 'overwrite',
                                                         user=user,
                                                         source=source,
                                                         save=save,
                                                         make_patch_fields_immutable=immutable,
                                                         leave_existing_values=operation.endswith('-empty'),
                                                         ignore_if_only_patching=ignore_if_only_patch)

                    patch_response.status = ClassificationPatchStatus.UPDATE if patch_response.modified_keys else ClassificationPatchStatus.NO_CHANGE

                    debug_timer.tick("Update Existing Record")
            else:
                record = record_ref.record

            if not record:
                if requested_delete:
                    patch_response.deleted = True
                    return patch_response
                else:
                    raise ClassificationProcessError('Record not found')

            if import_run:
                record.last_import_run = import_run
            if new_source_id:
                record.last_source_id = new_source_id

            data_request = data.pop('return_data', False)
            flatten = data_request == 'flat'
            include_data = data_request or flatten
            if include_data == 'changes':
                include_data = bool(patch_response.modified_keys)

            include_lab_config = data.pop('config', False)

            if data:
                patch_response.append_warning(code="unexpected_parameters", message=f"Unexpected parameters {data}")

            has_published = False
            if save:
                if requested_delete or requested_withdraw:
                    # deleting or withdrawing
                    record.check_can_write(user)

                    params = ClassificationJsonParams(
                        version=record_ref.version,
                        include_data=include_data,
                        current_user=user,
                        flatten=flatten,
                        include_lab_config=include_lab_config,
                        api_version=self.api_version)

                    patch_response.classification_json = record.as_json(params)

                    if record.share_level_enum >= ShareLevel.ALL_USERS or requested_withdraw:
                        # withdraw if shared or if only asking for withdraw
                        patch_response.status = ClassificationPatchStatus.WITHDRAWN if not record.withdrawn else ClassificationPatchStatus.NO_CHANGE
                        record.set_withdrawn(user=user, withdraw=True, reason=delete_reason)
                        patch_response.withdrawn = True
                    else:
                        patch_response.status = ClassificationPatchStatus.DELETED
                        patch_response.deleted = True
                        record.delete()

                    return patch_response

                elif share_level:
                    record.check_can_write(user)

                    if share_level < record.share_level_enum:
                        # raise ValueError('Record was already published at ' + record.share_level_enum.key + ', can no longer publish at ' + share_level.key)
                        patch_response.append_warning(code="shared_higher",
                                                      message=f"The record is already shared as '{record.share_level_enum.key}', re-sharing at that level instead of '{share_level}'")
                        patch_response.published = record.publish_latest(share_level=record.share_level_enum, user=user)
                        has_published = True
                        if patch_response.published:
                            patch_response.saved = True

                    elif not self.force_publish and record.has_errors():
                        patch_response.append_warning(code="share_failure", message="Cannot share record with errors")

                    else:
                        # 2ND PUBLISH
                        record.publish_latest(share_level=share_level, user=user)
                        has_published = True
                        patch_response.append_warning(code="shared",
                                                      message=f"Latest revision is now shared with {share_level}")

                if requested_undelete:
                    if record.withdrawn:
                        patch_response.status = ClassificationPatchStatus.UN_WITHDRAWN
                        record.set_withdrawn(user=user, withdraw=False)

                elif record.withdrawn:
                    # Withdrawn records are no longer automatically un-withdrawn
                    # record.set_withdrawn(user=user, withdraw=False)
                    patch_response.status = ClassificationPatchStatus.ALREADY_WITHDRAWN
                    patch_response.append_warning(code="withdrawn",
                                                  message="The record you updated is withdrawn and will not appear to users")

            else:
                patch_response.append_warning(code="test_mode", message="Test mode on, no changes have been saved")

            if save and is_new_record and not has_published:
                record.publish_latest(user=user)

            # Check to see if the current record matches the previous one exactly
            # then if the most recent record isn't marked as published, delete it.
            # - note we only compare against the very immediate previous record
            recent_modifications = list(ClassificationModification.objects.filter(classification=record.id).order_by('-created')[:2])
            resolution = 'closed'
            if recent_modifications:
                latest_mod = recent_modifications[0]
                if not latest_mod.is_last_edited:
                    raise ValueError('Assertion error - expected most recent record to be last edited')

                # detect if a value was changed and then changed back (if so, delete the pointless change)
                if len(recent_modifications) >= 2:
                    previous_mod = recent_modifications[1]

                    # if this latest record isn't published, and it matches the previous record exactly, we can delete it
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
                if not patch_response.saved and (import_run or new_source_id):
                    # we would be here if no data was patched
                    # update the classification with the new import_run, but don't update the modified date
                    # since nothing has otherwise changed
                    record.update_modified = False
                    record.save()

            json_data = record.as_json(ClassificationJsonParams(
                version=record_ref.version,
                include_data=include_data or not save,
                current_user=user,
                flatten=flatten,
                include_lab_config=include_lab_config,
                api_version=self.api_version
            ))

            debug_timer.tick("Generate Response")

            patch_response.classification_json = json_data
            return patch_response

        except ClassificationProcessError as ve:
            # expected error
            report_exc_info({"target": str(record_ref) if record_ref else "unknown"})
            patch_response.internal_error = ve
            return patch_response
        except Exception as e:
            # unexpected error
            report_exc_info({"target": str(record_ref) if record_ref else "unknown"})
            patch_response.internal_error = e
            return patch_response

    def finish(self):
        debug_timer: DebugTimer = self.debug_timer

        self.variant_resolver.process_queue()

        debug_timer.tick("Setup Async Variant Matching")

        if self.record_count > 1:
            time_taken = now() - self.start
            total_time = time_taken.total_seconds()
            time_per_record = total_time / self.record_count
            create_event(user=self.user, name="classification_import",
                         details=f"{self.record_count} records imported {self.new_record_count} new, avg record processing {time_per_record:.3f}s\n{debug_timer}")
