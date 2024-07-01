from classification.enums import ShareLevel, EvidenceKeyValueType, SpecialEKeys
from classification.models import Classification, ClassificationJsonParams, ClassificationModification, EvidenceKeyMap, \
    ImportedAlleleInfo
from classification.models.classification_json_definitions import ClassificationJsonAlleleDict, \
    ClassificationJsonAlleleRevolvedDict
from genes.hgvs import CHGVS


def get_allele_info_dict(classification: Classification) -> ClassificationJsonAlleleDict:
    allele_info_dict: ClassificationJsonAlleleDict = {}
    if allele_info := classification.allele_info:
        resolved_dict: ClassificationJsonAlleleRevolvedDict = {
            "allele_id": allele_info.allele_id,
            "allele_info_id": allele_info.id,
            "allele_info_status": allele_info.status,
            "status": allele_info.status,
            "include": allele_info.latest_validation.include if allele_info.latest_validation else None,
            "variant_coordinate": allele_info.variant_coordinate,
        }

        if (genome_build := classification.get_genome_build_opt()) and \
                (preferred_build := allele_info[genome_build]) and \
                (c_hgvs := preferred_build.c_hgvs_obj):
            resolved_dict.update(c_hgvs.to_json())
        elif c_hgvs_raw := classification.get(SpecialEKeys.C_HGVS):
            resolved_dict.update(CHGVS(c_hgvs_raw).to_json())

        include = False
        if latest_validation := allele_info.latest_validation:
            include = latest_validation.include

        resolved_dict["include"] = include
        if warning_icon := ImportedAlleleInfo.icon_for(status=allele_info.status, include=include):
            resolved_dict.update(warning_icon.as_json())

        allele_info_dict["resolved"] = resolved_dict

        genome_builds = {}
        for variant_info in allele_info.resolved_builds:
            genome_builds[variant_info.genome_build.name] = {
                'variant_id': variant_info.variant_id,
                SpecialEKeys.C_HGVS: variant_info.c_hgvs
            }

        if genome_builds:
            allele_info_dict["genome_builds"] = genome_builds

    return allele_info_dict

#
# @dataclass(frozen=True)
# class ClassificationVersions:
#     version: ClassificationModification
#     latest_modification: ClassificationModification
#     latest_published: ClassificationModification
#
#     @staticmethod
#     def instance_from(classification: Classification, params: ClassificationJsonParams):
#         version = params.version
#         latest_modification: Optional[ClassificationModification] = None
#         last_published_version: Optional[ClassificationModification] = None
#         if version and not isinstance(version, ClassificationModification):
#             version = classification.modification_at_timestamp(version)
#         if version:
#             if version.is_last_published:
#                 last_published_version = version
#             if version.is_last_edited:
#                 latest_modification = version
#         if not latest_modification:
#             latest_modification = classification.last_edited_version
#         if not last_published_version:
#             last_published_version = classification.last_published_version
#
#         return ClassificationVersions(
#             version=version,
#             latest_modification=latest_modification,
#             latest_published=last_published_version
#         )
#
#     def json_for_version(self, version: ClassificationModification, current_user: User) -> ClassificationJsonVersionDict:
#         version_data: ClassificationJsonVersionDict = {
#             "version": version.created.timestamp(),
#             "publish_level": version.share_level_enum.value,
#             "is_published": version.published,
#             "can_write": version.is_last_edited and version.classification.can_write(user_or_group=current_user)
#         }
#         return version_data
#
# def populate_classification_json_v3(classification: Classification, params: ClassificationJsonParams) -> dict:
#
#     current_user = params.current_user
#     include_data = params.include_data
#     version = params.version
#     flatten = params.flatten
#     include_lab_config = params.include_lab_config
#     include_messages = params.include_messages
#     strip_complicated = params.strip_complicated
#     api_version = params.api_version
#     lowest_share_level = classification.lowest_share_level(current_user)
#
#     latest_version_mode = version is None
#
#     versions = ClassificationVersions.instance_from(classification, params)
#
#     content: ClassificationJsonDictv3 = {
#         "id": classification.id,
#         "lab_record_id": classification.lab_record_id,
#         "cr_lab_id": classification.cr_lab_id,
#     }
#     version_data = versions.json_for_version(versions.version, current_user=params.current_user)
#     published_data = version_data if versions.version == versions.latest_published else versions.json_for_version(versions.latest_published, current_user=params.current_user)
#     last_edited_data = published_data if versions.latest_modification == versions.latest_published else versions.json_for_version(versions.latest_modification, current_user=params.current_user)
#
#     content["version"] = version_data
#     content["version_published"] = published_data
#     content["version_latest"] = last_edited_data
#     content["flag_collection"] = classification.flag_collection_safe.pk if classification.id else None
#
#     use_evidence = classification.evidence if latest_version_mode else versions.version.evidence
#     content["data"] = classification.get_visible_evidence(use_evidence, lowest_share_level)
#
#     if include_messages:
#         content["messages"] = Classification.validate_evidence(use_evidence)
#
#     return content


def populate_classification_json(classification: Classification, params: ClassificationJsonParams) -> dict:
    # if params.api_version == 3:
    #     return populate_classification_json_v3(classification, params)

    current_user = params.current_user
    include_data = params.include_data
    version = params.version
    flatten = params.flatten
    include_lab_config = params.include_lab_config
    include_messages = params.include_messages
    strip_complicated = params.strip_complicated
    api_version = params.api_version

    if version and not isinstance(version, ClassificationModification):
        version = classification.modification_at_timestamp(version)

    include_data_bool = include_data or flatten

    title = classification.lab.group_name + '/' + classification.lab_record_id
    # have version and last_edited as separate as we might be looking at an older version
    # but still want to know last edited
    latest_modification = None
    last_published_version = None

    if version:
        if version.is_last_published:
            last_published_version = version
        if version.is_last_edited:
            latest_modification = version

    if not latest_modification:
        latest_modification = classification.last_edited_version
    if not last_published_version:
        last_published_version = classification.last_published_version

    can_write = classification.can_write(user_or_group=current_user)

    if classification.share_level == ShareLevel.CURRENT_USER.key:
        last_published_version = None

    last_published_timestamp = None
    lowest_share_level = classification.lowest_share_level(current_user)
    if last_published_version:
        last_published_timestamp = last_published_version.created.timestamp()

    clinical_context = None
    if classification.clinical_context:
        clinical_context = classification.clinical_context.name

    content = {
        'id': classification.id,
        'lab_record_id': classification.lab_record_id,
        'institution_name': classification.lab.organization.name,
        'lab_id': classification.lab.group_name,
        'cr_lab_id': classification.cr_lab_id,
        'org_name': classification.lab.organization.shortest_name,
        'lab_name': classification.lab.name,
        'title': title,
        'publish_level': classification.share_level,
        'version_publish_level': version.share_level if version else classification.share_level,
        'version_is_published': version.published if version else None,
        'is_last_published': version.is_last_published if version else (True if bool(classification.id) else None),
        'published_version': last_published_timestamp,
        'can_write': can_write,
        'can_write_latest': can_write,
        'clinical_context': clinical_context,
        'data': classification.get_visible_evidence(classification.evidence, lowest_share_level),
        'resolved_condition': classification.condition_resolution_dict,
        'withdrawn': classification.withdrawn,
        'allele_origin_bucket': classification.allele_origin_bucket,
        'condition_text_match': classification.condition_text_record.pk if classification.condition_text_record else None
    }
    if latest_modification:
        content['version'] = latest_modification.created.timestamp()

    if classification.id:
        # only get flag collection if we've saved a record
        content['flag_collection'] = classification.flag_collection_safe.pk

    if last_published_version and latest_modification:
        content['has_changes'] = last_published_version.id != latest_modification.id

    if include_messages:
        content['messages'] = Classification.validate_evidence(classification.evidence)

    if can_write and latest_modification:
        content['last_edited'] = latest_modification.created.timestamp()

    if version is None or version.is_last_edited:
        # show the editable version
        pass
    else:
        version_timestamp = version.created.timestamp()
        versioned_data = version.evidence
        content['title'] = content['title'] + '.' + str(version_timestamp)
        content['version'] = version_timestamp
        content['data'] = classification.get_visible_evidence(versioned_data, lowest_share_level)
        if include_messages:
            content['messages'] = Classification.validate_evidence(versioned_data)
        content['can_write'] = False

    if include_lab_config:
        content['config'] = classification.evidence_key_overrides

    if include_messages:
        content['messages'] = (content['messages'] or []) + classification.current_state_validation(
            user=current_user).flat_messages

    # Summary
    data = content['data']

    content["allele"] = get_allele_info_dict(classification)

    if classification.sample:
        content["sample_id"] = classification.sample.pk
        content["sample_name"] = classification.sample.name

    if include_data is not None and not isinstance(include_data, bool):
        if isinstance(include_data, str):
            include_data = [include_data]
        data_copy = {}
        for key in include_data:
            if key in data:
                data_copy[key] = data[key]
            else:
                data_copy[key] = None
        data = data_copy
    else:
        # do a copy just in case, want to make sure
        # that there are no errant saves going on
        data = content['data'].copy()

    if data and params.remove_acmg_namespace:
        data_copy = {}
        for key, value in data.items():
            if key.startswith("acmg:"):
                key = key[5:]
            data_copy[key] = value
        data = data_copy

    if data and params.fix_data_types:
        # currently just fixes data types to multiselect array
        # should fix more data types and move the logic out of as_json, so it can be done for other datatypes
        e_keys = EvidenceKeyMap.instance()
        for key, blob in data.items():
            if e_keys.get(key).value_type == EvidenceKeyValueType.MULTISELECT:
                value = blob.get('value')
                if isinstance(value, str):
                    blob['value'] = [s.strip() for s in value.split(',')]

    for key, blob in data.items():
        if strip_complicated:
            remove_keys = set(blob.keys()) - {'value', 'note', 'explain', 'db_refs'}
            for remove_key in remove_keys:
                blob.pop(remove_key)
        elif not include_messages:
            blob.pop('validation', None)

    if flatten:
        data = Classification.flatten(data, ignore_none_values=True)

    if include_data_bool:
        content['data'] = data
    else:
        content.pop('data', None)

    if api_version >= 2:
        META_KEYS = [
            'id', 'lab_record_id', 'institution_name', 'lab_id', 'lab_name', 'title',
            'user', 'published_version', 'can_write', 'can_write_latest',
            'clinical_context', 'flag_collection', 'has_changes', 'version',
            'last_edited'
        ]
        meta = {}
        for meta_key in META_KEYS:
            if meta_key in content:
                meta[meta_key] = content.pop(meta_key)
        if 'publish_level' in content:
            content['publish'] = content.pop('publish_level')

        content['meta'] = meta
        content['id'] = classification.lab.group_name + '/' + classification.lab_record_id

    if params.hardcode_extra_data:
        for key, value in params.hardcode_extra_data.items():
            content[key] = value

    return content
