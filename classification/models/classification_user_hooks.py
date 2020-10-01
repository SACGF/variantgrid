from django.contrib.auth.models import User
from django.dispatch.dispatcher import receiver

from classification.enums.classification_enums import SpecialEKeys, \
    ValidationCode
from classification.models import PatchMeta
from classification.models.classification import Classification, \
    classification_validation_signal
from classification.models.classification_utils import ValidationMerger


@receiver(classification_validation_signal, sender=Classification)
def validate_variant_fields(sender, **kwargs) -> ValidationMerger:
    """ Checks the owner is valid """
    patch_meta: PatchMeta = kwargs.get('patch_meta')
    record: Classification = kwargs.get('classification')

    vm = ValidationMerger()

    if patch_meta.is_modified(SpecialEKeys.OWNER):
        # user_correct = False
        owner_text = patch_meta.get(SpecialEKeys.OWNER)
        user = User.objects.filter(username=owner_text).first()
        if user:
            # this check might be done before the record permissions have been set
            # so check for lab or ability to write to the record
            if not (user.groups.filter(id=record.lab.group.id).exists() or record.can_write(user=user)):
                message = f"User does not belong to lab {record.lab.name} so can't be owner of this record"
                vm.add_message(
                    key=SpecialEKeys.OWNER,
                    code=ValidationCode.NOT_LAB,
                    severity='warning',
                    message=message,
                )
            else:
                # user is valid for this lab (currently anyway)
                # assign the record to the user.
                record.user = user
                # user_correct = True
                if record.id:
                    record.fix_permissions()

        # removed this message as it was exposing owner names in Shariant where we didn't want to see them
        # if not user_correct:
        #     message = f'"{owner_text}" is not currently a user in this system, record is assigned to "{record.user.username}"'
        #     vm.add_message(
        #         key=SpecialEKeys.OWNER,
        #         code=ValidationCode.USER_NOT_FOUND,
        #         severity='info',
        #         message=message,
        #     )
    return vm
