from django.contrib.auth.models import User
from django.dispatch.dispatcher import receiver

from classification.enums.classification_enums import SpecialEKeys, ValidationCode
from classification.models import PatchMeta
from classification.models.classification import Classification, classification_validation_signal
from classification.models.classification_utils import ValidationMerger
from library.guardian_utils import user_in_group


@receiver(classification_validation_signal, sender=Classification)
def assign_classification_user(sender, **kwargs) -> ValidationMerger:  # pylint: disable=unused-argument
    """
    If the owner evidence key has been set, see if that corresponds to a user in the system, and re-assign the user value.
    Note that the user field doesn't get used much anymore
    """

    """ Checks the owner is valid """
    patch_meta: PatchMeta = kwargs.get('patch_meta')
    record: Classification = kwargs.get('classification')

    vm = ValidationMerger()

    if patch_meta.is_modified(SpecialEKeys.OWNER):
        owner_text = patch_meta.get(SpecialEKeys.OWNER)
        user = User.objects.filter(username=owner_text).first()
        if user:
            # this check might be done before the record permissions have been set
            # so check for lab or ability to write to the record
            if not (user_in_group(user, record.lab.group) or record.can_write(user_or_group=user)):
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
                if record.id:
                    record.fix_permissions()

        # No message when the owner isn't a user here, as it would expose owner names on Shariant
    return vm
