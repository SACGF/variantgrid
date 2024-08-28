from django.contrib.auth.models import User
from django.dispatch import receiver

from classification.enums import ShareLevel
from classification.models import classification_post_publish_signal, Classification, ClassificationModification
from classification.models.classification_grouping import ClassificationGrouping
from library.utils import DebugTimer


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              debug_timer: DebugTimer,
              **kwargs):  # pylint: disable=unused-argument
    ClassificationGrouping.assign_grouping_for_classification(classification)