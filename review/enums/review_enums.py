"""
Enums for reviews: QuestionValueType (the kind of answer a ReviewQuestion takes), and the ReviewMedium /
ReviewParticipants choices a Review records for how it was held and who took part. Stored codes: add values,
never change one.
"""
from django.db.models import TextChoices


class QuestionValueType(TextChoices):
    """
    The question type in a review topic, currently only supports disagreement, but we could support more in future
    """
    Disagreement = 'D', "Disagreement"


# TODO are we ReviewMedium or ReviewMethod or something else?
# TODO should these options be made configurable, or based on the topic?
class ReviewMedium(TextChoices):
    email = "email", "Email"
    phone = "phone", "Phone"
    video = "mtm", "Multidisciplinary Team Meeting"


class ReviewParticipants(TextChoices):
    curation = "curation", "Curation Scientists"
    clinicians = "clinicians", "Clinicians"
    external_experts = "external_experts", "External Experts"
