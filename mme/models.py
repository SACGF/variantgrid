from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone


class MMESubmissionStatus(models.TextChoices):
    DRAFT = "D", "Draft"
    SUBMITTED = "S", "Submitted"
    ERROR = "E", "Error"


class MMESubmission(models.Model):
    """ A single outbound submission of one Classification to one remote MME node.
        The MME 'patient' profile is synthesised from the classification: candidate
        variant from its allele, disorders from its condition-under-curation, and HPO
        features from a linked patient (if any). The unit is the classification, not
        the Patient, because classifications are the consented (share_level=PUBLIC),
        structured unit and are frequently not linked to a Patient at all. """
    classification = models.ForeignKey("classification.Classification", on_delete=models.CASCADE)
    node_id = models.CharField(max_length=64)          # key into settings.MME_NODES
    # Opaque, stable, non-PII id we send as patient.id (e.g. the classification's
    # lab record id / UUID).
    external_patient_id = models.CharField(max_length=255)
    status = models.CharField(max_length=1, choices=MMESubmissionStatus.choices,
                              default=MMESubmissionStatus.DRAFT)
    request_json = models.JSONField(null=True)          # exact profile we POSTed
    response_json = models.JSONField(null=True)         # raw /match response
    error = models.TextField(null=True, blank=True)
    created = models.DateTimeField(default=timezone.now)
    submitted = models.DateTimeField(null=True, blank=True)
    submitted_by = models.ForeignKey(User, null=True, on_delete=models.SET_NULL)

    class Meta:
        unique_together = ("classification", "node_id")

    def __str__(self):
        return f"MMESubmission(classification={self.classification_id}, node={self.node_id}, {self.status})"


class MMEMatchResult(models.Model):
    """ One candidate patient returned in a /match response (inbound or outbound). """
    submission = models.ForeignKey(MMESubmission, null=True, on_delete=models.CASCADE)
    score = models.FloatField()
    matched_patient_id = models.CharField(max_length=255)   # remote node's patient id
    contact_name = models.TextField(null=True)
    contact_href = models.TextField(null=True)
    patient_json = models.JSONField()                       # full returned patient obj
    created = models.DateTimeField(default=timezone.now)

    def __str__(self):
        return f"MMEMatchResult(patient={self.matched_patient_id}, score={self.score})"


class MMEInboundQuery(models.Model):
    """ Audit row for an inbound /match query served by us: who queried, when, what
        profile, and how many of our patients we returned. Keeps a compliance trail. """
    request_json = models.JSONField()
    num_results = models.IntegerField(default=0)
    created = models.DateTimeField(default=timezone.now)

    def __str__(self):
        return f"MMEInboundQuery({self.created}, num_results={self.num_results})"
