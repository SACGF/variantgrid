"""
Queues a relate of a patient's samples whenever which samples they have could have changed (#196): a sample
saved (eg linked to them) or two patients merged. A VCF's extract finishing relates its patients from
snpdb/tasks/somalier_tasks.py:somalier_vcf_id.
"""
from django.conf import settings
from django.db import transaction
from django.db.models.signals import post_save
from django.dispatch import receiver

from patients.models import Patient, patient_merged_signal
from snpdb.models.models_somalier import SomalierSampleExtract, get_sample_patient_ids
from snpdb.models.models_vcf import Sample
from snpdb.tasks.somalier_tasks import somalier_patient_relate


def _queue_patient_relate(patient_id: int):
    transaction.on_commit(somalier_patient_relate.si(patient_id).apply_async)


@receiver(post_save, sender=Sample)
def sample_somalier_patient_relate_handler(sender, instance: Sample, created: bool, **kwargs):
    """ The task skips a patient whose relate is already of these samples """
    if created or not settings.SOMALIER["enabled"]:
        return  # A new sample is related once its VCF is extracted
    if SomalierSampleExtract.objects.filter(sample=instance).exists():
        for patient_id in get_sample_patient_ids([instance]):
            _queue_patient_relate(patient_id)


@receiver(patient_merged_signal, sender=Patient)
def patient_merged_somalier_relate_handler(sender, patient: Patient, **kwargs):
    if settings.SOMALIER["enabled"]:
        _queue_patient_relate(patient.pk)
