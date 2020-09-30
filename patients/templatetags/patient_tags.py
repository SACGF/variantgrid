from django.template import Library

register = Library()

@register.inclusion_tag("patients/patient_modification_tag.html")
def show_patient_modification(patient_modification):
    return {'patient_modification': patient_modification}
