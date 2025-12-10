import django

vcf_import_success_signal = django.dispatch.Signal(providing_args=["vcf"])