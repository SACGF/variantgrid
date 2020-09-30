import django.dispatch

backend_vcf_import_start_signal = django.dispatch.Signal(providing_args=["backend_vcf"])
backend_vcf_import_success_signal = django.dispatch.Signal(providing_args=["backend_vcf"])
sequencing_run_created_signal = django.dispatch.Signal(providing_args=["sequencing_run"])
sequencing_run_sample_sheet_created_signal = django.dispatch.Signal(providing_args=["sample_sheet"])
sequencing_run_current_sample_sheet_changed_signal = django.dispatch.Signal(providing_args=["sequencing_run", "meaningfully_changed"])
