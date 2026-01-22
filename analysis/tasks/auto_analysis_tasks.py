import celery

from analysis.analysis_templates import auto_launch_analysis_templates_for_sample
from snpdb.models import VCF, Sample


@celery.shared_task
def auto_run_analyses_for_sample(sample_id: int, analysis_description: str, skip_already_analysed: bool):
    sample = Sample.objects.get(pk=sample_id)
    user = sample.vcf.user
    auto_launch_analysis_templates_for_sample(user, sample,
                                              analysis_description=analysis_description,
                                              skip_already_analysed=skip_already_analysed)


@celery.shared_task
def auto_run_analyses_for_vcf(vcf_id: int, analysis_description: str, skip_already_analysed: bool):
    vcf = VCF.objects.get(pk=vcf_id)

    for sample in vcf.sample_set.all():
        auto_launch_analysis_templates_for_sample(vcf.user, sample,
                                                  analysis_description=analysis_description,
                                                  skip_already_analysed=skip_already_analysed)
