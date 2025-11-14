from urllib import request

from analysis.models import CandidateSearchRun, ReanalysisCandidate, Analysis
from snpdb.models import ProcessingStatus, Variant
from variantgrid import celery


@celery.shared_task
def reanalysis_new_annotation_task(candidate_search_run_id):
    candidate_search_run = CandidateSearchRun.get_for_user(request.user, pk=candidate_search_run_id)
    candidate_search_run.status = ProcessingStatus.PENDING
    candidate_search_run.save()

    records = []

    try:
        # Get out config
        analysis_qs = Analysis.objects.none()

        # Start searching
        for analysis in analysis_qs:
            variant = Variant.objects.first()
            clinvar = None
            notes = "TODO"
            evidence = {

            }

            records.append(ReanalysisCandidate(
                search_run=candidate_search_run,
                analysis=analysis,
                variant=variant,
                clinvar=clinvar,
                notes=notes,
                evidence=evidence,
            ))

        if records:
            ReanalysisCandidate.objects.bulk_create(records, batch_size=1000)
        candidate_search_run.status = ProcessingStatus.SUCCESS
    except Exception as e:
        candidate_search_run.status = ProcessingStatus.ERROR

    candidate_search_run.save()
