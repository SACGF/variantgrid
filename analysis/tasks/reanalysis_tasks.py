from analysis.models import Analysis, Candidate
from analysis.tasks.abstract_candidate_search_task import AbstractCandidateSearchTask
from snpdb.models import Variant
from variantgrid.celery import app


class ReAnalysisNewAnnotationTask(AbstractCandidateSearchTask):
    def get_candidate_records(self, candidate_search_run):
        # Get out config
        analysis_qs = Analysis.objects.none()

        # Start searching
        records = []
        for analysis in analysis_qs:
            variant = Variant.objects.first()
            clinvar = None
            notes = "TODO"
            evidence = {

            }

            records.append(Candidate(
                search_run=candidate_search_run,
                analysis=analysis,
                variant=variant,
                clinvar=clinvar,
                notes=notes,
                evidence=evidence,
            ))
        return records

ReAnalysisNewAnnotationTask = app.register_task(ReAnalysisNewAnnotationTask())
