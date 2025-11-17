from celery import Task

from analysis.models import CandidateSearchRun, Candidate
from snpdb.models import ProcessingStatus


class AbstractCandidateSearchTask(Task):

    def get_candidate_records(self, candidate_search_run):
        raise NotImplementedError()

    def run(self, candidate_search_run_id):
        csr_qs = CandidateSearchRun.objects.filter(pk=candidate_search_run_id)
        csr_qs.update(status=ProcessingStatus.PROCESSING)
        candidate_search_run = csr_qs.get()

        try:
            records = self.get_candidate_records(candidate_search_run)
            if records:
                Candidate.objects.bulk_create(records, batch_size=1000)
            candidate_search_run.status = ProcessingStatus.SUCCESS
        except Exception as e:
            candidate_search_run.error_exception = str(e)
            candidate_search_run.status = ProcessingStatus.ERROR

        candidate_search_run.save()

