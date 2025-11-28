import itertools

from celery import Task

from analysis.models import CandidateSearchRun, Candidate
from library.log_utils import get_traceback, log_traceback
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
            log_traceback()
            candidate_search_run.error_exception = get_traceback()
            candidate_search_run.status = ProcessingStatus.ERROR

        candidate_search_run.save()

    @staticmethod
    def limit_sample_and_results(sample_records, max_samples, max_results):
        num_results = 0
        for sample, results in itertools.islice(sample_records, max_samples):
            max_remaining = max_results - num_results
            if max_remaining <= 0:
                break
            limited = results[:max_remaining]
            yield sample, limited
            num_results += len(limited)
            if num_results >= max_results:
                break

