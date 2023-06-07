from datetime import timedelta

from django.core.management import BaseCommand

from annotation.models import Citation, CitationFetchRequest
from annotation.models.models_citations import CitationSource
from library.utils import batch_iterator

error_citations = list(Citation.objects.filter(error__isnull=False))
for batch in batch_iterator(error_citations):
    CitationFetchRequest.fetch_all_now(batch, cache_age=timedelta(seconds=0))


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--load', action='store_true', default=False, help='Pre-emptively load all citations not yet loaded')
        parser.add_argument('--fix', action='store_true', default=False, help='Force re-load citations that appear to have an incorrect id')
        parser.add_argument('--reload', action='store_true', default=False, help='Attempts to reload all citations currently marked as in error')

    def handle(self, *args, **options):
        if options["reload"]:
            unloaded_citations = Citation.objects.filter(error__isnull=False)
            print(f"{unloaded_citations.count()} error citations")
            loaded = 0
            for batch in batch_iterator(unloaded_citations, batch_size=20):
                CitationFetchRequest.fetch_all_now(batch)
                loaded += len(batch)
                print(f"Loaded {loaded} citations")
            print(f"{unloaded_citations.count()} error citations remaining - hopefully these are just invalid IDs")

        if options["load"]:
            unloaded_citations = Citation.objects.filter(last_loaded__isnull=True)
            print(f"{unloaded_citations.count()} unloaded citations")
            loaded = 0
            for batch in batch_iterator(unloaded_citations, batch_size=20):
                CitationFetchRequest.fetch_all_now(batch)
                loaded += len(batch)
                print(f"Loaded {loaded} citations")

        if options["fix"]:
            bad_citations = list()
            for citation in Citation.objects.filter(old_id__isnull=False, data_json__isnull=False):
                if citation.source == CitationSource.NCBI_BOOKSHELF:
                    nbk_id = citation.data_json.get("RID")
                    if nbk_id != citation.index:
                        bad_citations.append(citation)
                        print(f"{citation} has data from wrong id {nbk_id}")
                elif citation.source == CitationSource.PUBMED_CENTRAL:
                    pmc_id = citation.data_json.get("PMC")
                    if pmc_id != citation.index:
                        bad_citations.append(citation)
                        print(f"{citation} has data from wrong id {pmc_id}")
                    elif not citation.authors and citation.data_json.get("CN"):
                        bad_citations.append(citation)
                        print(f"{citation} didn't have an author but did have corporate authors")
                    elif not citation.year and citation.data_json.get("DP"):
                        bad_citations.append(citation)
                        year = citation.data_json.get("DP")
                        print(f"{citation} didn't have year - maybe it can be parsed from {year}")

                elif citation.source == CitationSource.PUBMED:
                    pmid = citation.data_json.get("PMID")
                    if pmid != citation.index:
                        bad_citations.append(citation)
                        print(f"{citation} has data from wrong id {pmid}")
                    elif not citation.authors and citation.data_json.get("CN"):
                        bad_citations.append(citation)
                        print(f"{citation} didn't have an author but did have corporate authors")
                    elif not citation.year and citation.data_json.get("DP"):
                        bad_citations.append(citation)
                        year = citation.data_json.get("DP")
                        print(f"{citation} didn't have year - maybe it can be parsed from {year}")

            for citation in bad_citations:
                citation.data_json = None
                citation.last_loaded = None
                citation.save()
            CitationFetchRequest.fetch_all_now(bad_citations, cache_age=timedelta(seconds=0))
