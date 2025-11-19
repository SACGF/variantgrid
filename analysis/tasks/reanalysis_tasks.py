import logging
from collections import defaultdict

from analysis.models import Analysis, Candidate, AnalysisNode
from analysis.tasks.abstract_candidate_search_task import AbstractCandidateSearchTask
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models import AnnotationVersion
from snpdb.models import GenomeBuild, Sample
from variantgrid.celery import app


class ReAnalysisNewAnnotationTask(AbstractCandidateSearchTask):
    def get_candidate_records(self, candidate_search_run):
        # Get out config

        records = []
        for genome_build in GenomeBuild.builds_with_annotation():
            analyses_qs = Analysis.filter_for_user(candidate_search_run.user)
            analyses_qs = analyses_qs.filter(genome_build=genome_build,
                                             visible=True, template_type__isnull=True,
                                             annotation_version__clinvar_version__isnull=False)

            av_latest = AnnotationVersion.latest(genome_build)
            analysis_older_clinvar = analyses_qs.filter(annotation_version__clinvar_version__lt=av_latest.clinvar_version_id)
            # No point looking unless old ones exist
            if analysis_older_clinvar.exists():
                # Get all the samples - and the highest annotation version that they were analysed with
                # What we are looking for are samples where the LAST analysis is too old
                sample_analyses = defaultdict(set)
                nodes_qs = AnalysisNode.filter(analyses__in=analyses_qs, analysisnode_parent__isnull=True)
                for node in nodes_qs.select_related("analysis", "analysis__annotation_version").select_subclasses():
                    for sample in node.get_samples_from_node_only_not_ancestors():
                        sample_analyses[sample].add(node.analysis)

                # Double check to make sure we can see these samples not just analyses
                samples_qs = Sample.filter_for_user(candidate_search_run.user)
                samples_qs = samples_qs.filter(pk__in=sample_analyses)

                # It might be worth sorting them into analysis versions??
                samples_and_analyses_by_annotation_version = defaultdict(list)

                for sample in samples_qs:
                    if analyses := sample_analyses[sample]:
                        # highest annotation version, tiebreaker = last analysed
                        analysis = max(analyses, key=lambda a: (a.annotation_version, a.modified))
                        samples_and_analyses_by_annotation_version[analysis.annotation_version].append((sample, analysis))

                for annotation_version in sorted(samples_and_analyses_by_annotation_version):
                    if annotation_version.clinvar_version >= av_latest.clinvar_version:
                        logging.info("Skipped %s as uses same clinvar version", annotation_version)

                    qs_latest = get_variant_queryset_for_annotation_version(av_latest)
                    cv_patho_qs_latest = qs_latest.filter(clinvar__highest_pathogenicity=5)

                    qs = get_variant_queryset_for_annotation_version(annotation_version)
                    cv_patho_old = qs.filter(clinvar__highest_pathogenicity=5)
                    new_clinvar_patho_qs = cv_patho_qs_latest.difference(cv_patho_old)
                    # Test if there are any?? Skip if not?

                    if not new_clinvar_patho_qs.exists():
                        logging.info("Skipped %s as no new pathogenic clinvar variants", annotation_version)

                    # Faster if list??
                    new_patho_variants = new_clinvar_patho_qs.values_list("pk", flat=True)

                    for (sample, analysis) in samples_and_analyses_by_annotation_version[annotation_version]:
                        print(sample, analysis)
                        sample_qs = sample.get_variant_qs()
                        for variant in sample_qs.filter(pk__in=new_patho_variants):
                            # TODO: could clinvar for this version be more efficiently done by annotating the QS?
                            clinvar = variant.clinvar_set.filter(version=av_latest.clinvar_version).first()
                            notes = "TODO"
                            evidence = {
                            }
                            records.append(Candidate(
                                search_run=candidate_search_run,
                                analysis=analysis,
                                sample=sample,
                                variant=variant,
                                clinvar=clinvar,
                                notes=notes,
                                evidence=evidence,
                            ))

        return records

ReAnalysisNewAnnotationTask = app.register_task(ReAnalysisNewAnnotationTask())
