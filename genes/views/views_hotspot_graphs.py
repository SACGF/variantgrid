import uuid
from functools import cached_property
from typing import Tuple, List, Optional, Dict

from django.conf import settings
from django.db.models.aggregates import Count
from django.shortcuts import get_object_or_404
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from global_login_required import login_not_required

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models.models import AnnotationVersion, VariantAnnotation, VariantAnnotationVersion
from annotation.models.molecular_consequence_enums import MolecularConsequenceColors
from classification.enums import ShareLevel
from classification.models import ClassificationModification, Classification
from genes.models import Transcript, Gene, TranscriptVersion
from library.constants import MINUTE_SECS
from library.utils import segment
from snpdb.models import CohortGenotypeCollection, Cohort, VariantZygosityCountCollection
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant


class HotspotGraphView(TemplateView):
    template_name = "genes/hotspot_graph.html"
    has_graph_filter_toolbar = True

    def _get_title(self, transcript_description) -> str:
        return f"Het/Hom variants for {transcript_description}"

    def _get_y_label(self) -> str:
        return "# samples"

    def _get_variant_queryset(self, transcript_version):
        annotation_version = AnnotationVersion.latest(transcript_version.genome_build)
        qs = get_variant_queryset_for_annotation_version(annotation_version)
        return qs.filter(Variant.get_no_reference_q(),
                         varianttranscriptannotation__transcript_version=transcript_version,
                         varianttranscriptannotation__protein_position__isnull=False)

    def _get_values(self, transcript_version) -> Tuple[int, str, str, str, int]:
        """ :returns hgvs_p, pp, consequence, count, gnomad_af """
        qs = self._get_variant_queryset(transcript_version)
        qs, vzcc = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        qs = qs.filter(**{f"{vzcc.germline_counts_alias}__gt": 0})
        return qs.values_list("varianttranscriptannotation__hgvs_p",
                              "varianttranscriptannotation__protein_position",
                              "varianttranscriptannotation__consequence",
                              "variantannotation__gnomad_af",
                              vzcc.germline_counts_alias)

    @cached_property
    def genome_build(self):
        genome_build_name = self.kwargs["genome_build_name"]
        return GenomeBuild.get_name_or_alias(genome_build_name)

    def _prefer_canonical_with_diff_version(self) -> bool:
        """ Choose canonical transcripts (w/different version) over exact version non-canonical transcripts """
        return settings.VIEW_GENE_HOTSPOT_GRAPH_PREFER_CANONICAL_WITH_DIFF_VERSION

    def _transcript_url(self) -> str:
        return "gene_symbol_transcript_version_hotspot_graph"

    @cached_property
    def _pick_transcripts(self) -> Tuple[Optional[TranscriptVersion], str, List[Tuple[str, bool, str, str]], Dict]:
        """
            The transcript diagram is built via transcript_version.get_domains()
            which uses PfamSequenceIdentifier (linked to transcript and has a version number)

            The protein position of variants come from VariantAnnotation, which means we need to match the transcripts
            used between gene annotation release and PfamSequenceIdentifier

            Sometimes we have to choose between using canonical but bumping a version vs using non-canonical
        """

        transcript_id = self.kwargs.get("transcript_id")
        gene_id = self.kwargs.get("gene_id")
        gene_symbol_id = self.kwargs.get("gene_symbol")
        transcript_accession = self.kwargs.get("transcript_accession")

        transcript_version: Optional[TranscriptVersion] = None
        method = ""
        transcript_options = []
        transcript_urls = {}

        # Using gene annotation release will restrict it to the desired annotation consortium and versions
        vav = VariantAnnotationVersion.latest(self.genome_build)
        if gar := vav.gene_annotation_release:
            if transcript_id:
                lookup = f"Transcript: {transcript_id}"
                transcript = Transcript.objects.get(pk=transcript_id)
                tv_qs = gar.transcript_versions_for_transcript(transcript)
            elif gene_id or gene_symbol_id:
                if gene_id:
                    lookup = f"Gene: {gene_id}"
                    gene = get_object_or_404(Gene, identifier=gene_id)
                    tv_qs = gar.transcript_versions_for_gene(gene)
                else:
                    lookup = f"GeneSymbol: {gene_symbol_id}"
                    tv_qs = gar.transcript_versions_for_symbol(gene_symbol_id)
            else:
                raise ValueError("At least one of 'gene_symbol', 'gene_id' or 'transcript_id' must be in url kwargs")

            # Sort by canonical then choose higher version if any ties
            tv_list = list(sorted(tv_qs, key=lambda tv: (tv.canonical_score, tv.version), reverse=True))
            if tv_list:
                # First, we look for canonical in our gene annotation release
                canonical, non_canonical = segment(tv_list, filter_func=lambda tv: tv.canonical_score)
                if self._prefer_canonical_with_diff_version:
                    SEARCH_ORDER = [
                        ("canonical", canonical, True),
                        ("canonical", canonical, False),
                        ("non-canonical", non_canonical, True),
                        ("non-canonical", non_canonical, False),
                    ]
                else:
                    SEARCH_ORDER = [
                        ("canonical", canonical, True),
                        ("non-canonical", non_canonical, True),
                        ("canonical", canonical, False),
                        ("non-canonical", non_canonical, False),
                    ]

                for _description, tv_list, require_version_match in SEARCH_ORDER:
                    for tv in tv_list:
                        domains, domain_transcript_accession = tv.protein_domains_and_accession
                        if not domains.exists():
                            continue
                        version_match = tv.accession == domain_transcript_accession
                        found = version_match or not require_version_match
                        details = ""
                        if not version_match:
                            details = f"domain: {domain_transcript_accession}"

                        if found and transcript_version is None:
                            if transcript_accession:
                                if tv.accession == transcript_accession:
                                    transcript_version = tv
                                    method = f"Requested transcript accession {transcript_accession}"
                            else:
                                transcript_version = tv
                                method = f"Looked up: {lookup}"

                        # We loop through this twice so take last (so transcript_version has been set)
                        if not require_version_match:
                            active = tv == transcript_version
                            transcript_options.append((tv.accession, active, tv.canonical_tag, details))

                            kwargs = self.kwargs.copy()
                            kwargs["transcript_accession"] = tv.accession
                            transcript_urls[tv.accession] = reverse(self._transcript_url(), kwargs=kwargs)

        return transcript_version, method, transcript_options, transcript_urls

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        transcript_version, method, transcript_options, transcript_urls = self._pick_transcripts
        if transcript_version:
            variant_data = []
            for hgvs_p, protein_position, consequence, gnomad_af, count in self._get_values(transcript_version):
                consequence = MolecularConsequenceColors.CONSEQUENCE_LOOKUPS.get(consequence,
                                                                                 MolecularConsequenceColors.OTHER)
                protein_aa1 = ""
                if hgvs_p:
                    protein_aa3 = hgvs_p.split(":p.")[1]
                    protein_aa1 = VariantAnnotation.amino_acid_3_to_1(protein_aa3)

                pp = VariantAnnotation.protein_position_to_int(protein_position)
                variant_data.append((protein_aa1, pp, consequence, gnomad_af, count))

            domains, domain_transcript_accession = transcript_version.protein_domains_and_accession
            transcript_description = str(transcript_version.accession)
            if transcript_description != domain_transcript_accession:
                transcript_description += f" (prot. domains: {domain_transcript_accession})"

            molecular_consequence_colors = dict(MolecularConsequenceColors.HOTSPOT_COLORS)
            context.update({
                "transcript_version": transcript_version,
                "molecular_consequence_colors": molecular_consequence_colors,
                "domains": list(domains),
                "variant_data": variant_data,
                "uuid": uuid.uuid4(),
                "title": self._get_title(transcript_description),
                "lookup_method": method,
                "transcript_options": transcript_options,
                "transcript_urls": transcript_urls,
                "y_title": self._get_y_label(),
                "has_graph_filter_toolbar": self.has_graph_filter_toolbar,
            })
        return context


class ClassificationsHotspotGraphView(HotspotGraphView):
    has_graph_filter_toolbar = False

    def _get_title(self, transcript_description) -> str:
        return f"Classifications for {transcript_description}"

    def _get_y_label(self) -> str:
        return "# classifications"

    def _get_classifications(self):
        vcm_qs = ClassificationModification.latest_for_user(self.request.user, published=True)
        return Classification.objects.filter(classificationmodification__in=vcm_qs)

    def _get_values(self, transcript_version) -> Tuple[int, str, str, str, int]:
        """ :returns hgvs_p, pp, consequence, count, gnomad_af """
        qs = self._get_variant_queryset(transcript_version)
        vc_qs = self._get_classifications()
        # Need to join through Allele so we get classifications through all genome builds
        vc_column = "variantallele__allele__classification"  # variantallele__allele__classification
        qs = qs.filter(**{vc_column + "__in": vc_qs})
        count_column = "classifications_count"
        qs = qs.values("variantallele__allele").annotate(**{count_column: Count(vc_column)})
        qs = qs.filter(**{count_column + "__gt": 0})
        return qs.values_list("varianttranscriptannotation__hgvs_p",
                              "varianttranscriptannotation__protein_position",
                              "varianttranscriptannotation__consequence",
                              "variantannotation__gnomad_af",
                              count_column)

    def _prefer_canonical_with_diff_version(self) -> bool:
        return settings.VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS_PREFER_CANONICAL_WITH_DIFF_VERSION

    def _transcript_url(self) -> str:
        return "classifications_gene_symbol_transcript_version_hotspot_graph"


class CohortHotspotGraphView(HotspotGraphView):

    @cached_property
    def genome_build(self):
        return self.cohort.genome_build

    @cached_property
    def cohort(self):
        cohort_id = self.kwargs["cohort_id"]
        return Cohort.get_for_user(self.request.user, cohort_id)

    def _get_title(self, transcript_description) -> str:
        return f"{self.cohort} - {transcript_description}"

    def _get_values(self, transcript_version):
        """ :returns hgvs_p, pp, consequence, count, gnomad_af """
        qs = self._get_variant_queryset(transcript_version)
        qs = qs.filter(cohortgenotype__collection=self.cohort.cohort_genotype_collection)
        qs = qs.filter(variantannotation__protein_position__isnull=False)
        qs, count_column = CohortGenotypeCollection.annotate_all_counts(qs)
        qs = qs.filter(**{count_column + "__gt": 0})
        return qs.values_list("varianttranscriptannotation__hgvs_p",
                              "varianttranscriptannotation__protein_position",
                              "varianttranscriptannotation__consequence",
                              "variantannotation__gnomad_af",
                              count_column)


@method_decorator([cache_page(MINUTE_SECS), login_not_required], name='dispatch')
class PublicRUNX1HotspotGraphView(ClassificationsHotspotGraphView):
    """ RUNX1 would like a hotspot graph on the front page - but we don't want to expose
        every hotspot graph obviously, so make a special case one """

    @cached_property
    def genome_build(self):
        return GenomeBuild.get_name_or_alias("GRCh37")  # Hardcoded for server

    @cached_property
    def transcript(self) -> Optional[Transcript]:
        return Transcript.objects.get(pk="ENST00000300305")

    def _get_title(self, transcript_description) -> str:
        title = super()._get_title(transcript_description)
        return title + " (Germline)"

    def _get_classifications(self):
        # Modified from RUNX1_classified_damage - GeneCountType.get_classification_qs()
        kwargs = {"share_level__in": ShareLevel.same_and_higher(ShareLevel.ALL_USERS),
                  "published_evidence__allele_origin__value__in": ['germline', 'likely_germline'],
                  "is_last_published": True}
        vcm_qs = ClassificationModification.objects.filter(**kwargs)
        return Classification.objects.filter(pk__in=vcm_qs.values('classification'))
