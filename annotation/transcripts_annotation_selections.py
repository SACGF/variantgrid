import operator
from functools import cached_property
from typing import Any, Optional

from django.conf import settings
from django.forms.models import model_to_dict
from django.utils.timesince import timesince

from annotation.models import VEPSkippedReason, AnnotationStatus
from annotation.models.models import VariantAnnotation, AnnotationVersion, \
    InvalidAnnotationVersionError, VariantTranscriptAnnotation, AnnotationRun
from genes.hgvs import HGVSMatcher, HGVSException
from genes.models import TranscriptVersion, GnomADGeneConstraint, Transcript
from genes.models_enums import AnnotationConsortium
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


class VariantTranscriptSelections:
    ENSEMBL_TRANSCRIPT = "ensembl_transcript_accession"
    GENE_SYMBOL = "gene_symbol"
    GNOMAD_GENE_CONSTRAINT_OE_LOF_SUMMARY = "gnomad_gene_constraint_oe_lof_summary"
    GNOMAD_GENE_CONSTRAINT_METHOD = "gnomad_gene_constraint_method"
    GNOMAD_GENE_CONSTRAINT_URL = "gnomad_gene_constraint_url"
    REFSEQ_TRANSCRIPT = "refseq_transcript_accession"
    REPRESENTATIVE = "representative"
    NO_TRANSCRIPT = ''  # Can't be None as needs to be sortable

    def __init__(self, variant: Variant,
                 genome_build: GenomeBuild, annotation_version=None,
                 hide_other_annotation_consortium_transcripts=True):
        """ annotation_version defaults to latest version """

        if annotation_version is None:
            annotation_version = AnnotationVersion.latest(genome_build)

        self.genome_build = genome_build
        self.annotation_consortium = annotation_version.variant_annotation_version.annotation_consortium
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            self.other_annotation_consortium = AnnotationConsortium.ENSEMBL
        else:
            self.other_annotation_consortium = AnnotationConsortium.REFSEQ
        self.hide_other_annotation_consortium_transcripts = hide_other_annotation_consortium_transcripts
        self.variant_annotation: Optional[VariantAnnotation] = None
        self.gene_annotations = {}
        self.transcript_data = []  # See docstring in _populate for details
        self.warning_messages = []
        self.error_messages = []
        self.initial_transcript_id = self.NO_TRANSCRIPT
        self._populate(variant, annotation_version)
        self.other_annotation_consortium_transcripts_warning = None  # set in _add_other_annotation_consortium_transcripts
        self._add_other_annotation_consortium_transcripts(variant)
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            sort_order = [self.REFSEQ_TRANSCRIPT]
        else:
            sort_order = [self.ENSEMBL_TRANSCRIPT]
        self.transcript_data.sort(key=operator.itemgetter(*sort_order), reverse=True)
        self._set_initial_and_selected()

    def get_annotation_consortium_display(self):
        return AnnotationConsortium(self.annotation_consortium).label

    def get_other_annotation_consortium_display(self):
        return AnnotationConsortium(self.other_annotation_consortium).label

    @cached_property
    def variant_transcript_annotations_dict(self):
        return {d["transcript_id"]: d for d in self.transcript_data}

    def get_transcript_annotation(self, transcript_version: TranscriptVersion) -> dict[str, Any]:
        """ Try looking up Transcript accession or failing that, ID """
        vta = self.variant_transcript_annotations_dict
        try:
            return vta[transcript_version.accession]
        except KeyError:
            return vta[transcript_version.transcript_id]

    def _ac_key(self, annotation_consortium):
        if annotation_consortium == AnnotationConsortium.REFSEQ:
            ac_key = self.REFSEQ_TRANSCRIPT
        else:
            ac_key = self.ENSEMBL_TRANSCRIPT
        return ac_key

    def _populate(self, variant, annotation_version):
        vav = annotation_version.variant_annotation_version
        try:
            self.variant_annotation = variant.variantannotation_set.get(version=vav)
            if self.variant_annotation.vep_skipped_reason:
                annotation_error = "Unable to annotate variant"
                if self.variant_annotation.vep_skipped_reason != VEPSkippedReason.UNKNOWN:
                    annotation_error += ": " + self.variant_annotation.get_vep_skipped_reason_display()
                self.error_messages.append(annotation_error)

            representative_transcript = self.variant_annotation.transcript
            transcripts_list = list(variant.varianttranscriptannotation_set.filter(version=vav))

            for vsta in transcripts_list:
                t_data = self._get_transcript_data(vsta, representative_transcript)
                self.transcript_data.append(t_data)

                if vsta.gene_id and vsta.gene_id not in self.gene_annotations:
                    version_id = annotation_version.gene_annotation_version_id
                    gene_annotation = vsta.gene.geneannotation_set.filter(version_id=version_id).first()
                    if gene_annotation:
                        self.gene_annotations[vsta.gene_id] = model_to_dict(gene_annotation)
        except VariantAnnotation.DoesNotExist:
            # Probably due to variant being annotated - in that case show a warning message
            if ar := AnnotationRun.get_for_variant(variant, self.genome_build):
                if ar.status in (AnnotationStatus.ERROR, AnnotationStatus.FINISHED):
                    msg = f"Annotation for variant {variant} failed (Annotation run: {ar.get_status_display()})"
                    self.error_messages.append(msg)
                else:
                    msg = f"This variant has not yet been annotated. Last status: {ar.get_status_display()} ({timesince(ar.modified)} ago)"
                    self.warning_messages.append(msg)
            else:
                # Some other failure
                msg = f"Couldn't find VariantAnnotation for {variant}, VariantAnnotationVersion {vav}"
                self.error_messages.append(msg)
        except InvalidAnnotationVersionError as e:
            self.error_messages.append(str(e))

    def _get_transcript_data(self, obj: VariantTranscriptAnnotation, representative_transcript: Transcript) -> dict:
        data = model_to_dict(obj)
        # Use nice values if available
        for f in data:
            try:
                data[f] = getattr(obj, f"get_{f}_display")()
            except AttributeError:
                pass

        # Split/clean aggregate fields
        VEP_JOINED_FIELDS = ["domains"]
        for field in VEP_JOINED_FIELDS:
            field_value = data[field]
            if field_value:
                data[field] = field_value.replace("&", ", ")

        # These always need to be set for sorting
        data[self.REFSEQ_TRANSCRIPT] = self.NO_TRANSCRIPT
        data[self.ENSEMBL_TRANSCRIPT] = self.NO_TRANSCRIPT
        if obj.transcript:
            ac_key = self._ac_key(obj.transcript.annotation_consortium)
            data[ac_key] = obj.transcript_accession

        data[self.GENE_SYMBOL] = obj.symbol
        data["gene_id"] = obj.gene_id

        transcript_id = self.NO_TRANSCRIPT
        for col in settings.VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES:
            if transcript_id := data.get(col):
                break
        data["transcript_id"] = transcript_id

        transcript_canonical_score = 0
        if obj.transcript_version:
            data["tags"] = obj.transcript_version.tags
            data["protein_length"] = obj.transcript_version.protein_length
            transcript_canonical_score = obj.transcript_version.canonical_score
            ggc, ggc_method, ggc_url = GnomADGeneConstraint.get_for_transcript_version_with_method_and_url(obj.transcript_version)
            if ggc:
                data[self.GNOMAD_GENE_CONSTRAINT_OE_LOF_SUMMARY] = ggc.lof_oe_summary
                data[self.GNOMAD_GENE_CONSTRAINT_METHOD] = ggc_method
                data[self.GNOMAD_GENE_CONSTRAINT_URL] = ggc_url

        data["selected"] = False  # Will be set in _set_initial_and_selected
        data["canonical_score"] = transcript_canonical_score

        is_representative = bool(representative_transcript and representative_transcript == obj.transcript)
        data[self.REPRESENTATIVE] = is_representative
        return data

    def _set_initial_and_selected(self):
        # Initial set to what is passed in, otherwise use representative

        if self.transcript_data:
            # Use rep transcript to fall back on ties
            key_func = operator.itemgetter("canonical_score", self.REPRESENTATIVE)
            td_canonical = sorted(self.transcript_data, key=key_func, reverse=True)
            td_canonical = next(iter(td_canonical))
            selected_transcript_data = None
            if settings.VARIANT_TRANSCRIPT_USE_TRANSCRIPT_CANONICAL and td_canonical["canonical_score"]:
                selected_transcript_data = td_canonical
            else:
                # Use representative (VEP pick) that is used as VariantAnnotation
                for data in self.transcript_data:
                    if data[self.REPRESENTATIVE]:
                        selected_transcript_data = data
                        break

            if selected_transcript_data:
                selected_transcript_data["selected"] = True
                self.initial_transcript_id = selected_transcript_data["transcript_id"]

    def _add_other_annotation_consortium_transcripts(self, variant: Variant):
        """ VariantAnnotation is populated from VEP as either RefSeq or Ensembl
            Whichever one is used will have molecular consequences etc.

            We may have transcripts from consortium not used, so can add those. """
        ac_key = self._ac_key(self.annotation_consortium)
        other_ac_key = self._ac_key(self.other_annotation_consortium)

        svlen = variant.svlen or 0
        if abs(svlen) > settings.HGVS_MAX_SEQUENCE_LENGTH:
            other_ac = AnnotationConsortium(self.other_annotation_consortium)
            self.warning_messages.append(f"Not adding {other_ac.label} transcripts due to length")
            return

        gene_symbols = set()
        for ga in self.gene_annotations.values():
            if gene_symbol := ga.get("hgnc_symbol"):
                gene_symbols.add(gene_symbol)

        existing_other_transcripts = set()
        for td in self.transcript_data:
            if other_transcript_accession := td.get(other_ac_key):
                existing_other_transcripts.add(other_transcript_accession)
            # Just in case it's not the same as HGNC above
            gene_symbol = td.get("symbol")
            if gene_symbol:
                gene_symbols.add(gene_symbol)

        hgvs_matcher = HGVSMatcher(self.genome_build)
        kwargs = {
            "transcript__annotation_consortium": self.other_annotation_consortium,
            "genome_build": self.genome_build,
            "gene_version__gene_symbol_id__in": gene_symbols,
        }

        # Convert once to explicit, then pass this around
        variant_coordinate = variant.coordinate.as_external_explicit(self.genome_build)
        has_other_annotation_consortium_transcripts = False
        for transcript_version in TranscriptVersion.objects.filter(**kwargs).order_by("-version"):
            # Don't duplicate ones already available via RefSeq/Ensembl equivalence
            # and only take the highest version we have
            if transcript_version.transcript_id not in existing_other_transcripts:
                existing_other_transcripts.add(transcript_version.transcript_id)
                t_data = {
                    "gene_id": transcript_version.gene_version.gene_id,
                    self.GENE_SYMBOL: transcript_version.gene_version.gene_symbol_id,
                    ac_key: "",
                    other_ac_key: transcript_version.accession,
                    "transcript_id": transcript_version.accession,  # for transcript data keys
                    self.REPRESENTATIVE: False,
                    "consequence": "?",
                    "canonical_score": -1,
                    "hidden": self.hide_other_annotation_consortium_transcripts,
                }
                try:
                    # Transcript data may not be well formed
                    t_data["protein_length"] = transcript_version.protein_length
                except (HGVSException, KeyError):
                    pass

                try:
                    hgvs_variant = hgvs_matcher.variant_coordinate_to_hgvs_variant(variant_coordinate,
                                                                                   transcript_version.accession)
                    t_data["hgvs_c"] = hgvs_variant.format()
                except (HGVSException, KeyError):
                    pass
                self.transcript_data.append(t_data)
                has_other_annotation_consortium_transcripts = True

        if has_other_annotation_consortium_transcripts:
            ac_dict = dict(AnnotationConsortium.choices)
            self.other_annotation_consortium_transcripts_warning = f"""
This system only annotates {ac_dict[self.annotation_consortium]} transcripts. You may select
{ac_dict[self.other_annotation_consortium]} transcripts, but they will not be auto-populated with transcript annotation"""
