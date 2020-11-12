from django.conf import settings
from django.forms.models import model_to_dict
from lazy import lazy
import operator

from annotation.external_search_terms import get_variant_transcript_annotation_search_terms
from annotation.models import VEPSkippedReason
from annotation.models.models import VariantAnnotation, AnnotationVersion, \
    InvalidAnnotationVersionError, VariantTranscriptAnnotation
from genes.hgvs import HGVSMatcher
from genes.models import TranscriptVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


class VariantTranscriptSelections:
    ENSEMBL_TRANSCRIPT = "ensembl_transcript_accession"
    GENE_SYMBOL = "gene_symbol"
    REFSEQ_TRANSCRIPT = "refseq_transcript_accession"
    REPRESENTATIVE = "representative"
    SEARCH_TERMS = "search_terms"

    def __init__(self, variant: Variant,
                 genome_build: GenomeBuild, annotation_version=None,
                 initial_transcript_id=None, initial_transcript_column=REFSEQ_TRANSCRIPT,
                 add_other_annotation_consortium_transcripts=False):
        """ annotation_version defaults to latest version """

        if annotation_version is None:
            annotation_version = AnnotationVersion.latest(genome_build)

        self.genome_build = genome_build
        self.annotation_consortium = annotation_version.variant_annotation_version.annotation_consortium
        self.variant_annotation = None
        self.gene_annotations = {}
        self.transcript_data = []  # See docstring in _populate for details
        self.error_messages = []
        self.initial_transcript_id = None
        self._populate(variant, annotation_version, initial_transcript_id, initial_transcript_column)
        self.other_annotation_consortium_transcripts_warning = None  # set in _add_other_annotation_consortium_transcripts
        if add_other_annotation_consortium_transcripts:
            self._add_other_annotation_consortium_transcripts(variant)
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            sort_order = [self.REFSEQ_TRANSCRIPT]
        else:
            sort_order = [self.ENSEMBL_TRANSCRIPT]
        self.transcript_data.sort(key=operator.itemgetter(*sort_order), reverse=True)

    def get_annotation_consortium_display(self):
        return dict(AnnotationConsortium.CHOICES).get(self.annotation_consortium)

    @lazy
    def variant_transcript_annotations_dict(self):
        return {d["transcript_id"]: d for d in self.transcript_data}

    def get_transcript_annotation(self, transcript_version: TranscriptVersion):
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

    def _populate(self, variant, annotation_version, initial_transcript_id, initial_transcript_column):

        def get_transcript_data(obj: VariantTranscriptAnnotation, representative_transcript):
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
            data[self.REFSEQ_TRANSCRIPT] = ''
            data[self.ENSEMBL_TRANSCRIPT] = ''
            if obj.transcript:
                ac_key = self._ac_key(obj.transcript.annotation_consortium)
                data[ac_key] = obj.transcript_accession

            data[self.GENE_SYMBOL] = obj.symbol
            data["gene_id"] = obj.gene_id
            data["flagged_pathogenicity"] = obj.flagged_pathogenicity

            for col in settings.VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES:
                transcript_id = data.get(col)
                if transcript_id:
                    data["transcript_id"] = transcript_id
                    break

            is_representative = bool(representative_transcript and representative_transcript == obj.transcript)
            data[self.REPRESENTATIVE] = is_representative

            # Initial set to what is passed in, otherwise use representative
            if initial_transcript_id:
                transcript_id = data[initial_transcript_column]
                selected = initial_transcript_id == transcript_id
            else:
                selected = is_representative
            data["selected"] = selected
            if selected:
                self.initial_transcript_id = data["transcript_id"]

            data[self.SEARCH_TERMS] = get_variant_transcript_annotation_search_terms(obj)
            return data

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
                transcript = vsta.transcript
                if transcript:
                    t_data = get_transcript_data(vsta, representative_transcript)
                    self.transcript_data.append(t_data)

                    if vsta.gene_id not in self.gene_annotations:
                        version = annotation_version.ensembl_gene_annotation_version
                        gene_annotation = vsta.gene.ensemblgeneannotation_set.filter(version=version).first()
                        if gene_annotation:
                            self.gene_annotations[vsta.gene_id] = model_to_dict(gene_annotation)
        except VariantAnnotation.DoesNotExist:
            msg = f"Couldn't find VariantAnnotation for {variant}, VariantAnnotationVersion {vav}"
            self.error_messages.append(msg)
        except InvalidAnnotationVersionError as e:
            self.error_messages.append(str(e))

    def _add_other_annotation_consortium_transcripts(self, variant):
        """ VariantAnnotation is populated from VEP as either RefSeq or Ensembl
            Whichever one is used will have molecular consequences etc.

            We may have transcripts from consortium not used, so can add those. """
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            other_annotation_consortium = AnnotationConsortium.ENSEMBL
        else:
            other_annotation_consortium = AnnotationConsortium.REFSEQ

        ac_key = self._ac_key(self.annotation_consortium)
        other_ac_key = self._ac_key(other_annotation_consortium)

        gene_symbols = set()
        for ga in self.gene_annotations.values():
            gene_symbol = ga.get("hgnc_symbol")
            if gene_symbol:
                gene_symbols.add(gene_symbol)

        existing_other_transcripts = set()
        for td in self.transcript_data:
            other_transcript_accession = td.get(other_ac_key)
            if other_transcript_accession:
                existing_other_transcripts.add(other_transcript_accession)
            # Just in case it's not the same as HGNC above
            gene_symbol = td.get("symbol")
            if gene_symbol:
                gene_symbols.add(gene_symbol)

        hgvs_matcher = HGVSMatcher(self.genome_build)
        kwargs = {
            "transcript__annotation_consortium": other_annotation_consortium,
            "genome_build": self.genome_build,
            "gene_version__gene_symbol_id__in": gene_symbols,
        }

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
                    "consequence": "?"
                }
                try:
                    t_data["hgvs_c"] = hgvs_matcher.variant_to_c_hgvs(variant, transcript_version.accession)
                except (ValueError, KeyError):
                    pass
                self.transcript_data.append(t_data)
                has_other_annotation_consortium_transcripts = True

        if has_other_annotation_consortium_transcripts:
            ac_dict = dict(AnnotationConsortium.CHOICES)
            self.other_annotation_consortium_transcripts_warning = f"""
This system only annotates {ac_dict[self.annotation_consortium]} transcripts. You may select
{ac_dict[other_annotation_consortium]} transcripts, but they will not be auto-populated with transcript annotation"""
