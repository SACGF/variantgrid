import itertools

from django.conf import settings
from django.contrib.sites.models import Site
from typing import List, Tuple, Iterable, Optional

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.citations import get_citations
from annotation.models.damage_enums import FATHMMPrediction, \
    MutationTasterPrediction, Polyphen2Prediction, SIFTPrediction, \
    MutationAssessorPrediction
from annotation.models.models import VariantAnnotation, AnnotationVersion, GeneSymbolPubMedCount
from annotation.models.models_enums import VariantClass, ClinVarReviewStatus
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import VEP_SEPARATOR
from genes.hgvs import HGVSMatcher
from genes.models import TranscriptVersion, GnomADGeneConstraint
from genes.models_enums import AnnotationConsortium
from library.django_utils import get_choices_formatter
from library.log_utils import log_traceback
from seqauto.models import get_20x_gene_coverage
from snpdb.clingen_allele import get_clingen_allele_for_variant, ClinGenAlleleAPIException
from snpdb.models import Variant
from snpdb.models.models_clingen_allele import ClinGenAllele
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_enums import ColumnAnnotationLevel
from classification.enums import SubmissionSource, \
    SpecialEKeys
from classification.models.evidence_key import EvidenceKeyMap


AUTOPOPULATE_MERGE_KEYS = {SpecialEKeys.LITERATURE}


class AutopopulateData:

    def __init__(self, name: str):
        self.name = name
        self._data = {}
        self.message = None
        self.linked = []
        self.annotation_version = None
        self._all_data = {}

    def __contains__(self, key):
        if key in self._data:
            return True
        for other_data in self.linked:
            if key in other_data:
                return True
        return False

    def __getitem__(self, key: str, dvalue=None):
        return self.data.get(key, dvalue)

    def __setitem__(self, key: str, value):
        if value is not None:
            self._data[key] = value
            self._all_data = None

    def update(self, other: 'AutopopulateData'):
        if type(other) != type(self):
            raise ValueError(f'Can only update with another AutopopulateData, not {type(other)}')
        self.linked.append(other)
        self._all_data = None

    @property
    def data_direct(self) -> dict:
        """
        Returns data that is directly collected as part of this collection (use if you called flatten)
        """
        return self._data or {}

    @property
    def data(self) -> dict:
        if not self._all_data:
            all_data = {}
            all_data.update(self._data)
            for other_data in self.linked:
                for key, value in other_data.data.items():
                    if key in all_data and key in AUTOPOPULATE_MERGE_KEYS:
                        all_data[key] = '\n'.join([all_data[key], value])
                    else:
                        all_data[key] = value
            self._all_data = all_data
        return self._all_data

    @property
    def summary(self) -> str:
        sums = []
        ekeys = EvidenceKeyMap()
        flattened = self.flatten()
        sums.append(f'The following fields were auto-populated by {Site.objects.get_current().name}')
        for ap in flattened:
            key_data = [ekeys.get(key).pretty_label for key in ap._data.keys()]
            if key_data:
                key_list = ', '.join(key_data)
                sum_text = f'{ap.name} : {key_list}'
                if ap.message:
                    sum_text = sum_text + '\n*' + ap.message
                sums.append(sum_text)
        return '\n\n'.join(sums)

    def _flatten(self, evidence_list: List['AutopopulateData']):
        found = False
        for check_me in evidence_list:
            if check_me.name == self.name:
                check_me._data.update(self._data)
                found = True
                break
        if not found:
            evidence_list.append(self)

        for other_data in self.linked:
            other_data._flatten(evidence_list)

    def flatten(self) -> List['AutopopulateData']:
        evidence_list = []
        self._flatten(evidence_list)
        return evidence_list


def remove_leading_underscore(value):
    if value and value.startswith('_'):
        return value[1:]
    return value


def clinvar_pipe_formatter(value):
    if value:
        value = ', '.join([remove_leading_underscore(v).replace('_', ' ') for v in value.split('|')])
    return value


def vep_multi_to_list(value):
    return value.split(VEP_SEPARATOR)


def domain_to_pfam(value):
    domains_components = VariantAnnotation.get_domains_components(value)
    return ", ".join(domains_components["Pfam_domain"])


def pubmed_formatter(value):
    return "\n".join([f"PMID:{p}" for p in value.split(VEP_SEPARATOR)])


def ekey_from_vg_column_formatters():
    # Ekey options automatically split "," into lists - so need to use pipes
    clinvar_review_status_to_vcf_dict = {v: k.replace(",", "|") for k, v in ClinVarReviewStatus.VCF_MAPPINGS.items()}
    # key = EKey.key
    return {
        "clinvar_clinical_significance": clinvar_pipe_formatter,
        "clinvar_clinical_sources": clinvar_pipe_formatter,
        "clinvar_disease_database_name": clinvar_pipe_formatter,
        "clinvar_preferred_disease_name": clinvar_pipe_formatter,
        "clinvar_review_status": lambda crs: clinvar_review_status_to_vcf_dict.get(crs),
        "fathmm_pred_most_damaging": get_choices_formatter(FATHMMPrediction.CHOICES),
        "molecular_consequence": vep_multi_to_list,
        'mutation_assessor': get_choices_formatter(MutationAssessorPrediction.CHOICES),
        'mutation_taster': get_choices_formatter(MutationTasterPrediction.CHOICES),
        'polyphen2': get_choices_formatter(Polyphen2Prediction.CHOICES),
        'pfam_protein_domain': domain_to_pfam,
        "literature": pubmed_formatter,
        "sift": get_choices_formatter(SIFTPrediction.CHOICES),
        "variant_class": get_choices_formatter(VariantClass.choices),
    }


def get_clingen_allele_and_evidence_value_for_variant(genome_build: GenomeBuild, variant: Variant) -> Tuple[ClinGenAllele, str, str]:
    """ returns (clingen_allele, evidence_value) """
    message = None
    try:
        clingen_allele = get_clingen_allele_for_variant(genome_build, variant)
        evidence_value = str(clingen_allele)
    except ClinGenAlleleAPIException as cge:
        clingen_allele = None
        evidence_value = None
        message = str(cge)

    return clingen_allele, evidence_value, message


def get_evidence_fields_for_variant(genome_build: GenomeBuild, variant: Variant,
                                    refseq_transcript_accession, ensembl_transcript_accession,
                                    evidence_keys_list: list, annotation_version: AnnotationVersion) -> AutopopulateData:
    """ annotation_version is optional (defaults to latest for genome build) """

    data = AutopopulateData("basic variant")
    hgvs_matcher = HGVSMatcher(genome_build=genome_build)
    clingen_allele = None
    if variant:
        data[SpecialEKeys.G_HGVS] = hgvs_matcher.variant_to_g_hgvs(variant)

        clingen_allele, evidence_value, message = get_clingen_allele_and_evidence_value_for_variant(genome_build, variant)
        if message:
            data.message = message
        data[SpecialEKeys.CLINGEN_ALLELE_ID] = evidence_value
        data[SpecialEKeys.VARIANT_COORDINATE] = {
            'value': variant.full_string,
            'immutable': SubmissionSource.VARIANT_GRID
        }

    # Directly copy over fields according to EvidenceKey.variantgrid_column

    evidence_transcript_columns = {}  # These need to be copied off TranscriptVersion object
    evidence_variant_columns = {}  # Copied from query
    for evidence_key in evidence_keys_list:
        variantgrid_column = evidence_key.variantgrid_column
        if variantgrid_column:
            if variantgrid_column.annotation_level == ColumnAnnotationLevel.TRANSCRIPT_LEVEL:
                evidence = evidence_transcript_columns
            else:
                evidence = evidence_variant_columns

            evidence[evidence_key.key] = {'col': variantgrid_column.variant_column,
                                          'immutable': evidence_key.immutable}

    ekey_formatters = ekey_from_vg_column_formatters()
    # Populate from query 1st so can be overwritten
    if variant:
        data.update(get_evidence_fields_from_variant_query(variant, evidence_variant_columns, ekey_formatters, annotation_version))

    if refseq_transcript_accession or ensembl_transcript_accession:
        transcript_values = {"refseq_transcript_accession": refseq_transcript_accession,
                             "ensembl_transcript_accession": ensembl_transcript_accession}

        data.update(get_evidence_fields_from_transcript_data(genome_build, variant, clingen_allele, hgvs_matcher,
                                                             transcript_values, evidence_transcript_columns,
                                                             ekey_formatters, annotation_version))
    # Calculated fields
    try:
        clinvar = variant.clinvar_set.get(version=annotation_version.clinvar_version)

        try:
            clinvar_literature = get_literature(clinvar)
        except Exception:  # HTTPError but not sure where from (biopython urllib2?)
            clinvar_literature = "Warning: Error retrieving ClinVar citations abstract"

        # May already have existing literature from 'pubmed'
        # Autopopulate joins evidence together
        if clinvar_literature:
            data[SpecialEKeys.LITERATURE] = clinvar_literature
    except:
        pass

    return data


def transcript_autopopulate(transcript_version: TranscriptVersion):
    ac_source = transcript_version.transcript.get_annotation_consortium_display()
    return AutopopulateData(f"{ac_source} transcript")


def get_evidence_fields_from_transcript_data(
        genome_build: GenomeBuild,
        variant: Variant,
        clingen_allele: ClinGenAllele,
        hgvs_matcher: HGVSMatcher,
        transcript_values,
        evidence_transcript_columns,
        ekey_formatters,
        annotation_version) -> AutopopulateData:

    CONSORTIUM_EKEYS = {
        AnnotationConsortium.REFSEQ: [
            "refseq_transcript_accession",
            SpecialEKeys.REFSEQ_TRANSCRIPT_ID,
            SpecialEKeys.ENTREZ_GENE_ID
        ],
        AnnotationConsortium.ENSEMBL: [
            "ensembl_transcript_accession",
            SpecialEKeys.ENSEMBL_TRANSCRIPT_ID,
            SpecialEKeys.ENSEMBL_GENE_ID,
        ]
    }

    data = AutopopulateData("transcript")
    transcript_objects = {}
    for annotation_consortium, (transcript_key, transcript_ekey, gene_ekey) in CONSORTIUM_EKEYS.items():
        transcript_id = transcript_values[transcript_key]
        if transcript_id:
            try:
                transcript_version = TranscriptVersion.get(transcript_id, genome_build, annotation_consortium)
                transcript_objects[transcript_key] = transcript_version
                t_data = transcript_autopopulate(transcript_version)
                t_data[transcript_ekey] = transcript_version.accession
                t_data[gene_ekey] = transcript_version.gene_version.accession
                data.update(t_data)
            except TranscriptVersion.DoesNotExist:
                pass

    transcript_version = TranscriptVersion.get_preferred_transcript(transcript_objects)
    if transcript_version:
        data.update(get_evidence_fields_from_preferred_transcript(genome_build, variant, clingen_allele, hgvs_matcher,
                                                                  transcript_version, evidence_transcript_columns,
                                                                  ekey_formatters, annotation_version))
    return data


def get_evidence_fields_from_preferred_transcript(
        genome_build: GenomeBuild,
        variant: Variant,
        clingen_allele: ClinGenAllele,
        hgvs_matcher: HGVSMatcher,
        transcript_version: TranscriptVersion,
        evidence_transcript_columns,
        ekey_formatters,
        annotation_version) -> AutopopulateData:

    data = transcript_autopopulate(transcript_version)
    data[SpecialEKeys.GENE_SYMBOL] = transcript_version.gene_version.gene_symbol_id

    # Populate from TranscriptVersion data 1st (so can overwrite later)
    if variant:
        try:
            vts = VariantTranscriptSelections(variant, genome_build, annotation_version)
            transcript_data = vts.get_transcript_annotation(transcript_version)

            for evidence_key, transcript_config in evidence_transcript_columns.items():
                variant_column = transcript_config['col']
                immutable = transcript_config['immutable']
                # Getting out of dict directly, not joining via queryset
                transcript_column = variant_column.replace("variantannotation__", "")
                value = transcript_data.get(transcript_column)
                if value is not None:
                    set_evidence(data, evidence_key, value, immutable, ekey_formatters)

            gene_symbol = transcript_version.gene_version.gene_symbol
            data[SpecialEKeys.INTERNAL_SAMPLES_20X_COVERAGE] = get_20x_gene_coverage(gene_symbol)

            phastcons_dict = {
                "phastcons_30_way_mammalian": "30 way mammalian",
                "phastcons_46_way_mammalian": "46 way mammalian",
                "phastcons_100_way_vertebrate": "100 way vertebrate",
            }
            data[SpecialEKeys.PHASTCONS] = get_set_fields_summary(transcript_data, phastcons_dict, phastcons_dict)
            phylop_dict = {
                "phylop_30_way_mammalian": "30 way mammalian",
                "phylop_46_way_mammalian": "46 way mammalian",
                "phylop_100_way_vertebrate": "100 way vertebrate",
            }
            data[SpecialEKeys.PHYLOP] = get_set_fields_summary(transcript_data, phylop_dict, phylop_dict)
            if variant_annotation := vts.variant_annotation:
                data[SpecialEKeys.SEARCH_TERMS] = variant_annotation.get_search_terms()
                if settings.ANNOTATION_PUBMED_SEARCH_TERMS_ENABLED:
                    data[SpecialEKeys.PUBMED_SEARCH_TERMS] = variant_annotation.get_pubmed_search_terms()
        except:
            log_traceback()

    try:
        c_hgvs = hgvs_matcher.variant_to_c_hgvs_parts(variant, transcript_version.accession)
        if c_hgvs:
            data[SpecialEKeys.C_HGVS] = c_hgvs.full_c_hgvs
    except Exception as e:
        value_obj = {}
        data.message = 'Could not parse HGVS value %s' % str(e)
        data[SpecialEKeys.C_HGVS] = value_obj

    if clingen_allele:
        p_hgvs = clingen_allele.get_p_hgvs(transcript_version.accession, match_version=False)
    else:
        p_hgvs = None
    data[SpecialEKeys.P_HGVS] = p_hgvs

    gene_symbol_id = transcript_version.gene_version.gene_symbol_id
    gnomad_oe_lof_summary = get_gnomad_oe_lof_summary(transcript_version)
    if gnomad_oe_lof_summary:
        data[SpecialEKeys.GNOMAD_OE_LOF] = gnomad_oe_lof_summary

    gs_count = GeneSymbolPubMedCount.get_for_gene_symbol(gene_symbol_id)
    data[SpecialEKeys.PUBMED_GENE_SEARCH_COUNT] = {"value": gs_count.count,
                                                   "note": f"Retrieved {gs_count.modified.date()}"}
    return data


def set_evidence(data: AutopopulateData, evidence_key, value, immutable, ekey_formatters):
    formatter = ekey_formatters.get(evidence_key)
    if formatter:
        value = formatter(value)

    if immutable:
        data[evidence_key] = {'value': value, 'immutable': SubmissionSource.VARIANT_GRID}
    else:
        data[evidence_key] = value


def _get_mastermind_summary(variant_values: dict) -> Optional[str]:
    mastermind_summary: Optional[str] = None
    if mastermind_mmid3 := variant_values.get("variantannotation__mastermind_mmid3"):
        mastermind_fields = []
        for field, label in VariantAnnotation.MASTERMIND_FIELDS.items():
            value = variant_values.get(f"variantannotation__{field}")
            mastermind_fields.append(f"{label}: {value}")
        count_summary = ", ".join(mastermind_fields[:-1])
        mastermind_summary = f"Mastermind: {VariantAnnotation.get_mastermind_url(mastermind_mmid3)}\n"
        mastermind_summary += f"Search results: {count_summary}\n"
        mastermind_summary += mastermind_fields[-1]  # MMID3
    return mastermind_summary


def _get_spliceai_summary(variant_values: dict):
    """ Summarise like AG: 0.53@-5, AL: 0.1@35, DG: 0.0@-5, DL: 0.1@17 """
    spliceai_summary = None
    if variant_values.get("variantannotation__spliceai_pred_dp_ag") is not None:
        spliceai_template = []
        for k, (ds, dp) in VariantAnnotation.SPLICEAI_DS_DP.items():
            spliceai_template.append(f"{k}: %(variantannotation__{ds})f@%(variantannotation__{dp})d")
        spliceai_summary = ", ".join(spliceai_template) % variant_values
    return spliceai_summary


def get_evidence_fields_from_variant_query(
        variant: Variant,
        evidence_variant_columns,
        ekey_formatters,
        annotation_version) -> AutopopulateData:
    # Run a Variant query against the appropriate annotation versions etc
    # Then pull out values and populate
    data = AutopopulateData("variant annotations")
    # Retrieve fields for Mastermind/SpliceAI summaries
    va_fields_for_summaries = []
    va_fields_for_summaries.extend(VariantAnnotation.MASTERMIND_FIELDS)
    va_fields_for_summaries.extend(itertools.chain.from_iterable(VariantAnnotation.SPLICEAI_DS_DP.values()))

    qs = get_variant_queryset_for_annotation_version(annotation_version=annotation_version)
    qs = qs.filter(pk=variant.pk)
    columns = {f"variantannotation__{f}" for f in va_fields_for_summaries}
    columns.update([e['col'] for e in evidence_variant_columns.values()])
    values = qs.values(*columns).get()

    for evidence_key, variant_data in evidence_variant_columns.items():
        variant_column = variant_data['col']
        immutable = variant_data['immutable']
        value = values.get(variant_column)
        set_evidence(data, evidence_key, value, immutable, ekey_formatters)

    # VariantGridColumn 'pubmed' is transcript level so not set above
    data[SpecialEKeys.LITERATURE] = _get_mastermind_summary(values)
    data[SpecialEKeys.SPLICEAI] = _get_spliceai_summary(values)
    return data


def get_literature(clinvar):
    literature_references = None
    if clinvar:
        clinvar_citations = clinvar.get_citations()
        clinvar_citation_text_rows = []
        for citation in get_citations(clinvar_citations):
            citation_parts = [citation.authors_short]
            if citation.year:
                citation_parts.append(f"({citation.year})")
            if citation.journal_short:
                citation_parts.append(f"{citation.journal_short},")
            citation_parts.append(citation.title)
            citation_parts = [x for x in citation_parts if x is not None]
            cite_text = " ".join(citation_parts)
            clinvar_citation_text_rows.append(cite_text)
            clinvar_citation_text_rows.append(f"{citation.source}: {citation.citation_id}")
            clinvar_citation_text_rows.append('')

        if clinvar_citation_text_rows:
            literature_references = '\n'.join(clinvar_citation_text_rows)

    return literature_references


def get_gnomad_oe_lof_summary(transcript_version: TranscriptVersion):
    """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
        possible (transcript, gene, then symbol) """
    qs = GnomADGeneConstraint.objects.all()

    ggc = None
    if transcript_version.transcript.annotation_consortium == AnnotationConsortium.ENSEMBL:
        # May be able to get via transcript or gene
        ggc = qs.filter(transcript=transcript_version.transcript).first()
        if ggc is None:
            ggc = qs.filter(gene=transcript_version.gene).first()

    if ggc is None:  # Fall back on matching via symbol
        ggc = qs.filter(gene_symbol=transcript_version.gene_version.gene_symbol).first()

    oe_lof_summary = None
    if ggc:
        oe_lof_summary = ggc.oe_lof_summary
    return oe_lof_summary


def get_set_fields_summary(d: dict, fields: Iterable[str], field_labels: dict = None):
    values = {}
    for f in fields:
        v = d.get(f)
        if v:
            if field_labels:
                label = field_labels.get(f, f)
            else:
                label = f
            values[label] = v

    return ", ".join([f"{k}: {v}" for k, v in values.items()])
