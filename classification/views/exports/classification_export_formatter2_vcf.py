import urllib
from dataclasses import dataclass
from enum import Enum
from typing import List

from django.conf import settings
from django.contrib.sites.models import Site
from django.http import HttpRequest
from django.urls import reverse

from classification.enums import SpecialEKeys
from classification.models import EvidenceKey, ClassificationModification, EvidenceKeyMap
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2
from library.django_utils import get_url_from_view_path
from library.utils import delimited_row
from snpdb.models import GenomeBuildContig


class FormatDetailsVCFEncodingLevel(str, Enum):
    BASIC = "basic"
    FULL = "full"


@dataclass(frozen=True)
class FormatDetailsVCF:
    encoding_level: FormatDetailsVCFEncodingLevel
    db_prefix: str = settings.VARIANT_VCF_DB_PREFIX

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsVCF':
        encoding_level = FormatDetailsVCFEncodingLevel.BASIC
        if request.query_params.get('encoding') == 'full':
            encoding_level = FormatDetailsVCFEncodingLevel.FULL
        return FormatDetailsVCF(
            encoding_level=encoding_level
        )


@register_classification_exporter("vcf")
class ClassificationExportFormatter2VCF(ClassificationExportFormatter2):
    """
    Exports data in the format that Agilent's Alissa can import it
    """

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsVCF):
        self.format_details = format_details
        self.report_keys = [EvidenceKeyMap.cached_key(key) for key in [SpecialEKeys.CLINICAL_SIGNIFICANCE]]
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatter2VCF':
        classification_filter = ClassificationFilter.from_request(request)
        return ClassificationExportFormatter2VCF(
            classification_filter=classification_filter,
            format_details=FormatDetailsVCF.from_request(request)
        )

    def generate_info_for_key(self, ekey: EvidenceKey):
        return f'##INFO=<ID={ekey.key},Number=.,Type=String,Description="{ekey.pretty_label}">'

    def generate_value_for_key(self, ekey: EvidenceKey, record: ClassificationModification):
        return self.vcf_safe(record.get(ekey.key, None))

    def vcf_safe(self, value):
        if value is None:
            value = ''
        elif isinstance(value, list):
            value = '|'.join(value)
        else:
            value = str(value)

        if self.format_details.encoding_level == FormatDetailsVCFEncodingLevel.FULL:
            return urllib.parse.quote(value)
        else:
            value = value.replace('[', '').replace(']', '')  # square brackets used for org name but might cause issues
            value = value.replace(' ', '_').replace('\t', '_').replace('\n', '_')
            value = value.replace(';', ':')
            value = value.replace(',', '.')
            return value


    def header(self) -> List[str]:
        out = []
        out += [
            '##fileformat=VCFv4.1',
            '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">'
        ]
        genome_build = self.classification_filter.genome_build
        assembly = genome_build.name

        genome_build_str = '37'
        if self.genome_build.is_version(38):
            genome_build_str = '38'

        contig_field = f'classification__allele_info__grch{genome_build_str}__variant__locus__contig'
        contigs = self.classification_filter.cms_qs().values_list(contig_field, flat=True).order_by(contig_field).distinct()

        for gbc in GenomeBuildContig.objects.filter(genome_build=genome_build, contig__in=contigs).order_by('order').select_related(
                'contig'):
            contig = gbc.contig
            out += [f'##contig=<ID={contig.name},length={contig.length},assembly={assembly}>']
        # fixme - substitute application name

        out += [
            f'##INFO=<ID=count,Number=1,Type=Integer,Description="Number of records for this variant">',
            f'##INFO=<ID=url,Number=1,Type=String,Description="Shariant URL record">',
            f'##INFO=<ID=lab,Number=.,Type=String,Description="Source lab for the variant classifications">',
            f'##INFO=<ID=chgvs,Number=.,Type=String,Description="The c.hgvs for the variants in {self.genome_build.name}">',
            f'##INFO=<ID=multiple_clinical_significances,Number=0,Type=Flag,Description="Present if there are multiple clinical significances for the variant">',
            f'##INFO=<ID=discordant,Number=.,Type=Number,Description="If 1, indicates that the corresponding classification is in discordance">',
            f'##INFO=<ID=condition,Number=.,Type=String,Description="Condition Under Curation">'
        ]
        for e_key in self.report_keys:
            out += [self.generate_info_for_key(e_key)]

        site = Site.objects.get_current()
        out += [
            f'##readme=Generated from {site.name} using VariantGrid technology',
            f'##readme=For variants with multiple classifications, each info will be listed in corresponding order, e.g. the 2nd variant value is for the same record as the 2nd condition'
        ]

        # if len(self.error_message_ids):
        #     out += [f'##readme=Warning at least {len(self.error_message_ids)} records had issues with their variants and could not be included']

        out += [delimited_row(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], delimiter='\t')]
        return out


    def row(self, data: AlleleData) -> List[str]:
        if (vcms := data.cms) and (target_variant := data.cached_variant) and (locus := target_variant.locus) and (contig := locus.contig):
            cols = []

            # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
            cols.append(contig.name)  # CHROM
            cols.append(locus.position)  # POS
            cols.append(self.format_details.db_prefix + str(target_variant.id))  # ID
            cols.append(locus.ref)  # REF
            cols.append(target_variant.alt)  # ALT
            cols.append('.')  # QUAL
            cols.append('PASS')  # FILTER
            info = f'count={len(vcms)}'
            url = get_url_from_view_path(reverse("view_allele_from_variant", kwargs={"variant_id": target_variant.id}))
            info += f';url={self.vcf_safe(url)}'

            # add lab names, and also calculate if we have conflicting classifications and include that flag if appropriate
            cs_counts = {}
            lab_names = []
            c_hgvses = []
            discordances = []
            conditions = []

            for record in vcms:
                lab_names.append(self.vcf_safe(str(record.classification.lab)))
                cs = record.get(SpecialEKeys.CLINICAL_SIGNIFICANCE, None)
                if cs:
                    cs_counts[cs] = cs_counts.get(cs, 0) + 1
                c_hgvs = record.classification.allele_info[data.source.genome_build].c_hgvs_full or ''

                c_hgvses.append(c_hgvs)
                discordances.append('1' if self.classification_filter.is_discordant(record) else '0')
                conditions.append(self.vcf_safe(record.condition_text))

            info += f';lab=' + ','.join(lab_names) + ';chgvs=' + ','.join(c_hgvses)
            if len(cs_counts) > 1:
                info += ';multiple_clinical_significances'
            info += f';discordant=' + ','.join(discordances)
            info += f';condition=' + ','.join(conditions)

            # loop through all the keys we're choosing to print out
            for ekey in self.report_keys:
                info += f';{ekey.key}='
                parts = []
                for record in vcms:
                    parts.append(self.generate_value_for_key(ekey, record))
                info += ','.join(parts)

            cols.append(info)  # INFO

            self.row_count += 1
            return [delimited_row(cols, '\t')]
        else:
            return []

    def content_type(self) -> str:
        return "text/plain"

    def extension(self) -> str:
        return "vcf"
