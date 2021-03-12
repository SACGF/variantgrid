from django.conf import settings
from django.contrib.sites.models import Site
from django.urls.base import reverse
import urllib

from library.django_utils import get_url_from_view_path
from snpdb.models.models_genome import GenomeBuildContig
from classification.enums.classification_enums import SpecialEKeys
from classification.models.evidence_key import EvidenceKeyMap, EvidenceKey
from classification.models.classification import ClassificationModification
from classification.views.classification_export_utils import ExportFormatter, \
    AlleleGroup, VCFEncoding


class ExportFormatterVCF(ExportFormatter):
    """
    Formats as a VCF with INFO fields as RefSeq Transcript ID, clinical significance and condition as INFO
    """

    def __init__(self, encoding:str = VCFEncoding.BASIC, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.report_keys = [self.ekeys.get(key) for key in [SpecialEKeys.CLINICAL_SIGNIFICANCE]]
        self.db_prefix = settings.VARIANT_VCF_DB_PREFIX
        self.encoding = encoding

    def generate_info_for_key(self, ekey: EvidenceKey):
        return f'##INFO=<ID={ekey.key},Number=.,Type=String,Description="{ekey.pretty_label}">\n'

    def generate_value_for_key(self, ekey: EvidenceKey, record: ClassificationModification):
        return self.vcf_safe(record.get(ekey.key, None))

    def vcf_safe(self, value):
        if value is None:
            value = ''
        elif isinstance(value, list):
            value = '|'.join(value)
        else:
            value = str(value)

        if self.encoding == VCFEncoding.FULL:
            return urllib.parse.quote(value)
        else:
            value = value.replace(' ', '_').replace('\t', '_').replace('\n', '_')
            value = value.replace(';', ':')
            value = value.replace(',', '.')
            return value

    def header(self) -> str:
        out = '##fileformat=VCFv4.1\n'
        out += '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">\n'
        assembly = self.genome_build.name

        for gbc in GenomeBuildContig.objects.filter(genome_build=self.genome_build, contig__in=self.used_contigs).order_by('order').select_related('contig'):
            contig = gbc.contig
            out += f'##contig=<ID={contig.name},length={contig.length},assembly={assembly}>\n'
        # fixme - substitute application name
        out += f'##INFO=<ID=count,Number=1,Type=Integer,Description="Number of records for this variant">\n'
        out += f'##INFO=<ID=url,Number=1,Type=String,Description="Shariant URL record">\n'
        out += f'##INFO=<ID=lab,Number=.,Type=String,Description="Source lab for the variant classifications">\n'
        out += f'##INFO=<ID=chgvs,Number=.,Type=String,Description="The c.hgvs for the variants in {self.genome_build.name}">\n'
        out += f'##INFO=<ID=multiple_clinical_significances,Number=0,Type=Flag,Description="Present if there are multiple clinical significances for the variant">\n'
        out += f'##INFO=<ID=discordant,Number=.,Type=Number,Description="If 1, indicates that the corresponding classification is in discordance">\n'
        out += f'##INFO=<ID=condition,Number=.,Type=String,Description="Condition Under Curation">\n'
        for e_key in self.report_keys:
            out += self.generate_info_for_key(e_key)

        site = Site.objects.get_current()
        out += f'##readme=Generated from {site.name} using VariantGrid technology\n'
        out += f'##readme=For variants with multiple classifications, each info will be listed in corresponding order, e.g. the 2nd variant value is for the same record as the 2nd condition\n'

        if len(self.error_message_ids):
            out += f'##readme=Warning at least {len(self.error_message_ids)} records had issues with their variants and could not be included\n'

        out += ExportFormatter.write_single_row(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], delimiter='\t')
        return out

    def row(self, group: AlleleGroup) -> str:
        data = []
        vcms = group.data

        # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
        target_variant = group.target_variant
        locus = target_variant.locus
        contig = locus.contig
        data.append(contig.name)  # CHROM
        data.append(locus.position)  # POS
        data.append(self.db_prefix + str(target_variant.id))  # ID
        data.append(locus.ref)  # REF
        data.append(target_variant.alt)  # ALT
        data.append('.')  # QUAL
        data.append('PASS')  # FILTER
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
            lab_names.append(self.vcf_safe(record.classification.lab.name))
            cs = record.get(SpecialEKeys.CLINICAL_SIGNIFICANCE, None)
            if cs:
                cs_counts[cs] = cs_counts.get(cs, 0) + 1
            c_hgvs = group.liftover(record)

            c_hgvses.append(self.vcf_safe(c_hgvs.full_c_hgvs) if c_hgvs else '')
            discordances.append('1' if self.is_discordant(record.classification) else '0')
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
            info += ',' .join(parts)

        data.append(info)  # INFO

        return ExportFormatter.write_single_row(data, '\t')

    def filename(self) -> str:
        return self.generate_filename(extension='vcf')
