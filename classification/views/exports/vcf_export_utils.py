import inspect
import urllib
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional, Iterable, TypeAlias, Protocol
from django.conf import settings
from library.utils import ExportTweak, local_date_str_no_dash, get_decorated_methods
from snpdb.models import Variant, GenomeBuild, GenomeBuildContig, Contig, Allele


"""
This works much like ExportRow @export_column but for VCFs
This is for writing very customised VCFs, so the standard VCF Writer isn't always appropriate
It does NOT currently support adding FORMAT columns
"""


@dataclass(frozen=True)
class VCFExportTweak:
    """
    Used to export only some info fields form an ExportVCF
    """
    categories: dict[str, Any] = None


VCFExportTweak.DEFAULT = VCFExportTweak()


class VCFHeaderNumberSpecial(str, Enum):
    ONE_ALTERNATE_ALLELE = "A"
    ONE_VALUE_EACH_POSSIBLE_ALLELE = "R"
    ONE_VALUE_EACH_POSSIBLE_GENOTYPE = "G"
    UNBOUND = "."


class VCFHeaderClassStandard(str, Enum):
    INFO = "INFO"
    FILTER = "FILTER"
    FORMAT = "FORMAT"


class VCFHeaderType(str, Enum):
    Integer = "Integer"
    String = "String"
    Float = "Float"
    Flag = "Flag"
    Character = "Character"


VCFHeaderNumber: TypeAlias = int | VCFHeaderNumberSpecial


VCFHeaderClass: TypeAlias = str | VCFHeaderClassStandard


class FormatDetailsVCFEncodingLevel(str, Enum):
    BASIC = "basic"
    FULL = "full"


def vcf_value_encode(value, encoding_level: FormatDetailsVCFEncodingLevel = FormatDetailsVCFEncodingLevel.BASIC):
    """ encodes VCF info cell """
    if encoding_level == FormatDetailsVCFEncodingLevel.FULL:
        return urllib.parse.quote(value)
    else:
        value = value.replace('[', '').replace(']', '')  # square brackets used for org name but might cause issues
        value = value.replace(' ', '_').replace('\t', '_').replace('\n', '_')
        value = value.replace(';', ':')
        value = value.replace(',', '.')
        return value


@dataclass(frozen=True)
class VCFHeader:
    header_class: VCFHeaderClass  # INFO, FILTER, FORMAT, CUSTOM
    header_id: Optional[str] = None
    header_type: Optional[VCFHeaderType] = None
    number: Optional[VCFHeaderNumber] = None
    description: Optional[str] = None
    source: Optional[str] = None  # anyone use this?
    custom_attributes: Optional[dict[str, str]] = None
    value: Optional[str] = None

    def __str__(self):
        """
        Makes VCF header string
        """
        header_class = self.header_class.replace(" ", "_")
        if hasattr(header_class, "value"):
            header_class = header_class.value
        if self.value:
            return f"##{header_class}={self.value}"
        parts = [f"ID={self.header_id}"]
        if self.number is not None:
            number = self.number
            if hasattr(number, "value"):
                number = number.value
            parts.append(f"Number={number}")
        if header := self.header_type:
            if hasattr(header, "value"):
                header = header.value
            parts.append(f"Type={header}")
        if self.description is not None:
            description = self.description.replace("\"", "'")
            parts.append(f"Description=\"{description}\"")
        if custom_attributes := self.custom_attributes:
            for key, value in custom_attributes.items():
                parts.append(f"{key}={value}")

        return f"##{header_class}=<{','.join(parts)}>"

    def format_info_value(self, result: Any) -> Optional[str]:
        """
        Assuming the VCFHeader is an INFO field, this will fomrat the value e.g. "classification=3,4"
        Result type, and number thereof is validated
        :param result: a list or single value
        :return: text to be placed directly into the VCF info cell
        """
        if result is None or result == "" or (isinstance(result, list) and len(result) == 0):
            return None
        if self.header_type == VCFHeaderType.Flag:
            if isinstance(result, bool):
                if result:
                    return self.header_id
                else:
                    return None
            else:
                raise ValueError(f"VCF cell {self.header_id} did not return a bool")
        else:
            if not isinstance(result, list):
                result = [result]
            if isinstance(self.number, int):
                if len(result) != self.number:
                    raise ValueError(f"VCF cell {self.header_id} is meant to have {self.number} occurrences but value returned {len(result)})")
            complete_result = ",".join(self._format_info_value_item(item) for item in result)
            return f"{self.header_id}={complete_result}"

    def _format_info_value_item(self, value: Any):
        """
        Formats a single value (could be in a list of multiple values)
        """
        if self.header_type == VCFHeaderType.Float:
            if isinstance(value, float):
                return f"{value:.2f}"  # TODO, allow for definable number of decimal points
            else:
                raise ValueError(f"VCF cell {self.header_id} is meant to be float, received {value}")

        elif self.header_type == VCFHeaderType.Integer:
            if isinstance(value, int):
                return f"{value}"
            else:
                raise ValueError(f"VCF cell {self.header_id} is meant to be int, received {value}")

        elif self.header_type == VCFHeaderType.Character:
            if isinstance(value, str) and len(value) == 1:
                encoded = vcf_value_encode(value)
                if len(encoded) != 1:
                    raise ValueError(f"VCF cell {self.header_id} had VCF unsafe character returned {value}")
                return encoded
            else:
                raise ValueError(f"VCF cell {self.header_id} is meant to be character, received {value}")

        elif self.header_type == VCFHeaderType.String:
            if value is None:
                value = ""
            elif hasattr(value, "value"):
                value = value.value
            return vcf_value_encode(str(value))

        else:
            raise ValueError(f"VCF cell {self.header_id} was defined with an unexpected type {self.header_type}")


def export_vcf_info_cell(
        header_id: Optional[str] = None,
        number: VCFHeaderNumber = VCFHeaderNumberSpecial.UNBOUND,
        header_type: VCFHeaderType = VCFHeaderType.String,
        description: Optional[str] = None,
        categories: dict[Any, Any] = None):
    """
    Extend ExportRow and annotate methods with export_column.
    The order of defined methods determines the order that the results will appear in an export file
    :param header_id: ID as it's known in the VCF
    :param number: As in how many repeats, None for "." indicating unbound
    :param header_type: VCF data type
    :param description: Description as it appears in the header
    :param categories: If this export column is only valid in some contexts, provide it here
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            return method(*args, **kwargs)

        nonlocal header_id
        nonlocal number
        nonlocal description
        nonlocal categories

        if not header_id:
            header_id = wrapper.__name__
        if description and "$site_name" in description:
            description = description.replace("$site_name", settings.SITE_NAME)

        vcf_header = VCFHeader(
            header_class=VCFHeaderClassStandard.INFO,
            header_id=header_id,
            number=number,
            header_type=header_type,
            description=description
        )

        # have to cache the line number of the source method, otherwise we just get the line number of this wrapper
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.__name__ = method.__name__
        wrapper.is_vcf_export = True
        wrapper.vcf_header = vcf_header
        wrapper.categories = categories
        return wrapper
    return decorator


class VCFCell(Protocol):
    """
    A method annotated with @export_vcf_info_cell will wrap them method into conforming to VCFCell
    """
    vcf_header: VCFHeader
    categories: Optional[dict[str, Any]]
    line_number: int
    def __call__(self, *args) -> Any: ...


class ExportVCF:
    """
    Extend this class and have methods annotated with @vcf_info_cell to be able to make a VCF
    """

    @classmethod
    def get_export_methods(cls, export_tweak: VCFExportTweak = ExportTweak.DEFAULT) -> list[VCFCell]:
        return get_decorated_methods(cls, categories=export_tweak.categories, attribute="is_vcf_export")

    @classmethod
    def complete_header(
            cls,
            genome_build: GenomeBuild,
            contigs: Iterable[Contig],
            extras: list[VCFHeader] = None,
            export_tweak: VCFExportTweak = ExportTweak.DEFAULT) -> list[VCFHeader | str]:
        if not extras:
            extras = []
        return ([
            VCFHeader("fileformat", value="VCFv4.1"),
            VCFHeader("fileDate", value=local_date_str_no_dash()),
            VCFHeader("source", value=settings.SITE_NAME)
        ] +
            extras +
            [method.vcf_header for method in cls.get_export_methods(export_tweak=export_tweak)] +
            # [VCFHeader(header_class="ALT", header_id="NON_REF", description="Represents any possible alternative allele at this location")] + \
            ExportVCF.contig_headers(genome_build, contigs) +
            [("#" + "\t".join(["CHROM", "ID", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO"]))])  # NO SUPPORT FOR FORMAT, etc

    @staticmethod
    def contig_headers(genome_build: GenomeBuild, contigs: Iterable[Contig]) -> list[VCFHeader]:
        """
        Creates a list of contig headers
        """
        out: list[VCFHeader] = []
        assembly = genome_build.name
        for gbc in GenomeBuildContig.objects.filter(genome_build=genome_build, contig__in=contigs).order_by('order').select_related('contig'):
            contig = gbc.contig
            out.append(
                VCFHeader(
                    "contig",
                    header_id=contig.name,
                    custom_attributes={"length": contig.length, "assembly": assembly}
                )
            )
        return out

    def get_variant(self) -> Variant:
        raise NotImplementedError("ExportVCF needs to implement get_variant")

    def get_allele(self) -> Allele:
        # defaults to getting the allele from the variant, but can be overridden to make more efficient
        return self.get_variant().allele

    # noinspection PyMethodMayBeStatic
    def get_qual(self) -> str:
        return "."

    # noinspection PyMethodMayBeStatic
    def get_filter(self) -> str:
        return "PASS"

    def get_variant_id(self) -> str:
        # for the ID column, defaults to allele ID
        return f"a{self.get_allele().pk}"

    def vcf_row(self, export_tweak: ExportTweak = ExportTweak.DEFAULT) -> Optional[str]:
        if (variant := self.get_variant()) and (locus := variant.locus) and (contig := locus.contig):
            # standard columns
            chrom = contig.name
            pos = locus.position
            ref = locus.ref.seq
            alt = variant.alt.seq
            qual = self.get_qual()
            filter_val = self.get_filter()

            # info column
            info = []
            for method in self.__class__.get_export_methods(export_tweak=export_tweak):
                result = method(self)
                if cell := method.vcf_header.format_info_value(result):
                    info.append(cell)

            info_str = ";".join(str(info_cell) for info_cell in info)
            return "\t".join([chrom, self.get_variant_id(), str(pos), ref, alt, qual, filter_val, info_str])
