"""
The icon for each UploadedFileTypes value, as shown on the upload page and the upload pipeline page.

Entry point is file_type_icon_html, behind the {% file_type_icon %} tag (upload/templatetags/upload_tags.py) and the
upload poll JSON, so the server-rendered table and the rows the page's JS adds agree. Where the site already has an
icon for the concept the file wears the same one - the analysis node badges (BED / IntersectionNode, pedigree /
PedigreeNode, tags / TagNode, classifications / ClassificationsNode) and the preview icons on search results (patient
records / Cohort, ClinVar / the server status card, analysis / Analysis) - and the VCF sub-types are the VCF with a
badge in front. The drawn symbols (file-icon-*) are in uicore/templates/uicore/tags/svg_icon_sprite.html; sizing and
the badge are .file-type-icon / .file-type-badge in global.scss.
"""
from dataclasses import dataclass
from typing import Optional

from django.utils.html import format_html
from django.utils.safestring import SafeString

from library.preview_request import fa_icon_html, svg_symbol_icon_html
from upload.models.models_enums import UploadedFileTypes


@dataclass(frozen=True)
class FileTypeIcon:
    """ Exactly one of fa (FontAwesome classes) or symbol (an id in svg_icon_sprite.html) """
    fa: Optional[str] = None
    symbol: Optional[str] = None
    css: str = ""
    badge: Optional["FileTypeIcon"] = None

    def glyph_html(self) -> SafeString:
        if self.symbol:
            return svg_symbol_icon_html(self.symbol)
        return fa_icon_html(self.fa)


_VCF = FileTypeIcon(symbol="file-icon-vcf")
_GENE_LIST_G = FileTypeIcon(symbol="node-icon-gene-list", css="file-type-badge-gene-list")

FILE_TYPE_ICONS: dict[UploadedFileTypes, FileTypeIcon] = {
    UploadedFileTypes.ANALYSIS: FileTypeIcon(fa="fa-solid fa-diagram-project"),
    UploadedFileTypes.BED: FileTypeIcon(symbol="node-icon-intervals", css="file-type-icon-bed"),
    UploadedFileTypes.CLINVAR: FileTypeIcon(fa="fa-solid fa-earth-americas"),
    UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS: FileTypeIcon(symbol="file-icon-fusion"),
    UploadedFileTypes.GENE_LIST: FileTypeIcon(symbol="file-icon-list", css="file-type-icon-muted", badge=_GENE_LIST_G),
    UploadedFileTypes.GENE_COVERAGE: FileTypeIcon(symbol="file-icon-gene-coverage"),
    UploadedFileTypes.LIFTOVER: FileTypeIcon(symbol="file-icon-liftover"),
    UploadedFileTypes.MANUAL_VARIANT_ENTRY: FileTypeIcon(symbol=_VCF.symbol,
                                                         badge=FileTypeIcon(fa="fa-solid fa-pen", css="file-type-badge-manual")),
    UploadedFileTypes.PED: FileTypeIcon(symbol="node-icon-pedigree"),
    UploadedFileTypes.PATIENT_RECORDS: FileTypeIcon(fa="fa-solid fa-users"),
    UploadedFileTypes.VARIANT_CLASSIFICATIONS: FileTypeIcon(fa="fa-solid fa-clipboard-check"),
    UploadedFileTypes.VARIANT_TAGS: FileTypeIcon(fa="fa-solid fa-tags"),
    UploadedFileTypes.VCF: _VCF,
    UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY: FileTypeIcon(symbol=_VCF.symbol,
                                                             badge=FileTypeIcon(fa="fa-solid fa-v", css="file-type-badge-variant")),
    UploadedFileTypes.GENE_LEVEL_INSERT_VARIANTS_ONLY: FileTypeIcon(symbol=_VCF.symbol,
                                                                    badge=FileTypeIcon(fa="fa-solid fa-dna", css="file-type-badge-gene")),
    UploadedFileTypes.GENE_LEVEL_CNV_VCF: FileTypeIcon(symbol="file-icon-copy-number"),
    UploadedFileTypes.WIKI_GENE: FileTypeIcon(fa="fa-brands fa-wikipedia-w"),
    UploadedFileTypes.WIKI_VARIANT: FileTypeIcon(fa="fa-brands fa-wikipedia-w"),
}


def file_type_icon_html(file_type: Optional[str]) -> SafeString:
    """ The 40px icon block for a file type code - an empty block for a file whose type is not known yet """
    if not file_type:
        return SafeString('<div class="file-type-icon"></div>')
    file_type = UploadedFileTypes(file_type)
    icon = FILE_TYPE_ICONS[file_type]
    badge_html = ""
    if icon.badge:
        badge_html = format_html('<span class="file-type-badge {}">{}</span>', icon.badge.css, icon.badge.glyph_html())
    return format_html('<div class="file-type-icon {}" title="{}">{}{}</div>',
                       icon.css, file_type.label, icon.glyph_html(), badge_html)
