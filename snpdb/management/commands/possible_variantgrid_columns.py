from django.core.management.base import BaseCommand

from annotation.models import ClinVar, VariantAnnotation, GeneAnnotation, DBNSFPGeneAnnotation
from genes.models import Gene, TranscriptVersion, GeneVersion, GeneSymbolWiki, HGNC, UniProt
from library.django_utils import get_model_fields
from snpdb.models import Variant, Locus, VariantZygosityCount, Contig, ClinGenAllele, VariantWiki, VariantGridColumn


class Command(BaseCommand):
    """
        This looks to see whether we forgot to add any eg Annotation fields as VariantGrid columns
    """
    def handle(self, *args, **options):
        model_ignore_fields = {"id", "variant", "created", "modified", "version", "genome_build",
                               "import_source", "hgnc_import", "cached_web_resource"}

        models_by_path = {
            "": Variant,
            'clinvar': ClinVar,
            'global_variant_zygosity': VariantZygosityCount,
            'locus': Locus,
            'locus__contig': Contig,
            'variantallele__allele__clingen_allele': ClinGenAllele,
            'variantannotation': VariantAnnotation,
            'variantannotation__gene': Gene,
            'variantannotation__gene__geneannotation': GeneAnnotation,
            'variantannotation__gene__geneannotation__dbnsfp_gene': DBNSFPGeneAnnotation,
            'variantannotation__transcript_version': TranscriptVersion,
            'variantannotation__transcript_version__gene_version': GeneVersion,
            'variantannotation__transcript_version__gene_version__gene_symbol__genesymbolwiki': GeneSymbolWiki,
            'variantannotation__transcript_version__gene_version__hgnc': HGNC,
            'variantannotation__transcript_version__gene_version__hgnc__uniprot': UniProt,
            'variantwiki': VariantWiki,
        }

        # Manually reviewed - decided we don't want these
        ignored_columns = {
            "alt",
            "global_variant_zygosity__collection",
            "locus",
            "locus__contig",
            "locus__contig__assigned_molecule",
            "locus__contig__genbank_accession",
            "locus__contig__length",
            "locus__contig__molecule_type",
            "locus__contig__role",
            "locus__contig__ucsc_name",
            "locus__ref",
            "variantallele__allele__clingen_allele__api_response",
            "variantannotation__annotation_run",
            "variantannotation__gene",
            "variantannotation__gene__annotation_consortium",
            "variantannotation__gene__geneannotation__dbnsfp_gene",
            "variantannotation__gene__geneannotation__dbnsfp_gene__ensembl_transcript",
            "variantannotation__gene__geneannotation__dbnsfp_gene__gene_symbol",
            "variantannotation__gene__geneannotation__dbnsfp_gene__refseq_transcript",
            "variantannotation__gene__geneannotation__gene",
            "variantannotation__transcript",
            "variantannotation__transcript_version",
            "variantannotation__transcript_version__contig",
            "variantannotation__transcript_version__data",
            "variantannotation__transcript_version__gene_version",
            "variantannotation__transcript_version__gene_version__gene",
            "variantannotation__transcript_version__gene_version__gene_symbol",
            "variantannotation__transcript_version__gene_version__gene_symbol__genesymbolwiki__gene_symbol",
            "variantannotation__transcript_version__gene_version__gene_symbol__genesymbolwiki__last_edited_by",
            "variantannotation__transcript_version__gene_version__gene_symbol__genesymbolwiki__wiki_ptr",
            "variantannotation__transcript_version__gene_version__hgnc",
            "variantannotation__transcript_version__gene_version__hgnc__ensembl_gene_id",
            "variantannotation__transcript_version__gene_version__hgnc__gene_symbol",
            "variantannotation__transcript_version__gene_version__hgnc__refseq_ids",
            "variantannotation__transcript_version__gene_version__hgnc__status",
            "variantannotation__transcript_version__gene_version__hgnc__uniprot",
            "variantannotation__transcript_version__gene_version__hgnc__uniprot_ids",
            "variantannotation__transcript_version__gene_version__hgnc_identifier",
            "variantannotation__transcript_version__transcript",
            "variantannotation__vep_skipped_reason",
            "variantwiki__last_edited_by",
            "variantwiki__wiki_ptr",
        }

        possible_paths = set()
        for path_prefix, model in models_by_path.items():
            fields = get_model_fields(model, ignore_fields=model_ignore_fields)
            if path_prefix:
                prefix = path_prefix + "__"
            else:
                prefix = ""

            possible_paths.update([f"{prefix}{f}" for f in fields])

        vgc_qs = VariantGridColumn.objects.filter(queryset_field=True)
        existing_variantgrid_columns = set(vgc_qs.values_list("variant_column", flat=True))

        unused = possible_paths - ignored_columns - existing_variantgrid_columns
        if unused:
            print(f"{len(unused)} unused possible paths:")
            for path in unused:
                print(path)
        else:
            print(f"All possible paths ({len(possible_paths)}) already assigned VariantGridColumns")
