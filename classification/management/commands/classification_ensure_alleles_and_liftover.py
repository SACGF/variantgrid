from django.conf import settings
from django.core.management.base import BaseCommand

from classification.autopopulate_evidence_keys.evidence_from_variant import \
    get_clingen_allele_and_evidence_value_for_variant
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from library.git import Git
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import Variant, Allele
from snpdb.models.models_enums import ImportSource
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--add-clingen-allele', action='store_true',
                            help="Add 'clingen_allele_id' to Classifications missing it")

    def handle(self, *args, **options):
        script = __file__
        add_clingen_allele = options["add_clingen_allele"]
        for genome_build in GenomeBuild.builds_with_annotation():
            defaults = {"git_hash": Git(settings.BASE_DIR).hash}

            # Note: This deliberately only gets classifications where the submitting variant was against this genome build
            # ie we don't use Classification.get_variant_q_from_classification_qs() to get liftovers
            contigs_q = Variant.get_contigs_q(self.genome_build)
            variants_qs = Variant.objects.filter(contigs_q, importedalleleinfo__isnull=False)
            if variants_qs.exists():
                print(f"{genome_build} has variants - creating Allele/ClinGen + liftover")
                populate_clingen_alleles_for_variants(genome_build, variants_qs)

                allele_qs = Allele.objects.filter(variantallele__variant__in=variants_qs)
                create_liftover_pipelines(admin_bot(), allele_qs, ImportSource.COMMAND_LINE, genome_build)

                if add_clingen_allele:
                    # Patch those ClinGen alleles into the variant classifications
                    num_added_clingen_allele = 0
                    clingen_allele_key_null = "evidence__%s__isnull" % SpecialEKeys.CLINGEN_ALLELE_ID
                    for vc in Classification.objects.filter(variant__in=variants_qs,
                                                            **{clingen_allele_key_null: True}):
                        _, evidence_value, _ = get_clingen_allele_and_evidence_value_for_variant(genome_build,
                                                                                                 vc.variant)
                        vc.patch_value({SpecialEKeys.CLINGEN_ALLELE_ID: evidence_value},
                                       source=SubmissionSource.VARIANT_GRID)
                        vc.save()
                        num_added_clingen_allele += 1

                    print(f"Added {SpecialEKeys.CLINGEN_ALLELE_ID} to {num_added_clingen_allele} classifications")
