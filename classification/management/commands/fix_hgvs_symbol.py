from typing import Optional

from django.core.management import BaseCommand

from classification.enums import SpecialEKeys
from classification.models import Classification, ClassificationModification
from genes.hgvs import CHGVS


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--test', action='store_true', default=False)

    def handle(self, *args, **options):

        test = options["test"]
        for c_id, c_hgvs, gene_symbol in Classification.objects.filter(evidence__c_hgvs__value__isnull=False).values_list("id", "evidence__c_hgvs__value", "evidence__gene_symbol__value"):
            c_parts = CHGVS(c_hgvs)
            if not c_parts.gene and gene_symbol and gene_symbol not in c_hgvs and ":c." in c_hgvs:
                desired_chgvs = c_parts.with_gene_symbol(gene_symbol).full_c_hgvs

                c: Classification
                c = Classification.objects.get(pk=c_id)
                evidence = c.evidence
                evidence[SpecialEKeys.C_HGVS]["value"] = desired_chgvs
                if not test:
                    c.save(update_fields=["evidence"])

                changed_vcs = 0
                vc: ClassificationModification
                for vc in ClassificationModification.objects.filter(classification_id=c_id):
                    evidence_changed = False
                    if vc_evidence := vc.published_evidence:
                        if SpecialEKeys.C_HGVS in vc_evidence:
                            if vc_evidence[SpecialEKeys.C_HGVS].get("value") == c_hgvs:
                                vc_evidence[SpecialEKeys.C_HGVS]["value"] = desired_chgvs
                                evidence_changed = True

                    delta_changed = False
                    if vc_delta := vc.delta:
                        if SpecialEKeys.C_HGVS in vc_delta:
                            if vc_delta[SpecialEKeys.C_HGVS].get("value") == c_hgvs:
                                vc_delta[SpecialEKeys.C_HGVS]["value"] = desired_chgvs
                                delta_changed = True

                    if not test:
                        vc.save(update_fields=["published_evidence", "delta"])

                print(f"Classification classification {c_id} {c_hgvs} - > {desired_chgvs}")
