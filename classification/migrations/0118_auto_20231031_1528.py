# Generated by Django 4.2.5 on 2023-10-30 23:01
from typing import Any

from django.db import migrations

_horak_data = [
    {
        "code": "OVS1",
        "description": "Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a bona fide tumour suppressor gene.",
        "evidence_category": "Predictive Data",
        "score": "8",
        "sublabel": "Null variant in tumour suppressor"
    },
    {
        "code": "OS1",
        "description": "Same amino acid change as a previously established oncogenic variant (using this standard) regardless of nucleotide change. Example: Val→Leu caused by either G>C or G>T in the same codon.",
        "evidence_category": "Predictive Data",
        "score": "4",
        "sublabel": "Same amino acid change as a previously established oncogenic variant"
    },
    {
        "code": "OS2",
        "description": "Well-established in vitro or in vivo functional studies, supportive of an oncogenic effect of the variant.",
        "evidence_category": "Functional Data",
        "score": "4",
        "sublabel": "Well-established functional studies supportive of an oncogenic effect"
    },
    {
        "code": "OS3",
        "description": "Located in one of the hotspots in cancerhotspots.org with at least 50 samples with a somatic variant at the same amino acid position, and the same amino acid change count in cancerhotspots.org in at least 10 samples.",
        "evidence_category": "Cancer Hotspots",
        "score": "4",
        "sublabel": "Cancer hotspot with high frequency of recurrence"
    },
    {
        "code": "OM1",
        "description": "Located in a critical and well-established part of a functional domain (e.g., active site of an enzyme).",
        "evidence_category": "Predictive Data",
        "score": "2",
        "sublabel": "Located in a critical and well-established part of a functional domain"
    },
    {
        "code": "OM2",
        "description": "Protein length changes as a result of in-frame deletions/insertions in a known oncogene, or tumour suppressor gene or stop-loss variants in a known tumour suppressor gene.",
        "evidence_category": "Predictive Data",
        "score": "2",
        "sublabel": "Protein length changes as a result of in-frame deletions/insertions "
    },
    {
        "code": "OM3",
        "description": "Missense variant at an amino acid residue where a different missense variant determined to be oncogenic (using this standard) has been documented. Amino acid difference from reference amino acid should be greater or at least approximately the same as for missense change determined to be oncogenic.",
        "evidence_category": "Predictive Data",
        "score": "2",
        "sublabel": "Missense change at an amino acid residue where a different missense change determined to be oncogenic has been documented"
    },
    {
        "code": "OM4",
        "description": "Located in one of the hotspots in cancerhotspots.org with <50 samples with a somatic variant at the same amino acid position, and the same amino acid change count in cancerhotspots.org is at least 10.",
        "evidence_category": "Cancer Hotspots",
        "score": "2",
        "sublabel": "Cancer hotspot with moderate frequency of recurrence"
    },
    {
        "code": "OP1",
        "description": "All used lines of computational evidence support an oncogenic effect of a variant (conservation/evolutionary, splicing impact, etc.).",
        "evidence_category": "Computational Evidence",
        "score": "1",
        "sublabel": "All utilised lines of computational evidence support oncogenicity"
    },
    {
        "code": "OP2",
        "description": "Somatic variant in a gene in a malignancy with a single genetic etiology. Example: retinoblastoma is caused by bi-allelic RB1 inactivation.",
        "evidence_category": "Predictive Data",
        "score": "1",
        "sublabel": "Somatic variant in a gene in a malignancy with a single genetic etiology"
    },
    {
        "code": "OP3",
        "description": "Located in one of the hotspots in cancerhotspots.org and the particular amino acid change count in cancerhotspots.org is below 10.",
        "evidence_category": "Cancer Hotspots",
        "score": "1",
        "sublabel": "Cancer hotspots with low frequency of recurrence"
    },
    {
        "code": "OP4",
        "description": "Absent from controls (or at an extremely low frequency) in Genome Aggregation Database (gnomAD).",
        "evidence_category": "Population Data",
        "score": "1",
        "sublabel": "Absent in population databases"
    },
    {
        "code": "SBVS1",
        "description": "Minor allele frequency is >5% in Genome Aggregation Database (gnomAD) in any of 5 general continental populations: African, East Asian, European (Non-Finnish), Latino, and South Asian.",
        "evidence_category": "Population Data",
        "score": "-8",
        "sublabel": "MAF >5%"
    },
    {
        "code": "SBS1",
        "description": "Minor allele frequency is >1% in Genome Aggregation Database (gnomAD) in any of 5 general continental populations: African, East Asian, European (Non-Finnish), Latino, and South Asian.",
        "evidence_category": "Population Data",
        "score": "-4",
        "sublabel": "MAF >1%"
    },
    {
        "code": "SBS2",
        "description": "Well-established in vitro or in vivo functional studies show no oncogenic effects.",
        "evidence_category": "Functional Data",
        "score": "-4",
        "sublabel": "Well-established functional studies show no oncogenic effects"
    },
    {
        "code": "SBP1",
        "description": "All used lines of computational evidence suggest no effect of a variant (conservation/evolutionary, splicing impact, etc.).",
        "evidence_category": "Computational Evidence",
        "score": "-1",
        "sublabel": "All utilised lines of computational evidence suggest no impact of a variant"
    },
    {
        "code": "SBP2",
        "description": "A synonymous (silent) variant for which splicing prediction algorithms predict no effect on the splice consensus sequence nor the creation of a new splice site and the nucleotide is not highly conserved.",
        "evidence_category": "Predictive Data",
        "score": "-1",
        "sublabel": "Silent mutation (no predicted impact on splicing)"
    }
]

_category_map = {
    "Predictive Data": "CP",
    "Functional Data": "F",
    "Cancer Hotspots": "F",
    "Gene": "H",
    "Computational Evidence": "CP",
    "Population Data": "P"
}

_points_to_criteria = {
    -8: "BA",
    -4: "BS",
    -2: "BM",
    -1: "BP",
     1: "PP",
     2: "PM",
     4: "PS",
     8: "PVS"
}


def _install_horak_evidence_keys(apps, _schema_editor):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')

    for entry_dict in _horak_data:
        key = "horak:" + entry_dict.get("code").lower()
        label = entry_dict.get("code")
        sub_label = entry_dict.get("sublabel")
        description = entry_dict.get("description")
        evidence_category = _category_map.get(entry_dict.get("evidence_category")) or "H"
        value_type = "C"
        default_crit_evaluation = _points_to_criteria.get(int(entry_dict.get("score")))
        order = 10000
        EvidenceKey.objects.create(
            key=key,
            label=label,
            sub_label=sub_label,
            description=description,
            evidence_category=evidence_category,
            value_type=value_type,
            default_crit_evaluation=default_crit_evaluation,
            crit_uses_points=True,
            order=order
        )


def _update_assertion_methods(apps, _schema_editor):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')
    assertion_method = EvidenceKey.objects.get(key="assertion_method")
    options: list[dict[str, Any]] = assertion_method.options
    has_horak = False
    for option in options:
        key = option.get("key")
        if "acmg" in key.lower():
            option["namespaces"] = ["acmg"]
        elif "sherloc" in key:
            option["namespaces"] = ["sherloc"]
        if key == "horak":
            has_horak = True
    if not has_horak:
        options.append({
            "key": "horak",
            "index": 5,
            "label": "Horak (PMID:35101336)",
            "namespaces": ["horak"]
        })
    assertion_method.save()


def _remove_horak_evidence_keys(apps, _schema_editor):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')
    EvidenceKey.objects.filter(key__startswith="horak:").delete()


def _reverse_update_assertion_methods(apps, _schema_editor):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')
    assertion_method = EvidenceKey.objects.get(key="assertion_method")
    options: list[dict[str, Any]] = assertion_method.options
    updated_options = []
    for option in options:
        option.pop("namespaces", None)
        if option.get("key") != "horak":
            updated_options.append(option)
    assertion_method.options = updated_options
    assertion_method.save()


class Migration(migrations.Migration):
    dependencies = [
        ('classification', '0117_evidencekey_crit_allows_override_strengths_and_more'),
    ]

    operations = [
        migrations.RunPython(code=_install_horak_evidence_keys, reverse_code=_remove_horak_evidence_keys),
        migrations.RunPython(code=_update_assertion_methods, reverse_code=_reverse_update_assertion_methods)
    ]
