"""
`vg inspect sample <pk>`: a Sample (one genotype column of a VCF) - its VCF and build, patient, import
status, genotype stats, cohorts it belongs to, the active gene list, attached files, analyses that
use it and template runs launched for it.
"""
from typing import Any

from analysis.models.models_analysis import SampleAnalysisTemplateRun
from analysis.models.nodes.sources.sample_node import SampleNode
from genes.models.models_gene_list import ActiveSampleGeneList
from library.vg.inspect import capped, ref
from snpdb.models import CohortSample, Sample, SampleFilePath


def load(key: str) -> Sample:
    try:
        return Sample.objects.select_related("vcf__genome_build", "vcf__user", "patient").get(pk=int(key))
    except (ValueError, Sample.DoesNotExist) as e:
        raise LookupError(f"No Sample with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    sample = load(key)
    vcf = sample.vcf
    data: dict[str, Any] = {
        "id": sample.pk,
        "name": sample.name,
        "vcf_sample_name": sample.vcf_sample_name,
        "vcf": {**ref("vcf", vcf, vcf.name or ""), "build": vcf.genome_build_id, "user": vcf.user.username, "date": vcf.date.date().isoformat()},
        "import_status": sample.get_import_status_display(),
        "variants_type": sample.get_variants_type_display(),
        "patient": ref("patient", sample.patient, str(sample.patient)) if sample.patient_id else None,
        "no_dna_control": sample.no_dna_control,
        "url": sample.get_absolute_url(),
        "genotype_stats": _genotype_stats(sample),
    }
    active = ActiveSampleGeneList.objects.filter(sample=sample).select_related("sample_gene_list__gene_list").first()
    data["active_gene_list"] = ref("gene_list", active.sample_gene_list.gene_list, active.sample_gene_list.gene_list.name) if active else None
    if depth >= 2:
        data["cohorts"] = capped(CohortSample.objects.filter(sample=sample).select_related("cohort").order_by("cohort_id"), _cohort)
        data["files"] = capped(SampleFilePath.objects.filter(sample=sample).order_by("pk"), lambda f: {"type": f.get_file_type_display(), "path": f.file_path})
        data["analyses_using"] = capped(SampleNode.objects.filter(sample=sample).select_related("analysis").order_by("analysis_id"),
                                        lambda node: {**ref("analysis", node.analysis, node.analysis.name), "node": node.pk})
        data["template_runs"] = capped(SampleAnalysisTemplateRun.objects.filter(sample=sample).select_related("analysis_template_run__analysis").order_by("pk"),
                                       lambda run: ref("analysis", run.analysis_template_run.analysis, run.analysis_template_run.analysis.name))
    return data


def _genotype_stats(sample: Sample) -> dict[str, Any] | None:
    stats = sample.get_genotype_stats()
    if stats is None:
        return None
    return {"variants": stats.variant_count, "het": stats.het_count, "hom": stats.hom_count, "ref": stats.ref_count,
            "unknown": stats.unk_count, "snps": stats.snp_count, "insertions": stats.insertions_count, "deletions": stats.deletions_count}


def _cohort(cohort_sample: CohortSample) -> dict[str, Any]:
    cohort = cohort_sample.cohort
    return {**ref("cohort", cohort, cohort.name), "vcf_cohort": cohort.vcf_id is not None, "samples": cohort.sample_count,
            "packed_index": cohort_sample.cohort_genotype_packed_field_index}
