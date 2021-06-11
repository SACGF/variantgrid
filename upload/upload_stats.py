from collections import Counter, defaultdict
from django.db.models import ExpressionWrapper, F, fields
from django.db.models.aggregates import Sum, Count
from django.utils import timezone
from django.utils.datastructures import OrderedSet

import numpy as np
import pandas as pd
from snpdb.models import VCF
from upload.models import UploadStep


def get_upload_stats(uploadstep_qs, max_step_names: int = None):
    """ returns (num_upload_pipelines, total_times, time_per_kilo_variant) """
    step_name_durations, num_upload_pipelines = get_step_durations(uploadstep_qs, max_step_names=max_step_names)
    total_times, time_per_kilo_variant = get_upload_stats_from_durations(step_name_durations)
    return num_upload_pipelines, total_times, time_per_kilo_variant


def get_step_durations(uploadstep_qs, max_step_names: int = None):
    duration = ExpressionWrapper(F('end_date') - F('start_date'), output_field=fields.DurationField())
    ups_qs = uploadstep_qs.annotate(duration=duration).filter(duration__isnull=False)

    step_name_durations = defaultdict(lambda: ([], []))
    upload_pipelines = set()
    values_qs = ups_qs.values_list("upload_pipeline", "name", "items_processed", "duration")
    for upload_pipeline, name, items_processed, duration in values_qs:
        step_name_durations[name][0].append(items_processed)
        step_name_durations[name][1].append(duration.total_seconds())
        upload_pipelines.add(upload_pipeline)
    num_upload_pipelines = len(upload_pipelines)

    # Turn into numpy list
    step_name_durations = {k: (np.array(v[0]), np.array(v[1])) for k, v in step_name_durations.items()}

    # Collapse smaller steps into "other" bucket
    if max_step_names and len(step_name_durations) > max_step_names:
        longest_steps = (s[0] for s in reversed(sorted(step_name_durations.items(),
                                                       key=lambda i: i[1][1].sum())))
        collapse_steps = list(longest_steps)[max_step_names:]
        other_items_processed = []
        other_durations = []
        for s in collapse_steps:
            v = step_name_durations.pop(s)
            other_items_processed.append(v[0])
            other_durations.append(v[1])

        items_processed_array = np.concatenate(other_items_processed)
        durations_array = np.concatenate(other_durations)
        collapsed_steps_name = f"Other ({len(collapse_steps)} steps)"
        step_name_durations[collapsed_steps_name] = (items_processed_array, durations_array)

    return step_name_durations, num_upload_pipelines


def get_upload_stats_from_durations(step_name_durations):
    """ returns (total_times, time_per_kilo_variant) """

    total_times = {}
    time_per_kilo_variant = {}

    for name, (items_processed_array, durations_array) in step_name_durations.items():
        total_times[name] = durations_array.sum()

        non_zero = np.nonzero(items_processed_array)
        items_processed_array = items_processed_array[non_zero]
        if items_processed_array.any():
            durations_array = durations_array[non_zero]
            tpkv = durations_array / items_processed_array
            time_per_kilo_variant[name] = tpkv.tolist()

    return total_times, time_per_kilo_variant


def get_vcf_variant_upload_stats():
    """ df index of vcf_ids, cols = [cumulative_samples, total_variants, percent_known] """

    def get_totals_per_vcf(ups_name):
        qs = VCF.objects.filter(uploadedvcf__upload_pipeline__uploadstep__name=ups_name)
        qs = qs.annotate(total_items_processed=Sum("uploadedvcf__upload_pipeline__uploadstep__items_processed"))
        totals = {}
        for vcf_id, total_items_processed in qs.order_by("pk").values_list("pk", "total_items_processed"):
            totals[vcf_id] = total_items_processed

        return totals

    samples_per_vcf = {}
    for pk, num_samples in VCF.objects.annotate(num_samples=Count("sample")).values_list("pk", "num_samples"):
        samples_per_vcf[pk] = num_samples

    unknown = get_totals_per_vcf(UploadStep.CREATE_UNKNOWN_LOCI_AND_VARIANTS_TASK_NAME)
    total = get_totals_per_vcf(UploadStep.PROCESS_VCF_TASK_NAME)
    old_genotypes = get_totals_per_vcf("ObservedVariants SQL COPY")  # from old pipeline, obsolete as of May 2018
    new_genotypes = get_totals_per_vcf("CohortGenotypeCollection SQL COPY")

    data = {"num_samples": samples_per_vcf,
            "total": total,
            "unknown": unknown,
            "old_genotypes": old_genotypes,  # This copied each genotype individually
            "new_genotypes": new_genotypes}  # This has 1 genotype per variant - needs to be multiplied by samples

    df = pd.DataFrame(data=data)
    df = df.replace(np.NaN, 0)

    df["total_genotypes"] = df["old_genotypes"] + df["new_genotypes"] * df["num_samples"]
    df["unknown"] = df["unknown"] / 2  # We create ref and alt for unknowns
    df["cumulative_samples"] = np.cumsum(df["num_samples"])
    df["cumulative_variants"] = np.cumsum(df["total"])
    df["cumulative_genotypes"] = np.cumsum(df["total_genotypes"])
    known = df["total"] - df["unknown"]
    df["percent_known"] = 100.0 * known / df["total"]

    return df


def get_step_total_stats(upload_pipeline):
    """ Returns a list of tuples of (name, count) in start_date order """

    step_order = OrderedSet()
    total_seconds = Counter()
    total_items = Counter()
    total_steps = Counter()

    ups_qs = upload_pipeline.uploadstep_set.all()
    ups_qs = ups_qs.filter(parent_upload_step__isnull=True).order_by("start_date")  # skip sub-steps
    ups_values = ups_qs.values_list("name", "items_processed", "start_date", "end_date")
    if ups_values.exists():
        for name, items_processed, start_date, end_date in ups_values:
            step_order.add(name)
            total_steps[name] += 1
            total_items[name] += items_processed or 0
            if start_date and end_date:
                total_seconds[name] += (end_date - start_date).total_seconds()

    step_total_stats_list = []
    for step_name in step_order:
        step_total_stats_list.append((step_name, total_steps[step_name], total_items[step_name], total_seconds[step_name]))

    return step_total_stats_list


def get_step_order_and_step_start_end_lines(upload_pipeline):
    """  step_start_end_lines = dict of lists of lines (which are list of start/end tuples) """

    ups_qs = upload_pipeline.uploadstep_set.all()
    now = timezone.now()

    ups_qs = ups_qs.filter(upload_pipeline=upload_pipeline,
                           name__isnull=False,
                           parent_upload_step__isnull=True,  # Skip sub-steps
                           start_date__isnull=False)
    ups_qs = ups_qs.order_by("start_date")

    step_order = OrderedSet()
    step_start_end_lines = defaultdict(list)
    pipeline_start = None

    for name, status, start_date, end_date in ups_qs.values_list("name", "status", "start_date", "end_date"):
        if end_date is None:
            end_date = now  # Still running...

        step_order.add(name)
        if pipeline_start is None:
            pipeline_start = start_date

        start_secs = (start_date - pipeline_start).total_seconds()
        end_secs = (end_date - pipeline_start).total_seconds()

        ssel = step_start_end_lines[name]
        found = False
        if ssel:
            # Put it in 1st spot it fits
            for line in ssel:
                last = line[-1]
                if start_secs > last[1]:
                    line.append((start_secs, end_secs, status))
                    found = True
                    break

        if not found:
            ssel.append([(start_secs, end_secs, status)])

    return step_order, step_start_end_lines
