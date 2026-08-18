"""
Tag stats page - @see https://github.com/SACGF/variantgrid/issues/1751

The page is a skeleton of cards, each lazy loading its own JSON endpoint so first paint is instant. Every
endpoint caches its response in Redis keyed by a hash of its parameters, and reports when it was calculated
so users can see how stale a card is.

Counts are labelled by what they count. Tagging repeats a lot (a known artefact gets re-tagged in every
analysis it turns up in) so "tag events", "distinct (variant, tag)" and "distinct variants" are all very
different numbers.
"""
import hashlib
import json
from collections import Counter, defaultdict
from datetime import timedelta

from django.conf import settings
from django.contrib.auth.models import User
from django.core.cache import cache
from django.db.models import Count, QuerySet
from django.db.models.functions import TruncMonth, TruncQuarter
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, render
from django.utils.timezone import localtime, now

from analysis.models import VariantTag
from annotation.models import VariantAnnotationVersion
from library.constants import DAY_SECS
from snpdb.models import GenomeBuild, Lab, Tag, UserSettings, Variant
from variantopedia.forms import TagStatsGenesForm, TagStatsTagForm, TagStatsTagsForm

TAG_STATS_CACHE_TIMEOUT = DAY_SECS
DEFAULT_TOP_TAGS = 8
DEFAULT_TOP_GENES = 10
MEGA_ARTEFACT_VARIANTS = 40
OTHER = "other"
GENE_SYMBOL_FIELD = "variant__variantannotation__transcript_version__gene_version__gene_symbol_id"


def _cache_key(card: str, genome_build: GenomeBuild, params: dict) -> str:
    params_hash = hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()[:16]
    return f"tag_stats/{settings.CACHE_VERSION}/{card}/{genome_build.name}/{params_hash}"


def _cached_json(card: str, genome_build: GenomeBuild, params: dict, calculate) -> JsonResponse:
    cache_key = _cache_key(card, genome_build, params)
    data = cache.get(cache_key)
    if data is None:
        data = calculate()
        data["calculated"] = localtime(now()).isoformat()
        cache.set(cache_key, data, timeout=TAG_STATS_CACHE_TIMEOUT)
    return JsonResponse(data)


def _tags_qs(genome_build: GenomeBuild) -> QuerySet[VariantTag]:
    return VariantTag.get_for_build(genome_build)


def _top_names(counter: Counter, top_n: int) -> list[str]:
    return [name for name, _ in counter.most_common(top_n)]


def _grouped_series(counts: dict[tuple, int], buckets: list, top_n: int) -> list[dict]:
    """ counts keyed by (bucket, name) -> a plotly series per top name, everything else summed into "other" """
    totals = Counter()
    for (_bucket, name), count in counts.items():
        totals[name] += count
    top = _top_names(totals, top_n)

    series_counts = defaultdict(Counter)
    for (bucket, name), count in counts.items():
        series_counts[name if name in top else OTHER][bucket] += count

    names = [n for n in top if n in series_counts]
    if OTHER in series_counts:
        names.append(OTHER)
    return [{"name": name, "counts": [series_counts[name].get(b, 0) for b in buckets]} for name in names]


def tag_stats(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    artefact_tag = _default_artefact_tag(genome_build)
    context = {
        "genome_build": genome_build,
        "genes_form": TagStatsGenesForm(prefix="genes"),
        "co_occurrence_form": TagStatsTagsForm(prefix="co-occurrence"),
        "re_tagged_form": TagStatsTagForm(prefix="re-tagged", initial={"tag": artefact_tag}),
        "gene_time_form": TagStatsTagForm(prefix="gene-time", initial={"tag": _most_used_tag(genome_build)}),
    }
    return render(request, 'variantopedia/tag_stats.html', context)


def _default_artefact_tag(genome_build: GenomeBuild) -> str:
    """ Artefact tagging is what drives the re-tagging numbers - fall back to the most used tag """
    if artefact := Tag.objects.filter(pk__iexact="artefact").first():
        return artefact.pk
    return _most_used_tag(genome_build)


def _most_used_tag(genome_build: GenomeBuild) -> str:
    qs = _tags_qs(genome_build).values("tag_id").annotate(count=Count("pk")).order_by("-count")
    if top := qs.first():
        return top["tag_id"]
    return ""


def tag_stats_headline(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    def calculate():
        qs = _tags_qs(genome_build)
        year_ago = now() - timedelta(days=365)
        return {
            "tag_events": qs.count(),
            "distinct_variant_tags": qs.values("variant_id", "tag_id").distinct().count(),
            "distinct_variants": qs.values("variant_id").distinct().count(),
            "analyses": qs.exclude(analysis__isnull=True).values("analysis_id").distinct().count(),
            "taggers_this_year": qs.filter(created__gte=year_ago).values("user_id").distinct().count(),
        }

    return _cached_json("headline", genome_build, {}, calculate)


def tag_stats_over_time(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    top_n = _get_int(request, "top_n", DEFAULT_TOP_TAGS)

    def calculate():
        qs = _tags_qs(genome_build).annotate(month=TruncMonth("created"))
        counts = {}
        totals = Counter()
        months = set()
        for row in qs.values("month", "tag_id").annotate(count=Count("pk")):
            month = row["month"].strftime("%Y-%m")
            months.add(month)
            counts[(month, row["tag_id"])] = row["count"]
            totals[row["tag_id"]] += row["count"]

        sorted_months = sorted(months)
        return {
            "months": sorted_months,
            "series": _grouped_series(counts, sorted_months, top_n),
            "totals": dict(totals.most_common()),
        }

    return _cached_json("over_time", genome_build, {"top_n": top_n}, calculate)


def tag_stats_for_user(request, genome_build_name, user_id=None):
    """ Shared by the stats page's "Your tagging" card and the user page section """
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    user = get_object_or_404(User, pk=user_id) if user_id else request.user

    def calculate():
        qs = _tags_qs(genome_build).filter(user=user)
        tag_counts = qs.values("tag_id").annotate(count=Count("pk")).order_by("-count")

        thirty_days_ago = (now() - timedelta(days=30)).date()
        daily = Counter(created.date() for created in
                        qs.filter(created__date__gte=thirty_days_ago).values_list("created", flat=True))
        days = [thirty_days_ago + timedelta(days=i) for i in range(31)]

        return {
            "username": str(user),
            "tag_events": qs.count(),
            "distinct_variants": qs.values("variant_id").distinct().count(),
            "top_tags": [{"tag": r["tag_id"], "count": r["count"]} for r in tag_counts[:DEFAULT_TOP_TAGS]],
            "recent_days": [d.isoformat() for d in days],
            "recent_counts": [daily.get(d, 0) for d in days],
        }

    return _cached_json("user", genome_build, {"user_id": user.pk}, calculate)


def tag_stats_by_lab(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    def calculate():
        qs = _tags_qs(genome_build)
        user_counts = qs.values("user_id", "user__username").annotate(count=Count("pk")).order_by("-count")

        labs_by_group_name = {lab.group_name: lab for lab in Lab.objects.exclude(group_name__isnull=True)}
        lab_names_by_user_id = defaultdict(list)
        for user_id, group_name in User.objects.filter(groups__name__in=labs_by_group_name).values_list(
                "pk", "groups__name"):
            if lab := labs_by_group_name.get(group_name):
                lab_names_by_user_id[user_id].append(str(lab))

        users = []
        lab_counts = Counter()
        for row in user_counts:
            lab_names = lab_names_by_user_id.get(row["user_id"], [])
            users.append({"user": row["user__username"], "labs": lab_names, "count": row["count"]})
            for lab_name in lab_names:
                lab_counts[lab_name] += row["count"]

        return {
            "users": users,
            "labs": [{"lab": lab_name, "count": count} for lab_name, count in lab_counts.most_common()],
        }

    return _cached_json("by_lab", genome_build, {}, calculate)


def _with_gene_symbol(genome_build: GenomeBuild, tags_qs: QuerySet[VariantTag]) -> QuerySet[VariantTag]:
    """ Restrict to tags whose variant has a representative transcript in the latest annotation, so
        GENE_SYMBOL_FIELD can be grouped on """
    vav = VariantAnnotationVersion.latest(genome_build)
    return tags_qs.filter(variant__variantannotation__version=vav,
                          variant__variantannotation__transcript_version__isnull=False)


def tag_stats_genes(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    gene_symbols = _get_list(request, "gene_symbols")
    tag_ids = _get_list(request, "tags")
    top_n = _get_int(request, "top_n", DEFAULT_TOP_GENES)

    def calculate():
        tags_qs = _tags_qs(genome_build)
        if tag_ids:
            tags_qs = tags_qs.filter(tag__in=tag_ids)
        qs = _with_gene_symbol(genome_build, tags_qs)

        counts = {}
        gene_totals = Counter()
        for row in qs.values(GENE_SYMBOL_FIELD, "tag_id").annotate(count=Count("pk")):
            gene_symbol = row[GENE_SYMBOL_FIELD]
            if gene_symbols and gene_symbol not in gene_symbols:
                continue
            counts[(gene_symbol, row["tag_id"])] = row["count"]
            gene_totals[gene_symbol] += row["count"]

        genes = gene_symbols or _top_names(gene_totals, top_n)
        genes = [g for g in genes if gene_totals.get(g)]
        return {
            "genes": genes,
            "series": _grouped_series({k: v for k, v in counts.items() if k[0] in genes}, genes, DEFAULT_TOP_TAGS),
        }

    params = {"gene_symbols": sorted(gene_symbols), "tags": sorted(tag_ids), "top_n": top_n}
    return _cached_json("genes", genome_build, params, calculate)


def tag_stats_co_occurrence(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    def calculate():
        qs = _tags_qs(genome_build).filter(allele__isnull=False)
        tags_by_allele = defaultdict(set)
        for allele_id, tag_id in qs.values_list("allele_id", "tag_id").distinct():
            tags_by_allele[allele_id].add(tag_id)

        pair_counts = Counter()
        tag_totals = Counter()
        for tag_ids in tags_by_allele.values():
            sorted_tags = sorted(tag_ids)
            for tag_id in sorted_tags:
                tag_totals[tag_id] += 1
            for i, tag_a in enumerate(sorted_tags):
                for tag_b in sorted_tags[i + 1:]:
                    pair_counts[(tag_a, tag_b)] += 1

        tags = _top_names(tag_totals, DEFAULT_TOP_TAGS * 2)
        tags.sort()
        matrix = [[pair_counts.get(tuple(sorted((tag_a, tag_b))), 0) if tag_a != tag_b else None
                   for tag_b in tags] for tag_a in tags]
        return {
            "tags": tags,
            "matrix": matrix,
            "top_pairs": [{"tags": list(pair), "alleles": count} for pair, count in pair_counts.most_common(20)],
        }

    return _cached_json("co_occurrence", genome_build, {}, calculate)


def tag_stats_re_tagged(request, genome_build_name):
    """ Variants that get tagged over and over - the case for excluding them in analysis templates """
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    tag_id = request.GET.get("tag") or _default_artefact_tag(genome_build)

    def calculate():
        qs = _tags_qs(genome_build).filter(tag=tag_id)
        top_qs = qs.values("variant_id").annotate(count=Count("pk")).order_by("-count")[:MEGA_ARTEFACT_VARIANTS]
        top = list(top_qs)
        variant_qs = Variant.objects.filter(pk__in=[r["variant_id"] for r in top])
        variants_by_id = {v.pk: v for v in variant_qs.select_related("locus__contig", "locus__ref", "alt")}

        variants = []
        for row in top:
            variant = variants_by_id[row["variant_id"]]
            variants.append({
                "variant_id": variant.pk,
                "variant": str(variant),
                "count": row["count"],
            })

        tag_events = qs.count()
        top_events = sum(r["count"] for r in top)
        return {
            "tag": tag_id,
            "variants": variants,
            "tag_events": tag_events,
            "distinct_variants": qs.values("variant_id").distinct().count(),
            "top_events": top_events,
        }

    return _cached_json("re_tagged", genome_build, {"tag": tag_id}, calculate)


def tag_stats_tag_genes_over_time(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    tag_id = request.GET.get("tag") or _most_used_tag(genome_build)
    top_n = _get_int(request, "top_n", DEFAULT_TOP_GENES)

    def calculate():
        qs = _with_gene_symbol(genome_build, _tags_qs(genome_build).filter(tag=tag_id))
        qs = qs.annotate(quarter=TruncQuarter("created"))

        counts = {}
        quarters = set()
        for row in qs.values("quarter", GENE_SYMBOL_FIELD).annotate(count=Count("pk")):
            quarter = _quarter_label(row["quarter"])
            quarters.add(quarter)
            counts[(quarter, row[GENE_SYMBOL_FIELD])] = row["count"]

        sorted_quarters = sorted(quarters)
        return {
            "tag": tag_id,
            "quarters": sorted_quarters,
            "series": _grouped_series(counts, sorted_quarters, top_n),
        }

    return _cached_json("tag_genes_over_time", genome_build, {"tag": tag_id, "top_n": top_n}, calculate)


def _quarter_label(quarter_start) -> str:
    return f"{quarter_start.year}-Q{(quarter_start.month - 1) // 3 + 1}"


def _get_int(request, param: str, default: int) -> int:
    try:
        return int(request.GET[param])
    except (KeyError, ValueError):
        return default


def _get_list(request, param: str) -> list[str]:
    return [v for v in request.GET.getlist(param) if v]
