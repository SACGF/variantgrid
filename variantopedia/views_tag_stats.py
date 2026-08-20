"""
Tag stats page - @see https://github.com/SACGF/variantgrid/issues/1751

The page is a skeleton of cards, each lazy loading its own JSON endpoint so first paint is instant. Every
endpoint caches its response in Redis keyed by a hash of its parameters, and reports when it was calculated
so users can see how stale a card is.

Everything counts on allele, so the numbers don't depend on which builds a tag happens to have been lifted
over to - there's no genome build switch here. Tags that failed liftover have no allele at all, and are
reported as a separate count rather than being silently dropped.

Every card counts only the tags the requesting user can see, so the cache is keyed by user as well as by
the card's parameters.

Counts are labelled by what they count. Tagging repeats a lot (a known artefact gets re-tagged in every
analysis it turns up in) so "tag events", "distinct (allele, tag)" and "distinct alleles" are all very
different numbers.

An allele origin filter narrows every card to germline or somatic tags, so a curator can drop the half of
the vocabulary that isn't theirs. Tags marked as applying to both stay in either view, which is most of the
tagging. The filter is part of each card's cache key, and the page picks its starting value from the user's
allele origin settings.
"""
import hashlib
import json
from collections import Counter, defaultdict
from datetime import timedelta

from django.conf import settings
from django.contrib.auth.models import User
from django.core.cache import cache
from django.db.models import Count, F, QuerySet
from django.db.models.functions import TruncDate, TruncMonth, TruncQuarter
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, render
from django.utils.timezone import localtime, now

from analysis.models import VariantTag
from annotation.models import VariantAnnotationVersion
from library.constants import DAY_SECS
from snpdb.models import AlleleOriginFilterDefault, GenomeBuild, Lab, Tag, UserSettings, Variant, VariantAllele
from variantopedia.forms import TagStatsGenesForm, TagStatsTagForm, TagStatsTagsForm

TAG_STATS_CACHE_TIMEOUT = DAY_SECS
DEFAULT_TOP_TAGS = 8
DEFAULT_TOP_GENES = 10
TOP_USERS = 10
MEGA_ARTEFACT_VARIANTS = 40
OTHER = "other"
GENE_SYMBOL_FIELD = "variant__variantannotation__transcript_version__gene_version__gene_symbol_id"


def _cache_key(card: str, user: User, params: dict) -> str:
    params_hash = hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()[:16]
    return f"tag_stats/{settings.CACHE_VERSION}/{card}/{user.pk}/{params_hash}"


def _cached_json(card: str, user: User, allele_origin_filter: AlleleOriginFilterDefault, params: dict,
                 calculate) -> JsonResponse:
    """ The origin filter is folded in here rather than by each card, so no card can serve one origin's
        numbers to another origin's request """
    cache_key = _cache_key(card, user, {**params, "allele_origin": allele_origin_filter.value})
    data = cache.get(cache_key)
    if data is None:
        data = calculate()
        data["calculated"] = localtime(now()).isoformat()
        cache.set(cache_key, data, timeout=TAG_STATS_CACHE_TIMEOUT)
    return JsonResponse(data)


def _get_allele_origin_filter(request) -> AlleleOriginFilterDefault:
    """ The page always sends its radio, having resolved the initial value server side. Elsewhere (the user
        page's tagging section) there's nothing to widen the filter back out with, so absent means all """
    try:
        return AlleleOriginFilterDefault(request.GET["allele_origin"])
    except (KeyError, ValueError):
        return AlleleOriginFilterDefault.SHOW_ALL


def _tags_qs(user: User, allele_origin_filter: AlleleOriginFilterDefault) -> QuerySet[VariantTag]:
    """ Tags in analyses the user can't see are left out of every card, the same way the grids do it """
    qs = VariantTag.filter_for_user(user)
    if allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
        qs = qs.filter(tag__allele_origin_bucket__in=allele_origin_filter.buckets)
    return qs


def _tag_qs(allele_origin_filter: AlleleOriginFilterDefault) -> QuerySet[Tag]:
    """ Germline and Somatic each keep the tags marked as applying to both """
    qs = Tag.objects.all()
    if allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
        qs = qs.filter(allele_origin_bucket__in=allele_origin_filter.buckets)
    return qs


def _with_allele(tags_qs: QuerySet[VariantTag]) -> QuerySet[VariantTag]:
    """ Alleles are assigned asynchronously (@see analysis.tasks.variant_tag_tasks._liftover_variant_tag) so a
        tag has none if that's still to run, or if liftover failed """
    return tags_qs.exclude(allele__isnull=True)


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


def tag_stats(request):
    allele_origin_filter = _initial_allele_origin_filter(request.user)
    artefact_tag = _default_artefact_tag(request.user, allele_origin_filter)
    context = {
        # Only used to link out to the tagged variants grid, which shows coordinates so is build specific
        "genome_build": UserSettings.get_genome_build_or_default(request.user),
        "allele_origin_filter": allele_origin_filter.value,
        "allele_origin_choices": AlleleOriginFilterDefault.choices,
        "genes_form": TagStatsGenesForm(prefix="genes"),
        "co_occurrence_form": TagStatsTagsForm(prefix="co-occurrence"),
        "re_tagged_form": TagStatsTagForm(prefix="re-tagged", initial={"tag": artefact_tag}),
        "gene_time_form": TagStatsTagForm(prefix="gene-time",
                                          initial={"tag": _most_used_tag(request.user, allele_origin_filter)}),
    }
    return render(request, 'variantopedia/tag_stats.html', context)


def _initial_allele_origin_filter(user: User) -> AlleleOriginFilterDefault:
    """ Somebody who set an origin focus and asked for it to be applied automatically starts there - the same
        bargain classifications offer, driven by the same two settings """
    if allele_origin := UserSettings.get_for_user(user).default_allele_origin:
        return AlleleOriginFilterDefault(allele_origin)
    return AlleleOriginFilterDefault.SHOW_ALL


def _default_artefact_tag(user: User, allele_origin_filter: AlleleOriginFilterDefault) -> str:
    """ Artefact tagging is what drives the re-tagging numbers - fall back to the most used tag """
    if artefact := _tag_qs(allele_origin_filter).filter(pk__iexact="artefact").first():
        return artefact.pk
    return _most_used_tag(user, allele_origin_filter)


def _most_used_tag(user: User, allele_origin_filter: AlleleOriginFilterDefault) -> str:
    qs = _tags_qs(user, allele_origin_filter).values("tag_id").annotate(count=Count("pk")).order_by("-count")
    if top := qs.first():
        return top["tag_id"]
    return ""


def tag_stats_headline(request):
    allele_origin_filter = _get_allele_origin_filter(request)

    def calculate():
        qs = _tags_qs(request.user, allele_origin_filter)
        with_allele_qs = _with_allele(qs)
        year_ago = now() - timedelta(days=365)
        return {
            "tag_events": qs.count(),
            "distinct_allele_tags": with_allele_qs.values("allele_id", "tag_id").distinct().count(),
            "distinct_alleles": with_allele_qs.values("allele_id").distinct().count(),
            "no_allele": qs.filter(allele__isnull=True).count(),
            "analyses": qs.exclude(analysis__isnull=True).values("analysis_id").distinct().count(),
            "taggers_this_year": qs.filter(created__gte=year_ago).values("user_id").distinct().count(),
        }

    return _cached_json("headline", request.user, allele_origin_filter, {}, calculate)


def tag_stats_over_time(request):
    top_n = _get_int(request, "top_n", DEFAULT_TOP_TAGS)
    allele_origin_filter = _get_allele_origin_filter(request)

    def calculate():
        qs = _tags_qs(request.user, allele_origin_filter).annotate(month=TruncMonth("created"))
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

    return _cached_json("over_time", request.user, allele_origin_filter, {"top_n": top_n}, calculate)


def tag_stats_for_user(request, user_id=None):
    """ Shared by the stats page's "Your tagging" card and the user page section """
    user = get_object_or_404(User, pk=user_id) if user_id else request.user
    allele_origin_filter = _get_allele_origin_filter(request)

    def calculate():
        qs = _tags_qs(request.user, allele_origin_filter).filter(user=user)
        tag_counts = qs.values("tag_id").annotate(count=Count("pk")).order_by("-count")

        thirty_days_ago = (now() - timedelta(days=30)).date()
        daily = {row["day"]: row["count"] for row in
                 qs.filter(created__date__gte=thirty_days_ago).annotate(day=TruncDate("created"))
                 .values("day").annotate(count=Count("pk"))}
        days = [thirty_days_ago + timedelta(days=i) for i in range(31)]

        return {
            "username": str(user),
            "tag_events": qs.count(),
            "distinct_alleles": _with_allele(qs).values("allele_id").distinct().count(),
            "no_allele": qs.filter(allele__isnull=True).count(),
            "top_tags": [{"tag": r["tag_id"], "count": r["count"]} for r in tag_counts[:DEFAULT_TOP_TAGS]],
            "recent_days": [d.isoformat() for d in days],
            "recent_counts": [daily.get(d, 0) for d in days],
        }

    return _cached_json("user", request.user, allele_origin_filter, {"user_id": user.pk}, calculate)


def tag_stats_by_lab(request):
    allele_origin_filter = _get_allele_origin_filter(request)

    def calculate():
        qs = _tags_qs(request.user, allele_origin_filter)
        user_counts = qs.values("user_id", "user__username").annotate(count=Count("pk")).order_by("-count")
        thirty_days_ago = now() - timedelta(days=30)
        recent_counts = (qs.filter(created__gte=thirty_days_ago)
                         .values("user_id", "user__username").annotate(count=Count("pk")).order_by("-count"))

        labs_by_group_name = {lab.group_name: lab for lab in Lab.objects.exclude(group_name__isnull=True)}
        lab_names_by_user_id = defaultdict(list)
        for user_id, group_name in User.objects.filter(groups__name__in=labs_by_group_name).values_list(
                "pk", "groups__name"):
            if lab := labs_by_group_name.get(group_name):
                lab_names_by_user_id[user_id].append(str(lab))

        def user_rows(counts) -> list[dict]:
            return [{"user": row["user__username"],
                     "labs": lab_names_by_user_id.get(row["user_id"], []),
                     "count": row["count"]} for row in counts]

        # Lab totals still count everyone - only the user tables are capped
        lab_counts = Counter()
        for row in user_counts:
            for lab_name in lab_names_by_user_id.get(row["user_id"], []):
                lab_counts[lab_name] += row["count"]

        return {
            "top_users": user_rows(user_counts[:TOP_USERS]),
            "top_users_recent": user_rows(recent_counts[:TOP_USERS]),
            "labs": [{"lab": lab_name, "count": count} for lab_name, count in lab_counts.most_common()],
        }

    return _cached_json("by_lab", request.user, allele_origin_filter, {"top_users": TOP_USERS}, calculate)


def _with_gene_symbol(tags_qs: QuerySet[VariantTag]) -> QuerySet[VariantTag]:
    """ Restrict to tags whose variant has a representative transcript in the latest annotation for the build it
        was tagged in, so GENE_SYMBOL_FIELD can be grouped on without counting a lifted over tag twice """
    latest_versions = list(VariantAnnotationVersion.latest_for_all_builds())
    # A variant on a contig shared between builds has annotation in both, so pin the join to the tag's own build
    return tags_qs.filter(variant__variantannotation__version__in=latest_versions,
                          variant__variantannotation__version__genome_build=F("genome_build"),
                          variant__variantannotation__transcript_version__isnull=False)


def tag_stats_genes(request):
    gene_symbols = _get_list(request, "gene_symbols")
    tag_ids = _get_list(request, "tags")
    top_n = _get_int(request, "top_n", DEFAULT_TOP_GENES)
    allele_origin_filter = _get_allele_origin_filter(request)

    def calculate():
        tags_qs = _tags_qs(request.user, allele_origin_filter)
        if tag_ids:
            tags_qs = tags_qs.filter(tag__in=tag_ids)
        qs = _with_gene_symbol(tags_qs)

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
    return _cached_json("genes", request.user, allele_origin_filter, params, calculate)


def tag_stats_co_occurrence(request):
    allele_origin_filter = _get_allele_origin_filter(request)

    def calculate():
        qs = _with_allele(_tags_qs(request.user, allele_origin_filter))
        # Self join on allele - tag_id__gt gives each unordered pair once, and keeps a tag off its own diagonal.
        # Both conditions are in the one filter() so they apply to the same joined tag
        pairs_qs = (qs.filter(allele__varianttag__in=qs, allele__varianttag__tag_id__gt=F("tag_id"))
                    .values("tag_id", pair_tag_id=F("allele__varianttag__tag_id"))
                    .annotate(alleles=Count("allele_id", distinct=True)))
        # SQL emits each pair once, but orders the two tags by the DB collation - re-sort so keys match lookups
        pair_counts = Counter({tuple(sorted((row["tag_id"], row["pair_tag_id"]))): row["alleles"]
                               for row in pairs_qs})
        tag_totals = Counter({row["tag_id"]: row["alleles"] for row in
                              qs.values("tag_id").annotate(alleles=Count("allele_id", distinct=True))})

        tags = _top_names(tag_totals, DEFAULT_TOP_TAGS * 2)
        tags.sort()
        matrix = [[pair_counts.get(tuple(sorted((tag_a, tag_b))), 0) if tag_a != tag_b else None
                   for tag_b in tags] for tag_a in tags]
        return {
            "tags": tags,
            "matrix": matrix,
            "top_pairs": [{"tags": list(pair), "alleles": count} for pair, count in pair_counts.most_common(20)],
        }

    return _cached_json("co_occurrence", request.user, allele_origin_filter, {}, calculate)


def _variants_by_allele_id(allele_ids: list[int], genome_build: GenomeBuild) -> dict[int, Variant]:
    """ A variant to show per allele - the user's build where it's been lifted over, otherwise whatever we have """
    variants_by_allele_id = {}
    va_qs = VariantAllele.objects.filter(allele_id__in=allele_ids)
    for variant_allele in va_qs.select_related("variant__locus__contig", "variant__locus__ref", "variant__alt"):
        if variant_allele.allele_id not in variants_by_allele_id or variant_allele.genome_build_id == genome_build.pk:
            variants_by_allele_id[variant_allele.allele_id] = variant_allele.variant
    return variants_by_allele_id


def tag_stats_re_tagged(request):
    """ Alleles that get tagged over and over - the case for excluding them in analysis templates """
    allele_origin_filter = _get_allele_origin_filter(request)
    tag_id = request.GET.get("tag") or _default_artefact_tag(request.user, allele_origin_filter)
    genome_build = UserSettings.get_genome_build_or_default(request.user)

    def calculate():
        qs = _with_allele(_tags_qs(request.user, allele_origin_filter)).filter(tag=tag_id)
        top = list(qs.values("allele_id").annotate(count=Count("pk")).order_by("-count")[:MEGA_ARTEFACT_VARIANTS])
        variants_by_allele_id = _variants_by_allele_id([r["allele_id"] for r in top], genome_build)

        alleles = []
        for row in top:
            allele_id = row["allele_id"]
            variant = variants_by_allele_id.get(allele_id)
            alleles.append({
                "allele_id": allele_id,
                "variant": str(variant) if variant else f"Allele {allele_id}",
                "count": row["count"],
            })

        return {
            "tag": tag_id,
            "alleles": alleles,
            "tag_events": _tags_qs(request.user, allele_origin_filter).filter(tag=tag_id).count(),
            "distinct_alleles": qs.values("allele_id").distinct().count(),
            "top_events": sum(r["count"] for r in top),
        }

    return _cached_json("re_tagged", request.user, allele_origin_filter,
                        {"tag": tag_id, "genome_build": genome_build.name}, calculate)


def tag_stats_tag_genes_over_time(request):
    allele_origin_filter = _get_allele_origin_filter(request)
    tag_id = request.GET.get("tag") or _most_used_tag(request.user, allele_origin_filter)
    top_n = _get_int(request, "top_n", DEFAULT_TOP_GENES)

    def calculate():
        qs = _with_gene_symbol(_tags_qs(request.user, allele_origin_filter).filter(tag=tag_id))
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

    return _cached_json("tag_genes_over_time", request.user, allele_origin_filter,
                        {"tag": tag_id, "top_n": top_n}, calculate)


def _quarter_label(quarter_start) -> str:
    return f"{quarter_start.year}-Q{(quarter_start.month - 1) // 3 + 1}"


def _get_int(request, param: str, default: int) -> int:
    try:
        return int(request.GET[param])
    except (KeyError, ValueError):
        return default


def _get_list(request, param: str) -> list[str]:
    return [v for v in request.GET.getlist(param) if v]
