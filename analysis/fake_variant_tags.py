"""
Fake VariantTag data, so the tag stats page has enough to look at.
Run via 'manage.py create_fake_data tags' - @see https://github.com/SACGF/variantgrid/issues/1751

Everything created is obviously fake ('fake-artefact', 'fake_alice_somatic', 'Fake Somatic Unit') so it can't be
mistaken for real curation, and '--delete' takes it all away again.

The shape of the data matters more than the volume - the cards only look interesting if the numbers behave like
real tagging does: a handful of tags doing most of the work, tags introduced part way through, users coming and
going, genes with a long tail, and variants (artefacts especially) that get re-tagged over and over.

Tags are created with bulk_create, so none of the VariantTag signal handlers run - they're all about analyses and
liftover, and these tags have no analysis and get their allele assigned here.
"""
import random
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timedelta

from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.core.management.base import CommandError
from django.db import transaction
from django.utils.timezone import get_current_timezone, now
from guardian.models import GroupObjectPermission

from analysis.models import VariantTag
from analysis.models.enums import TagLocation
from annotation.fake_data import get_variant_ids_by_gene, zipf_weight
from library.guardian_utils import all_users_group
from snpdb.models import Allele, AlleleConversionTool, AlleleOrigin, GenomeBuild, GlobalSettings, Lab, Organization, \
    Tag, TagColor, TagColorsCollection, VariantAllele

FAKE_PREFIX = "fake"
FAKE_TAG_COLORS = "Fake demo tag colors"
SOMATIC = "somatic"
GERMLINE = "germline"
BOTH = "both"
BATCH_SIZE = 5000


@dataclass(frozen=True)
class FakeTag:
    tag_id: str
    rgb: str
    group: str
    weight: float
    introduced: float = 0.0  # How far through the 5 years the tag started being used
    genes: tuple[str, ...] = ()  # Tagged onto these genes far more often than the rest
    partners: tuple[str, ...] = ()  # Tags that end up on the same variant


FAKE_TAGS = [
    FakeTag("fake-artefact", "#bab0ac", SOMATIC, 22, genes=("TTN", "KIT", "EGFR"),
            partners=("fake-low-vaf", "fake-technical-check")),
    FakeTag("fake-low-vaf", "#9d755d", SOMATIC, 9, partners=("fake-artefact", "fake-clonal-haem")),
    FakeTag("fake-somatic-driver", "#e45756", SOMATIC, 11, genes=("KRAS", "NRAS", "BRAF", "TP53"),
            partners=("fake-actionable-tier1", "fake-resistance")),
    FakeTag("fake-actionable-tier1", "#f58518", SOMATIC, 7, genes=("EGFR", "BRAF", "KIT"),
            partners=("fake-somatic-driver",)),
    FakeTag("fake-resistance", "#b279a2", SOMATIC, 4, introduced=0.35, genes=("EGFR", "KIT", "FLT3"),
            partners=("fake-somatic-driver",)),
    FakeTag("fake-clonal-haem", "#eeca3b", SOMATIC, 6, introduced=0.5, genes=("DNMT3A", "TET2", "ASXL1", "SF3B1"),
            partners=("fake-low-vaf",)),
    FakeTag("fake-fusion-check", "#ff9da6", SOMATIC, 3, introduced=0.7, genes=("RUNX1", "JAK2"),
            partners=("fake-technical-check",)),
    FakeTag("fake-germline-pathogenic", "#4c78a8", GERMLINE, 10, genes=("BRCA1", "BRCA2", "MLH1", "MSH2"),
            partners=("fake-secondary-finding", "fake-review-me")),
    FakeTag("fake-vus-review", "#54a24b", GERMLINE, 14, genes=("ATM", "CHEK2", "PALB2", "TTN"),
            partners=("fake-review-me", "fake-carrier")),
    FakeTag("fake-carrier", "#72b7b2", GERMLINE, 8, genes=("CFTR", "DMD", "RYR1"),
            partners=("fake-vus-review",)),
    FakeTag("fake-secondary-finding", "#6a51a3", GERMLINE, 4, introduced=0.25, genes=("LDLR", "MYBPC3", "APC"),
            partners=("fake-germline-pathogenic",)),
    FakeTag("fake-review-me", "#495057", BOTH, 8, partners=("fake-vus-review", "fake-germline-pathogenic")),
    FakeTag("fake-technical-check", "#17becf", BOTH, 5, introduced=0.15,
            partners=("fake-artefact", "fake-fusion-check")),
]

# Real gene symbols (they need annotated variants to join against) in rough order of how often they're tagged
SOMATIC_GENES = ["TP53", "KRAS", "NRAS", "BRAF", "EGFR", "PIK3CA", "KIT", "JAK2", "DNMT3A", "TET2",
                 "ASXL1", "RUNX1", "IDH1", "IDH2", "FLT3", "NPM1", "SF3B1", "CALR", "PTEN", "SMAD4"]
GERMLINE_GENES = ["BRCA1", "BRCA2", "ATM", "MSH2", "MLH1", "MSH6", "PMS2", "PALB2", "CHEK2", "APC",
                  "CFTR", "TTN", "MYBPC3", "RYR1", "DMD", "FBN1", "NF1", "LDLR", "SCN1A", "COL4A5"]


@dataclass(frozen=True)
class FakeLab:
    name: str
    group_name: str
    city: str


FAKE_LABS = {
    SOMATIC: FakeLab("Fake Somatic Unit", f"{FAKE_PREFIX}_health/{FAKE_PREFIX}_somatic_unit", "Faketown"),
    GERMLINE: FakeLab("Fake Germline Unit", f"{FAKE_PREFIX}_health/{FAKE_PREFIX}_germline_unit", "Faketown"),
}


@dataclass(frozen=True)
class FakeUser:
    username: str
    group: str
    weight: float
    active: tuple[float, float] = (0.0, 1.0)  # Fraction of the 5 years they were around for


FAKE_USERS = [
    FakeUser("fake_alice_somatic", SOMATIC, 10),
    FakeUser("fake_bruno_somatic", SOMATIC, 7),
    FakeUser("fake_chen_somatic", SOMATIC, 6, active=(0.2, 1.0)),
    FakeUser("fake_dina_somatic", SOMATIC, 4, active=(0.0, 0.6)),
    FakeUser("fake_eli_somatic", SOMATIC, 3, active=(0.55, 1.0)),
    FakeUser("fake_freya_germline", GERMLINE, 10),
    FakeUser("fake_gus_germline", GERMLINE, 6),
    FakeUser("fake_hana_germline", GERMLINE, 5, active=(0.3, 1.0)),
    FakeUser("fake_ivan_germline", GERMLINE, 4, active=(0.0, 0.75)),
    FakeUser("fake_jo_germline", GERMLINE, 3, active=(0.6, 1.0)),
    FakeUser("fake_rotating_registrar", BOTH, 2, active=(0.8, 1.0)),
    FakeUser("fake_import_bot", BOTH, 3, active=(0.1, 1.0)),
]

# Every variant gets one, then the extras are handed out by weight - a few variants soak up hundreds of them
MEGA_ARTEFACT_VARIANTS = 40
MEGA_ARTEFACT_WEIGHT = 400
LOCATION_WEIGHTS = {
    TagLocation.VARIANT_DETAILS: 6,
    TagLocation.GENE_PAGE: 2,
    TagLocation.EXTERNAL: 1,
}


@dataclass
class FakeVariant:
    variant_id: int
    group: str
    gene_symbol: str
    primary_tag: FakeTag
    events: int = 1
    weight: float = 1.0
    allele_id: int = field(default=0)


class FakeVariantTags:
    HELP = "Variant tags, with the re-tagging and long tails the tag stats page is built to show"

    def __init__(self, stdout):
        self.stdout = stdout

    @staticmethod
    def add_arguments(parser):
        parser.add_argument("--genome-build", default="GRCh37")
        parser.add_argument("--somatic-events", type=int, default=100_000)
        parser.add_argument("--germline-events", type=int, default=100_000)
        parser.add_argument("--events-per-variant", type=float, default=3.0)
        parser.add_argument("--years", type=int, default=5)
        parser.add_argument("--seed", type=int, default=1751)
        parser.add_argument("--also-tag-as", help="Comma separated existing usernames to mix into the taggers, "
                                                  "so their 'Your tagging' card has something in it")

    def create(self, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        random.seed(options["seed"])
        self.stdout.write(f"Creating fake tags in {genome_build}")

        users = self._create_users(options["also_tag_as"])
        self._create_tags()
        months = _months(options["years"])

        fake_variants = []
        for group, genes, num_events in ((SOMATIC, SOMATIC_GENES, options["somatic_events"]),
                                         (GERMLINE, GERMLINE_GENES, options["germline_events"])):
            num_variants = int(num_events / options["events_per_variant"])
            group_variants = self._pick_variants(genome_build, group, genes, num_variants)
            _allocate_events(group_variants, num_events)
            fake_variants.extend(group_variants)
            self.stdout.write(f"{group}: {num_events} tag events over {len(group_variants)} variants")

        self._create_alleles(genome_build, fake_variants)
        variant_tags = _build_variant_tags(genome_build, fake_variants, users, months)
        self._save_variant_tags(variant_tags)

    def _create_users(self, also_tag_as: str) -> list[tuple[User, FakeUser]]:
        organization, _ = Organization.objects.get_or_create(
            group_name=f"{FAKE_PREFIX}_health", defaults={"name": "Fake Health Network", "short_name": "FakeHealth"})

        labs_by_group = {}
        for group, fake_lab in FAKE_LABS.items():
            lab, _ = Lab.objects.get_or_create(group_name=fake_lab.group_name,
                                               defaults={"name": fake_lab.name, "city": fake_lab.city,
                                                         "organization": organization})
            labs_by_group[group] = Group.objects.get_or_create(name=lab.group_name)[0]

        users = []
        for fake_user in FAKE_USERS:
            user, created = User.objects.get_or_create(username=fake_user.username,
                                                       defaults={"first_name": "Fake", "is_active": False})
            if created:
                user.groups.add(all_users_group())
                if lab_group := labs_by_group.get(fake_user.group):
                    user.groups.add(lab_group)
            users.append((user, fake_user))

        for username in filter(None, (also_tag_as or "").split(",")):
            user = User.objects.filter(username=username.strip()).first()
            if not user:
                raise CommandError(f"--also-tag-as: no such user '{username}'")
            users.append((user, FakeUser(user.username, BOTH, 6)))

        return users

    def _create_tags(self):
        collection = _tag_colors_collection()
        for fake_tag in FAKE_TAGS:
            tag, _ = Tag.objects.get_or_create(pk=fake_tag.tag_id)
            TagColor.objects.update_or_create(collection=collection, tag=tag, defaults={"rgb": fake_tag.rgb})

    def _pick_variants(self, genome_build: GenomeBuild, group: str, genes: list[str],
                       num_variants: int) -> list[FakeVariant]:
        """ Real variants, so the gene cards have real symbols to group on and the links all work """
        variant_ids_by_gene = get_variant_ids_by_gene(genome_build, genes)
        tags = [t for t in FAKE_TAGS if t.group in (group, BOTH)]

        fake_variants = []
        for i, gene_symbol in enumerate(genes):
            available = variant_ids_by_gene.get(gene_symbol, [])
            wanted = int(num_variants * zipf_weight(i) / sum(zipf_weight(j) for j in range(len(genes))))
            if len(available) < wanted:
                self.stdout.write(f"{gene_symbol}: only {len(available)} annotated variants, wanted {wanted}")
                wanted = len(available)
            for variant_id in random.sample(available, wanted):
                weights = [t.weight * (6 if gene_symbol in t.genes else 1) for t in tags]
                primary_tag = random.choices(tags, weights=weights)[0]
                weight = random.lognormvariate(0, 1.1)
                fake_variants.append(FakeVariant(variant_id, group, gene_symbol, primary_tag, weight=weight))

        artefacts = [fv for fv in fake_variants if fv.primary_tag.tag_id == "fake-artefact"]
        for fake_variant in random.sample(artefacts, min(len(artefacts), MEGA_ARTEFACT_VARIANTS)):
            fake_variant.weight = MEGA_ARTEFACT_WEIGHT
        return fake_variants

    def _create_alleles(self, genome_build: GenomeBuild, fake_variants: list[FakeVariant]):
        """ Tags carry an allele - the co-occurrence card and the tagged variants grid both join through it """
        alleles = Allele.objects.bulk_create([Allele() for _ in fake_variants], batch_size=BATCH_SIZE)
        variant_alleles = []
        for fake_variant, allele in zip(fake_variants, alleles):
            fake_variant.allele_id = allele.pk
            variant_alleles.append(VariantAllele(variant_id=fake_variant.variant_id, genome_build=genome_build,
                                                 allele=allele, origin=AlleleOrigin.IMPORTED_TO_DATABASE,
                                                 allele_linking_tool=AlleleConversionTool.SAME_CONTIG))
        VariantAllele.objects.bulk_create(variant_alleles, batch_size=BATCH_SIZE)
        self.stdout.write(f"Created {len(alleles)} alleles")

    def _save_variant_tags(self, variant_tags: list[VariantTag]):
        with _keeping_our_timestamps():
            VariantTag.objects.bulk_create(variant_tags, batch_size=BATCH_SIZE)
        self.stdout.write(f"Created {len(variant_tags)} variant tags")

        # bulk_create skips GuardianPermissionsAutoInitialSaveMixin.save() - everyone can see (and clean up)
        # fake data, which is all the permission it needs
        content_type = ContentType.objects.get_for_model(VariantTag)
        group = all_users_group()
        permissions = Permission.objects.filter(content_type=content_type,
                                                codename__in=[VariantTag.get_read_perm(), VariantTag.get_write_perm()])
        num_permissions = 0
        for permission in permissions:
            for i in range(0, len(variant_tags), BATCH_SIZE):
                GroupObjectPermission.objects.bulk_create(
                    [GroupObjectPermission(group=group, permission=permission, content_type=content_type,
                                           object_pk=str(vt.pk)) for vt in variant_tags[i:i + BATCH_SIZE]])
                num_permissions += len(variant_tags[i:i + BATCH_SIZE])
        self.stdout.write(f"Created {num_permissions} group permissions")

    @transaction.atomic
    def delete(self, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        tag_ids = [ft.tag_id for ft in FAKE_TAGS]
        variant_tags_qs = VariantTag.objects.filter(tag__in=tag_ids)
        allele_ids = set(variant_tags_qs.values_list("allele_id", flat=True))

        content_type = ContentType.objects.get_for_model(VariantTag)
        object_pks = [str(pk) for pk in variant_tags_qs.values_list("pk", flat=True)]
        deleted_permissions = 0
        for i in range(0, len(object_pks), BATCH_SIZE):
            deleted, _ = GroupObjectPermission.objects.filter(
                content_type=content_type, object_pk__in=object_pks[i:i + BATCH_SIZE]).delete()
            deleted_permissions += deleted
        deleted_tags, _ = variant_tags_qs.delete()

        # Only the alleles we made - ie no ClinGen record, and nothing else left pointing at them
        alleles_qs = Allele.objects.filter(pk__in=allele_ids, clingen_allele__isnull=True,
                                           flag_collection__isnull=True, varianttag__isnull=True)
        VariantAllele.objects.filter(allele__in=alleles_qs, genome_build=genome_build).delete()
        deleted_alleles, _ = alleles_qs.delete()

        Tag.objects.filter(pk__in=tag_ids).delete()  # Cascades to TagColor
        TagColorsCollection.objects.filter(name=FAKE_TAG_COLORS).delete()
        User.objects.filter(username__in=[fu.username for fu in FAKE_USERS]).delete()
        Lab.objects.filter(group_name__in=[fl.group_name for fl in FAKE_LABS.values()]).delete()
        Organization.objects.filter(group_name=f"{FAKE_PREFIX}_health").delete()
        Group.objects.filter(name__startswith=f"{FAKE_PREFIX}_health").delete()
        self.stdout.write(f"Deleted {deleted_tags} variant tags, {deleted_permissions} permissions, "
                          f"{deleted_alleles} alleles, and the fake tags/users/labs")


@contextmanager
def _keeping_our_timestamps():
    """ created/modified are auto_now_add/auto_now, which would stamp all 5 years of tagging as happening now """
    created = VariantTag._meta.get_field("created")
    modified = VariantTag._meta.get_field("modified")
    created.auto_now_add = False
    modified.auto_now = False
    try:
        yield
    finally:
        created.auto_now_add = True
        modified.auto_now = True


def _tag_colors_collection() -> TagColorsCollection:
    global_settings = GlobalSettings.objects.first()
    if global_settings and global_settings.tag_colors:
        return global_settings.tag_colors

    collection = TagColorsCollection.objects.filter(name=FAKE_TAG_COLORS).first()
    if collection is None:
        collection = TagColorsCollection(name=FAKE_TAG_COLORS)
        collection.save(assign_permissions=False)  # Global collection - no user to assign to
    if global_settings:
        global_settings.tag_colors = collection
        global_settings.save()
    return collection


def _allocate_events(fake_variants: list[FakeVariant], num_events: int):
    """ Everything is tagged once, then the rest are dealt out by weight - the tail is where re-tagging lives """
    extra = num_events - len(fake_variants)
    if extra <= 0:
        return
    weights = [fv.weight for fv in fake_variants]
    for fake_variant in random.choices(fake_variants, weights=weights, k=extra):
        fake_variant.events += 1


def _months(years: int) -> list[datetime]:
    """ Start of each month in the period, oldest first """
    today = now().astimezone(get_current_timezone()).date()
    month_starts = []
    year, month = today.year, today.month
    for _ in range(years * 12):
        month_starts.append(datetime(year, month, 1, tzinfo=get_current_timezone()))
        year, month = (year - 1, 12) if month == 1 else (year, month - 1)
    return sorted(month_starts)


def _month_weights(months: list[datetime]) -> list[float]:
    """ Tagging grew over the period, with a quiet patch over the holidays and a bit of noise """
    weights = []
    for i, month in enumerate(months):
        weight = 0.35 + 0.65 * i / len(months)
        if month.month in (1, 12):
            weight *= 0.6
        weights.append(weight * random.uniform(0.8, 1.2))
    return weights


def _cumulative(weights: list[float]) -> list[float]:
    total = 0.0
    cumulative = []
    for weight in weights:
        total += weight
        cumulative.append(total)
    return cumulative


def _random_working_datetime(month_start: datetime, latest: datetime) -> datetime:
    """ A weekday during working hours in the month - the current month only runs up until now """
    days = max(1, min(28, (latest - month_start).days))
    day = month_start + timedelta(days=random.randrange(days), hours=random.randrange(8, 18),
                                  minutes=random.randrange(60), seconds=random.randrange(60))
    while day.weekday() >= 5:
        day -= timedelta(days=1)
    return day


def _build_variant_tags(genome_build: GenomeBuild, fake_variants: list[FakeVariant],
                        users: list[tuple[User, FakeUser]], months: list[datetime]) -> list[VariantTag]:
    latest = now()
    month_weights = _month_weights(months)
    tags_by_id = {ft.tag_id: ft for ft in FAKE_TAGS}

    # A tag can only be used from the month it was introduced, so each has its own slice of the timeline
    month_indexes = {month: i for i, month in enumerate(months)}
    months_by_tag = {}
    for fake_tag in FAKE_TAGS:
        start = int(fake_tag.introduced * len(months))
        months_by_tag[fake_tag.tag_id] = (months[start:], _cumulative(month_weights[start:]))

    # Users come and go, so who was around depends on which month we land in
    users_by_month = {}
    for i in range(len(months)):
        fraction = i / len(months)
        active = [(user, fake_user) for user, fake_user in users
                  if fake_user.active[0] <= fraction <= fake_user.active[1]]
        users_by_month[i] = (active, _cumulative([fu.weight for _, fu in active]))

    locations = list(LOCATION_WEIGHTS)
    location_cumulative = _cumulative(list(LOCATION_WEIGHTS.values()))
    group_tags = {group: [ft for ft in FAKE_TAGS if ft.group in (group, BOTH)] for group in (SOMATIC, GERMLINE)}
    group_tag_cumulative = {group: _cumulative([ft.weight for ft in tags]) for group, tags in group_tags.items()}

    variant_tags = []
    for fake_variant in fake_variants:
        primary = fake_variant.primary_tag
        for _ in range(fake_variant.events):
            roll = random.random()
            if roll < 0.72 or not primary.partners:
                fake_tag = primary
            elif roll < 0.9:
                fake_tag = tags_by_id[random.choice(primary.partners)]
            else:
                fake_tag = random.choices(group_tags[fake_variant.group],
                                          cum_weights=group_tag_cumulative[fake_variant.group])[0]

            tag_months, tag_cumulative = months_by_tag[fake_tag.tag_id]
            month = random.choices(tag_months, cum_weights=tag_cumulative)[0]
            active_users, user_cumulative = users_by_month[month_indexes[month]]
            user, _ = random.choices(active_users, cum_weights=user_cumulative)[0]

            created = _random_working_datetime(month, latest)
            variant_tags.append(VariantTag(variant_id=fake_variant.variant_id, allele_id=fake_variant.allele_id,
                                           genome_build=genome_build, tag_id=fake_tag.tag_id, user=user,
                                           location=random.choices(locations, cum_weights=location_cumulative)[0],
                                           created=created, modified=created))
    return variant_tags
