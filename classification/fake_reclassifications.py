"""
Fake germline classifications with a curation history, so the reclassification analytics page has
something to render. Run via 'manage.py create_fake_data reclassifications' - @see
https://github.com/SACGF/variantgrid/issues/1790

Everything created is obviously fake (labs and users are all 'Fake ...' / 'fake_...', and every record's
lab_record_id starts with 'fake-recl-') so it can't be mistaken for real curation, and '--delete' takes it
all away again.

The page reads ReclassificationEvent, which is derived from published modification history rather than
written directly, so this builds the history the derivation expects - a chain of published
ClassificationModifications per record, each with its own clinical_significance, curation_date and ACMG
criteria - and then hands it to ReclassificationEventBuilder like any real import would.

An allele is what the classification listing groups on, so every record gets one. Germline records get an
allele each; somatic records are dealt out in clusters that share an allele and a lab, which is exactly what
a classification grouping is - so the listing has groupings holding several records to render. Somatic
records are curated on tiers rather than the germline scale, and the reclassification page leaves them out
of its timelines the same way it leaves out any somatic record.

The shape matters more than the volume, because every chart on the page is about how curation behaves:
records mostly move one bucket at a time, VUS is where nearly all the movement is, P and B hardly ever
budge, most reviews confirm the call rather than change it, labs revisit at different rates, and the
criteria that flip are the ones that would actually justify the direction travelled.
"""
import random
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta
from typing import Optional

from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.db import transaction
from django.utils.timezone import get_current_timezone, now
from guardian.models import GroupObjectPermission

from annotation.fake_data import get_variant_ids_by_gene, zipf_weight
from classification.enums import (AlleleOriginBucket, ClinicalSignificance, CriteriaEvaluation, ShareLevel,
                                  SpecialEKeys, SubmissionSource)
from classification.models import (Classification, ClassificationModification, ClassificationSummaryCalculator,
                                   EvidenceKeyMap, ImportedAlleleInfo, ReclassificationEvent,
                                   ReclassificationEventBuilder, ResolvedVariantInfo)
from classification.models.classification_grouping import (AlleleOriginGrouping, ClassificationGrouping,
                                                           ClassificationGroupingEntry)
from classification.models.classification_variant_info_models import ImportedAlleleInfoStatus
from library.guardian_utils import all_users_group
from snpdb.models import (Allele, AlleleConversionTool, AlleleOrigin, GenomeBuild, GenomeBuildPatchVersion, Lab,
                          Organization, VariantAllele)

FAKE_PREFIX = "fake"
FAKE_ORGANIZATION_GROUP = f"{FAKE_PREFIX}_curation_network"
LAB_RECORD_PREFIX = "fake-recl-"
SOMATIC_LAB_RECORD_PREFIX = f"{LAB_RECORD_PREFIX}som-"
BATCH_SIZE = 2000

# The page's survival curve takes its cohort five years back and follows it forward, so the history has to
# start comfortably before that and the initial classifications have to be front loaded
SURVIVAL_ORIGIN_YEARS_BACK = 5
INITIAL_SPREAD = 0.35
""" Initial classifications land in this fraction of the period, so there's a cohort to follow """


@dataclass(frozen=True)
class FakeCurationLab:
    name: str
    weight: float
    """ Share of the classifications """
    review_rate: float
    """ Multiplier on how often this lab revisits a record - drives the spread on the by-lab chart """
    move_rate: float = 1.0
    """ Multiplier on how often one of its reviews changes the call """
    external: bool = False

    @property
    def group_name(self) -> str:
        slug = self.name.lower().replace(" ", "_")
        return f"{FAKE_ORGANIZATION_GROUP}/{slug}"


FAKE_LABS = [
    FakeCurationLab("Fake Diagnostic Genetics", weight=10, review_rate=1.0),
    FakeCurationLab("Fake Molecular Pathology", weight=6, review_rate=1.7, move_rate=0.8),
    FakeCurationLab("Fake Rare Disease Service", weight=4, review_rate=0.5, move_rate=1.4),
    FakeCurationLab("Fake Partner Laboratory", weight=3, review_rate=0.8, external=True),
]

FAKE_CURATORS = [
    "fake_curator_nadia",
    "fake_curator_omar",
    "fake_curator_priya",
    "fake_curator_quinn",
    "fake_curator_reuben",
]


@dataclass(frozen=True)
class FakeSignificanceBehaviour:
    """ How a record sitting at one clinical significance behaves when a curator comes back to it """
    weight: float
    """ Share of the initial classifications that start here """
    move_chance: float
    """ Chance a review moves the call rather than confirming it """
    benign_bias: float
    """ Of the moves, the share heading towards benign """


# VUS is where the curation effort goes and where nearly all the movement is; the ends of the scale are
# where labs are most confident, so they hardly ever move and mostly get confirmed instead
SIGNIFICANCE_BEHAVIOUR = {
    ClinicalSignificance.PATHOGENIC: FakeSignificanceBehaviour(8, move_chance=0.06, benign_bias=0.95),
    ClinicalSignificance.LIKELY_PATHOGENIC: FakeSignificanceBehaviour(10, move_chance=0.38, benign_bias=0.30),
    ClinicalSignificance.VUS: FakeSignificanceBehaviour(55, move_chance=0.30, benign_bias=0.58),
    ClinicalSignificance.LIKELY_BENIGN: FakeSignificanceBehaviour(15, move_chance=0.32, benign_bias=0.72),
    ClinicalSignificance.BENIGN: FakeSignificanceBehaviour(12, move_chance=0.05, benign_bias=0.05),
}

SIGNIFICANCE_WEIGHTS = {significance: behaviour.weight for significance, behaviour in SIGNIFICANCE_BEHAVIOUR.items()}

# Ordered benign to pathogenic, so a step of one is a move to the adjacent bucket
SIGNIFICANCE_SCALE = [
    ClinicalSignificance.BENIGN,
    ClinicalSignificance.LIKELY_BENIGN,
    ClinicalSignificance.VUS,
    ClinicalSignificance.LIKELY_PATHOGENIC,
    ClinicalSignificance.PATHOGENIC,
]

# Real gene symbols (they need annotated variants to join against), longest tail first - the VUS burden
# chart is only interesting if a few genes carry far more of it than the rest
GENES = ["TTN", "BRCA2", "BRCA1", "ATM", "APC", "MSH6", "NF1", "PMS2", "MLH1", "MSH2",
         "PALB2", "CHEK2", "JAK2", "RYR1", "CFTR", "DMD", "FBN1", "MYBPC3", "SCN1A", "LDLR",
         "COL4A5", "TP53", "PTEN", "STK11", "CDH1", "RB1", "VHL", "RET", "SDHB", "MUTYH", "BMPR1A"]

# Somatic curation runs on tiers instead, and most of what a lab sees sits in the tiers it can't act on
SOMATIC_TIER_WEIGHTS = {"tier_1": 8, "tier_1_or_2": 12, "tier_2": 20, "tier_3": 45, "tier_4": 15}
# Ordered least to most clinically actionable - tier_1_or_2 is off the scale, a review resolves it instead
SOMATIC_TIER_SCALE = ["tier_4", "tier_3", "tier_2", "tier_1"]
SOMATIC_MOVE_CHANCE = 0.25
""" Chance a somatic review moves the tier - it takes a new trial or guideline, not just a re-read """
SOMATIC_TOWARDS_ACTIONABLE = 0.7
""" Of the tier moves, the share heading towards tier 1 - evidence accumulates more often than it evaporates """
SOMATIC_REVIEWS_PER_RECORD = 1.4

# The AMP level a tier gets reported at, which together with the tier is what the somatic sort is built from
TIER_AMP_LEVELS = {
    "tier_1": ["amp:level_a", "amp:level_b"],
    "tier_2": ["amp:level_c", "amp:level_d"],
}
AMP_CONTEXTS = ["therapeutic", "diagnostic", "prognostic"]

# The oncogenicity call that goes with a tier - the actionable tiers are the ones that are oncogenic
ONCOGENICITY_FOR_TIER = {
    "tier_1": {"O": 8, "LO": 2},
    "tier_1_or_2": {"O": 5, "LO": 5},
    "tier_2": {"O": 3, "LO": 7},
    "tier_3": {"VUS": 1},
    "tier_4": {"LB": 6, "B": 4},
}

# Somatic curation is against the tumour rather than an inherited condition
FAKE_TUMOUR_TYPES = ["Fake melanoma", "Fake colorectal carcinoma", "Fake acute myeloid leukaemia",
                     "Fake non-small cell lung carcinoma", "Fake glioblastoma", "Fake breast carcinoma"]

PATHOGENIC_CRITERIA = {
    "acmg:pvs1": CriteriaEvaluation.PATHOGENIC_VERY_STRONG,
    "acmg:ps1": CriteriaEvaluation.PATHOGENIC_STRONG,
    "acmg:ps3": CriteriaEvaluation.PATHOGENIC_STRONG,
    "acmg:pm1": CriteriaEvaluation.PATHOGENIC_MODERATE,
    "acmg:pm2": CriteriaEvaluation.PATHOGENIC_MODERATE,
    "acmg:pm5": CriteriaEvaluation.PATHOGENIC_MODERATE,
    "acmg:pp1": CriteriaEvaluation.PATHOGENIC_SUPPORTING,
    "acmg:pp3": CriteriaEvaluation.PATHOGENIC_SUPPORTING,
}
BENIGN_CRITERIA = {
    "acmg:ba1": CriteriaEvaluation.BENIGN_STANDALONE,
    "acmg:bs1": CriteriaEvaluation.BENIGN_STRONG,
    "acmg:bs3": CriteriaEvaluation.BENIGN_STRONG,
    "acmg:bs4": CriteriaEvaluation.BENIGN_STRONG,
    "acmg:bp1": CriteriaEvaluation.BENIGN_SUPPORTING,
    "acmg:bp4": CriteriaEvaluation.BENIGN_SUPPORTING,
    "acmg:bp7": CriteriaEvaluation.BENIGN_SUPPORTING,
}
# A criterion can be met more weakly or more strongly than its default, which is what tells the evidence
# chart's "strengthened" and "weakened" bars apart from a plain on/off
WEAKER_STRENGTH = {
    CriteriaEvaluation.PATHOGENIC_VERY_STRONG: CriteriaEvaluation.PATHOGENIC_STRONG,
    CriteriaEvaluation.PATHOGENIC_STRONG: CriteriaEvaluation.PATHOGENIC_MODERATE,
    CriteriaEvaluation.PATHOGENIC_MODERATE: CriteriaEvaluation.PATHOGENIC_SUPPORTING,
    CriteriaEvaluation.BENIGN_STANDALONE: CriteriaEvaluation.BENIGN_STRONG,
    CriteriaEvaluation.BENIGN_STRONG: CriteriaEvaluation.BENIGN_SUPPORTING,
}
STRONGER_STRENGTH = {weaker: stronger for stronger, weaker in WEAKER_STRENGTH.items()}

# The handful of non-criteria keys a curator actually rewrites when they revisit a record. The evidence
# chart groups everything outside the criteria, so these are what fills the non-criteria rows
NARRATIVE_KEYS = ["literature", "interpretation_summary", "gnomad_af", "condition", "variant_type",
                  "molecular_consequence", "patient_phenotype", "case_details"]
SOMATIC_NARRATIVE_KEYS = ["somatic:summary_interpretation", "literature", "interpretation_summary",
                          "molecular_consequence"]

REVIEWS_PER_RECORD = 2.6
""" Average times a record is revisited, before the lab's own review_rate is applied """


@dataclass
class FakeRecord:
    """ One classification and the curation history that will be published for it """
    lab: FakeCurationLab
    lab_obj: Lab
    variant_id: int
    gene_symbol: str
    started: datetime
    somatic: bool = False
    significance: Optional[str] = None
    """ Germline: the initial call the history walks forward from """
    tier: Optional[str] = None
    """ Somatic: the initial tier the history walks forward from """
    steps: list['FakeStep'] = field(default_factory=list)
    allele_info_id: Optional[int] = None
    classification_id: Optional[int] = None


@dataclass
class FakeAlleleGroup:
    """ The records that resolved to the one allele - a somatic allele carries a lab's whole cluster """
    records: list[FakeRecord]

    @property
    def variant_id(self) -> int:
        return self.records[0].variant_id

    @property
    def gene_symbol(self) -> str:
        return self.records[0].gene_symbol


@dataclass
class FakeStep:
    """ One published modification - the call as it stood, and the evidence behind it """
    significance: str
    published_on: datetime
    curated_on: date
    evidence: dict


class FakeReclassifications:
    HELP = ("Classifications with a curation history - germline for the reclassification analytics page, "
            "somatic clustered onto shared alleles")

    def __init__(self, stdout):
        self.stdout = stdout

    @staticmethod
    def add_arguments(parser):
        parser.add_argument("--genome-build", default="GRCh38")
        parser.add_argument("--classifications", type=int, default=500,
                            help="Germline classifications, one allele each")
        parser.add_argument("--somatic-classifications", type=int, default=200,
                            help="Somatic classifications, clustered onto shared alleles")
        parser.add_argument("--somatic-per-grouping", type=int, default=5,
                            help="Average somatic classifications sharing an allele at the one lab")
        parser.add_argument("--somatic-spread", type=int, default=3,
                            help="How far either side of that average a cluster can land")
        parser.add_argument("--years", type=int, default=7,
                            help="History spans this many years back from today")
        parser.add_argument("--adjacent-percent", type=float, default=85,
                            help="Share of moves that go to the neighbouring bucket rather than further")
        parser.add_argument("--seed", type=int, default=1790)

    def create(self, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        random.seed(options["seed"])
        self.stdout.write(f"Creating {options['classifications']} germline and about "
                          f"{options['somatic_classifications']} somatic fake classifications in {genome_build}")

        labs, curators = self._create_labs_and_curators()
        groups = self._pick_groups(genome_build, labs, options)
        records = [record for group in groups for record in group.records]
        _build_histories(records, options["adjacent_percent"] / 100)

        self._create_allele_infos(genome_build, groups)
        self._create_classifications(records, curators)
        modifications = self._create_modifications(records, curators)
        self._assign_permissions(records, modifications)
        self._summarise_and_group(records)
        self._build_timelines(records)

    def _create_labs_and_curators(self) -> tuple[dict[str, Lab], list[User]]:
        organization, _ = Organization.objects.get_or_create(
            group_name=FAKE_ORGANIZATION_GROUP,
            defaults={"name": "Fake Curation Network", "short_name": "FakeCuration"})

        labs = {}
        lab_groups = []
        for fake_lab in FAKE_LABS:
            lab, _ = Lab.objects.get_or_create(
                group_name=fake_lab.group_name,
                defaults={"name": fake_lab.name, "city": "Faketown", "organization": organization,
                          "external": fake_lab.external})
            labs[fake_lab.name] = lab
            lab_groups.append(Group.objects.get_or_create(name=lab.group_name)[0])

        curators = []
        for username in FAKE_CURATORS:
            user, created = User.objects.get_or_create(username=username,
                                                       defaults={"first_name": "Fake", "is_active": False})
            if created:
                user.groups.add(all_users_group(), *lab_groups)
            curators.append(user)
        return labs, curators

    def _pick_groups(self, genome_build: GenomeBuild, labs: dict[str, Lab], options: dict) -> list[FakeAlleleGroup]:
        """ Real variants, so the gene chart has real symbols and the records link somewhere sensible """
        variant_ids_by_gene = get_variant_ids_by_gene(genome_build, GENES, without_alleles=True)
        years = options["years"]
        period_start = _period_start(years)
        initial_window = INITIAL_SPREAD * years * 365
        per_grouping = max(1, options["somatic_per_grouping"])
        somatic_clusters = round(options["somatic_classifications"] / per_grouping)
        total_weight = sum(zipf_weight(i) for i in range(len(GENES)))

        groups = []
        for index, gene_symbol in enumerate(GENES):
            gene_share = zipf_weight(index) / total_weight
            available = variant_ids_by_gene.get(gene_symbol, [])
            germline_wanted = round(options["classifications"] * gene_share)
            wanted = germline_wanted + round(somatic_clusters * gene_share)
            if len(available) < wanted:
                self.stdout.write(f"{gene_symbol}: only {len(available)} annotated variants, wanted {wanted}")
                wanted = len(available)
                germline_wanted = min(germline_wanted, wanted)
            chosen = random.sample(available, wanted)

            for variant_id in chosen[:germline_wanted]:
                groups.append(FakeAlleleGroup([
                    _germline_record(labs, gene_symbol, variant_id, _started(period_start, initial_window))]))
            for variant_id in chosen[germline_wanted:]:
                groups.append(_somatic_cluster(labs, gene_symbol, variant_id, period_start, initial_window,
                                               _cluster_size(per_grouping, options["somatic_spread"])))

        somatic_groups = [group for group in groups if group.records[0].somatic]
        somatic_records = sum(len(group.records) for group in somatic_groups)
        self.stdout.write(f"{len(groups) - len(somatic_groups)} germline records, {somatic_records} somatic "
                          f"records in {len(somatic_groups)} clusters, over {len(variant_ids_by_gene)} genes")
        return groups

    def _create_allele_infos(self, genome_build: GenomeBuild, groups: list[FakeAlleleGroup]):
        """ The gene chart reads gene_symbol off the allele info, and a record only reaches the
            classification listing once it has an allele to be grouped under """
        patch_version = GenomeBuildPatchVersion.get_or_create(genome_build.name)
        alleles = Allele.objects.bulk_create([Allele() for _ in groups], batch_size=BATCH_SIZE)
        VariantAllele.objects.bulk_create([
            VariantAllele(variant_id=group.variant_id, genome_build=genome_build, allele=allele,
                          origin=AlleleOrigin.IMPORTED_TO_DATABASE,
                          allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
            for group, allele in zip(groups, alleles)], batch_size=BATCH_SIZE)

        allele_infos = ImportedAlleleInfo.objects.bulk_create([
            ImportedAlleleInfo(imported_c_hgvs=f"{group.gene_symbol}:c.{group.variant_id}A>G",
                               imported_genome_build_patch_version=patch_version,
                               allele=allele,
                               status=ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS)
            for group, allele in zip(groups, alleles)], batch_size=BATCH_SIZE)

        variant_infos = ResolvedVariantInfo.objects.bulk_create([
            ResolvedVariantInfo(allele_info=allele_info, genome_build=genome_build,
                                variant_id=group.variant_id, gene_symbol_id=group.gene_symbol,
                                c_hgvs=allele_info.imported_c_hgvs)
            for group, allele_info in zip(groups, allele_infos)], batch_size=BATCH_SIZE)

        build_field = "grch38" if genome_build.name == "GRCh38" else "grch37"
        for group, allele_info, variant_info in zip(groups, allele_infos, variant_infos):
            setattr(allele_info, build_field, variant_info)
            for record in group.records:
                record.allele_info_id = allele_info.pk
        ImportedAlleleInfo.objects.bulk_update(allele_infos, [build_field], batch_size=BATCH_SIZE)
        self.stdout.write(f"Created {len(allele_infos)} alleles and allele infos")

    def _create_classifications(self, records: list[FakeRecord], curators: list[User]):
        classifications = []
        for index, record in enumerate(records):
            latest = record.steps[-1]
            prefix = SOMATIC_LAB_RECORD_PREFIX if record.somatic else LAB_RECORD_PREFIX
            classifications.append(Classification(
                user=random.choice(curators), lab=record.lab_obj,
                lab_record_id=f"{prefix}{index + 1}",
                allele_info_id=record.allele_info_id,
                evidence=latest.evidence,
                clinical_significance=latest.significance,
                allele_origin_bucket=AlleleOriginBucket.SOMATIC if record.somatic else AlleleOriginBucket.GERMLINE,
                share_level=ShareLevel.ALL_USERS.key,
                created=record.started, modified=latest.published_on))

        with _keeping_our_timestamps(Classification):
            Classification.objects.bulk_create(classifications, batch_size=BATCH_SIZE)
        for record, classification in zip(records, classifications):
            record.classification_id = classification.pk
        self.stdout.write(f"Created {len(classifications)} classifications")

    def _create_modifications(self, records: list[FakeRecord],
                              curators: list[User]) -> list[ClassificationModification]:
        modifications = []
        for record in records:
            curator = random.choice(curators)
            for position, step in enumerate(record.steps):
                modifications.append(ClassificationModification(
                    classification_id=record.classification_id, user=curator,
                    source=SubmissionSource.API, published=True, published_evidence=step.evidence,
                    delta=step.evidence if position == 0 else {},
                    clinical_significance=step.significance,
                    share_level=ShareLevel.ALL_USERS.key,
                    is_last_published=position == len(record.steps) - 1,
                    is_last_edited=position == len(record.steps) - 1,
                    created=step.published_on, modified=step.published_on))

        with _keeping_our_timestamps(ClassificationModification):
            ClassificationModification.objects.bulk_create(modifications, batch_size=BATCH_SIZE)
        self.stdout.write(f"Created {len(modifications)} published modifications")
        return modifications

    def _assign_permissions(self, records: list[FakeRecord], modifications: list[ClassificationModification]):
        """ bulk_create skips the mixin's save(), and everyone can see (and clean up) fake data """
        group = all_users_group()
        total = 0
        for klass, objects in ((Classification, [r.classification_id for r in records]),
                               (ClassificationModification, [m.pk for m in modifications])):
            content_type = ContentType.objects.get_for_model(klass)
            permissions = Permission.objects.filter(content_type=content_type,
                                                    codename__in=[klass.get_read_perm(), klass.get_write_perm()])
            for permission in permissions:
                for i in range(0, len(objects), BATCH_SIZE):
                    GroupObjectPermission.objects.bulk_create(
                        [GroupObjectPermission(group=group, permission=permission, content_type=content_type,
                                               object_pk=str(pk)) for pk in objects[i:i + BATCH_SIZE]])
                    total += len(objects[i:i + BATCH_SIZE])
        self.stdout.write(f"Created {total} group permissions")

    def _summarise_and_group(self, records: list[FakeRecord]):
        """ The classification listing renders groupings, and both they and the summary they read are
            derived from the published modifications the same way an import would leave them """
        modifications_qs = ClassificationModification.objects \
            .filter(classification_id__in=[r.classification_id for r in records], is_last_published=True) \
            .select_related("classification", "classification__allele_info")

        classifications = []
        for modification in modifications_qs:
            classification = modification.classification
            classification.summary = ClassificationSummaryCalculator(modification).cache_dict()
            classifications.append(classification)
        Classification.objects.bulk_update(classifications, ["summary"], batch_size=BATCH_SIZE)

        for classification in classifications:
            ClassificationGrouping.assign_grouping_for_classification(classification)
        for grouping in ClassificationGrouping.objects.filter(dirty=True).iterator():
            grouping.update()
        for allele_origin_grouping in AlleleOriginGrouping.objects.filter(dirty=True).iterator():
            allele_origin_grouping.update()

        groupings = ClassificationGroupingEntry.objects.filter(classification__in=classifications) \
            .values("grouping").distinct().count()
        self.stdout.write(f"Summarised {len(classifications)} classifications into {groupings} groupings")

    def _build_timelines(self, records: list[FakeRecord]):
        """ Same call the analytics page makes, so what it renders is what a real import would have left """
        classification_qs = Classification.objects.filter(pk__in=[r.classification_id for r in records])
        events = ReclassificationEventBuilder.rebuild(classification_qs)
        reclassifications = ReclassificationEvent.reclassifications_qs() \
            .filter(classification__in=classification_qs).count()
        self.stdout.write(f"Built {events} timeline events, {reclassifications} of them reclassifications")

    @transaction.atomic
    def delete(self, **_options):
        classifications_qs = Classification.objects.filter(lab_record_id__startswith=LAB_RECORD_PREFIX)
        classification_ids = list(classifications_qs.values_list("pk", flat=True))
        modification_ids = list(ClassificationModification.objects
                                .filter(classification__in=classifications_qs).values_list("pk", flat=True))
        allele_info_ids = [pk for pk in classifications_qs.values_list("allele_info_id", flat=True) if pk]
        allele_ids = [pk for pk in ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids)
                      .values_list("allele_id", flat=True) if pk]

        deleted_permissions = 0
        for klass, object_ids in ((Classification, classification_ids),
                                  (ClassificationModification, modification_ids)):
            content_type = ContentType.objects.get_for_model(klass)
            object_pks = [str(pk) for pk in object_ids]
            for i in range(0, len(object_pks), BATCH_SIZE):
                deleted, _ = GroupObjectPermission.objects.filter(
                    content_type=content_type, object_pk__in=object_pks[i:i + BATCH_SIZE]).delete()
                deleted_permissions += deleted

        ReclassificationEvent.objects.filter(classification__in=classifications_qs).delete()
        classifications_qs.delete()  # cascades to modifications
        ResolvedVariantInfo.objects.filter(allele_info_id__in=allele_info_ids).delete()
        ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids).delete()
        VariantAllele.objects.filter(allele_id__in=allele_ids).delete()
        Allele.objects.filter(pk__in=allele_ids).delete()  # cascades the groupings built on them

        User.objects.filter(username__in=FAKE_CURATORS).delete()
        Lab.objects.filter(group_name__in=[fl.group_name for fl in FAKE_LABS]).delete()
        Organization.objects.filter(group_name=FAKE_ORGANIZATION_GROUP).delete()
        Group.objects.filter(name__startswith=FAKE_ORGANIZATION_GROUP).delete()
        self.stdout.write(f"Deleted {len(classification_ids)} classifications, {len(modification_ids)} "
                          f"modifications, {deleted_permissions} permissions, {len(allele_info_ids)} allele "
                          f"infos, {len(allele_ids)} alleles, and the fake labs/users")


@contextmanager
def _keeping_our_timestamps(klass):
    """ created/modified are auto_now_add/auto_now, which would stamp years of curation as happening now """
    created = klass._meta.get_field("created")
    modified = klass._meta.get_field("modified")
    created.auto_now_add = False
    modified.auto_now = False
    try:
        yield
    finally:
        created.auto_now_add = True
        modified.auto_now = True


def _period_start(years: int) -> datetime:
    return now().astimezone(get_current_timezone()) - timedelta(days=years * 365)


def _started(period_start: datetime, initial_window: float) -> datetime:
    """ Front loaded, so most of the catalogue already exists when the survival cohort is taken """
    return period_start + timedelta(days=random.random() ** 1.6 * initial_window,
                                    hours=random.randrange(8, 18))


def _pick_lab(labs: dict[str, Lab]) -> tuple[FakeCurationLab, Lab]:
    fake_lab = random.choices(FAKE_LABS, weights=[lab.weight for lab in FAKE_LABS])[0]
    return fake_lab, labs[fake_lab.name]


def _germline_record(labs: dict[str, Lab], gene_symbol: str, variant_id: int, started: datetime) -> FakeRecord:
    fake_lab, lab_obj = _pick_lab(labs)
    return FakeRecord(lab=fake_lab, lab_obj=lab_obj, variant_id=variant_id, gene_symbol=gene_symbol,
                      started=started, significance=_weighted_choice(SIGNIFICANCE_WEIGHTS))


def _somatic_cluster(labs: dict[str, Lab], gene_symbol: str, variant_id: int, period_start: datetime,
                     initial_window: float, size: int) -> FakeAlleleGroup:
    """ One lab's cases on the one allele - a grouping is (allele, lab, bucket), so the cluster lands in one """
    fake_lab, lab_obj = _pick_lab(labs)
    return FakeAlleleGroup([
        FakeRecord(lab=fake_lab, lab_obj=lab_obj, variant_id=variant_id, gene_symbol=gene_symbol,
                   started=_started(period_start, initial_window), somatic=True,
                   tier=_weighted_choice(SOMATIC_TIER_WEIGHTS))
        for _ in range(size)])


def _cluster_size(average: int, spread: int) -> int:
    """ How many records share the one grouping - never fewer than one """
    return max(1, random.randint(average - spread, average + spread))


def _weighted_choice(weights: dict[str, float]) -> str:
    values = list(weights)
    return random.choices(values, weights=[weights[value] for value in values])[0]


def _evidence_significance(significance: str) -> str:
    """ Evidence holds the evidence key's own value ('P'), the model fields hold the ClinicalSignificance code """
    return ClinicalSignificance.SHORT_LABELS[significance]


def _clinical_significance_code(evidence_value: str) -> str:
    """ The way back, which is how an oncogenicity call reaches the model fields """
    vg_codes = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("vg")
    return vg_codes[evidence_value]


def _somatic_call(tier: str) -> tuple[str, str]:
    """ The oncogenicity a tier gets called at, as the evidence value and as the code the fields store """
    value = _weighted_choice(ONCOGENICITY_FOR_TIER[tier])
    return value, _clinical_significance_code(value)


def _move_tier(tier: str) -> str:
    """ A review resolves tier_1_or_2, otherwise the tier steps to its neighbour """
    if tier not in SOMATIC_TIER_SCALE:
        return random.choice(["tier_1", "tier_2"])
    index = SOMATIC_TIER_SCALE.index(tier)
    direction = 1 if random.random() < SOMATIC_TOWARDS_ACTIONABLE else -1
    moved = index + direction
    if not 0 <= moved < len(SOMATIC_TIER_SCALE):
        moved = index - direction
    return SOMATIC_TIER_SCALE[moved]


def _build_histories(records: list[FakeRecord], adjacent_share: float):
    """ Walks each record forward from its initial call, deciding at each review whether it moves """
    latest = now().astimezone(get_current_timezone())
    step_weights = _step_size_weights(adjacent_share)

    for record in records:
        if record.somatic:
            _build_somatic_history(record, latest)
        else:
            _build_germline_history(record, latest, step_weights)


def _build_germline_history(record: FakeRecord, latest: datetime, step_weights: dict[int, float]):
    evidence = _initial_evidence(record)
    curated_on = record.started.date()
    record.steps = [FakeStep(record.significance, record.started, curated_on, evidence)]

    num_reviews = _poisson(REVIEWS_PER_RECORD * record.lab.review_rate)
    for published_on in _review_dates(record.started, latest, num_reviews):
        significance = record.steps[-1].significance
        behaviour = SIGNIFICANCE_BEHAVIOUR[significance]
        moves = random.random() < behaviour.move_chance * record.lab.move_rate
        if moves:
            significance = _move(significance, behaviour.benign_bias, step_weights)
        # the curation date is what tells a re-evaluation from an untouched republish, so it always
        # advances even when the call holds
        curated_on = (published_on - timedelta(days=random.randrange(1, 21))).date()
        evidence = _reviewed_evidence(evidence, record.steps[-1].significance, significance, curated_on)
        record.steps.append(FakeStep(significance, published_on, curated_on, evidence))


def _build_somatic_history(record: FakeRecord, latest: datetime):
    """ Somatic curation moves on the tier, and moves less often - it takes new evidence for the tumour
        type rather than another read of the same literature """
    tier = record.tier
    significance_value, significance = _somatic_call(tier)
    curated_on = record.started.date()
    evidence = _initial_somatic_evidence(record, tier, significance_value)
    record.steps = [FakeStep(significance, record.started, curated_on, evidence)]

    num_reviews = _poisson(SOMATIC_REVIEWS_PER_RECORD * record.lab.review_rate)
    for published_on in _review_dates(record.started, latest, num_reviews):
        if random.random() < SOMATIC_MOVE_CHANCE * record.lab.move_rate:
            tier = _move_tier(tier)
            significance_value, significance = _somatic_call(tier)
        curated_on = (published_on - timedelta(days=random.randrange(1, 21))).date()
        evidence = _reviewed_somatic_evidence(evidence, tier, significance_value, curated_on)
        record.steps.append(FakeStep(significance, published_on, curated_on, evidence))


def _review_dates(started: datetime, latest: datetime, num_reviews: int) -> Iterator[datetime]:
    """ Reviews spread over what is left of the record's life, stopping if one would land after today """
    published_on = started
    remaining = (latest - started).days
    for _ in range(num_reviews):
        gap = max(45, int(random.expovariate(1 / (remaining / (num_reviews + 1)))))
        published_on = published_on + timedelta(days=gap, hours=random.randrange(8, 18))
        if published_on >= latest:
            return
        yield published_on


def _step_size_weights(adjacent_share: float) -> dict[int, float]:
    """ Mostly to the neighbouring bucket, and the further it goes the rarer it gets """
    rest = max(0.0, 1 - adjacent_share)
    return {1: adjacent_share, 2: rest * 0.8, 3: rest * 0.15, 4: rest * 0.05}


def _move(significance: str, benign_bias: float, step_weights: dict[int, float]) -> str:
    index = SIGNIFICANCE_SCALE.index(significance)
    direction = -1 if random.random() < benign_bias else 1
    sizes = list(step_weights)
    for size in random.choices(sizes, weights=[step_weights[s] for s in sizes], k=len(sizes)):
        moved = index + direction * size
        if 0 <= moved < len(SIGNIFICANCE_SCALE):
            return SIGNIFICANCE_SCALE[moved]
    # nothing in that direction fits, so it goes the other way instead of standing still
    return SIGNIFICANCE_SCALE[max(0, min(len(SIGNIFICANCE_SCALE) - 1, index - direction))]


def _poisson(mean: float) -> int:
    """ Enough of a Poisson draw for review counts - most records get a couple, a few get many """
    count = 0
    total = random.expovariate(1)
    while total < mean:
        count += 1
        total += random.expovariate(1)
    return count


def _initial_evidence(record: FakeRecord) -> dict:
    evidence = {
        SpecialEKeys.ALLELE_ORIGIN: {"value": "germline"},
        SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": _evidence_significance(record.significance)},
        SpecialEKeys.CURATION_DATE: {"value": record.started.date().isoformat()},
        "gene_symbol": {"value": record.gene_symbol},
        "condition": {"value": f"Fake condition for {record.gene_symbol}"},
        "interpretation_summary": {"value": f"Fake curation of a {record.gene_symbol} variant"},
    }
    for key, strength in _criteria_for(record.significance).items():
        evidence[key] = {"value": strength}
    return evidence


def _criteria_for(significance: str) -> dict[str, str]:
    """ Criteria that would plausibly land on this call - more pathogenic ones the higher up the scale """
    index = SIGNIFICANCE_SCALE.index(significance)
    num_pathogenic = (0, 0, 1, 2, 3)[index]
    num_benign = (3, 2, 1, 0, 0)[index]
    criteria = {}
    for pool, wanted in ((PATHOGENIC_CRITERIA, num_pathogenic), (BENIGN_CRITERIA, num_benign)):
        for key in random.sample(sorted(pool), min(wanted, len(pool))):
            criteria[key] = pool[key]
    if not criteria:
        criteria[random.choice(sorted(PATHOGENIC_CRITERIA))] = CriteriaEvaluation.NOT_MET
    return criteria


def _reviewed_evidence(previous: dict, was: str, now_significance: str, curated_on: date) -> dict:
    """
    The evidence as it stood after a review. Where the call moved, the criteria that change are the ones
    pointing the way it went - that's what the evidence chart splits into applied / unapplied /
    strengthened / weakened.
    """
    evidence = {key: dict(value) for key, value in previous.items()}
    evidence[SpecialEKeys.CLINICAL_SIGNIFICANCE] = {"value": _evidence_significance(now_significance)}
    evidence[SpecialEKeys.CURATION_DATE] = {"value": curated_on.isoformat()}

    delta = SIGNIFICANCE_SCALE.index(now_significance) - SIGNIFICANCE_SCALE.index(was)
    if delta:
        towards_pathogenic = delta > 0
        gaining = PATHOGENIC_CRITERIA if towards_pathogenic else BENIGN_CRITERIA
        losing = BENIGN_CRITERIA if towards_pathogenic else PATHOGENIC_CRITERIA
        for _ in range(min(abs(delta), 2)):
            _apply_criterion(evidence, gaining)
            if random.random() < 0.4:
                _drop_criterion(evidence, losing)

    # a curator rewrites some of the narrative whether or not the call moved
    for key in random.sample(NARRATIVE_KEYS, random.randrange(1, 4)):
        evidence[key] = {"value": f"Fake {key} revised {curated_on.isoformat()}"}
    return evidence


def _initial_somatic_evidence(record: FakeRecord, tier: str, significance_value: str) -> dict:
    evidence = {
        SpecialEKeys.ALLELE_ORIGIN: {"value": "somatic"},
        SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": significance_value},
        SpecialEKeys.CURATION_DATE: {"value": record.started.date().isoformat()},
        "gene_symbol": {"value": record.gene_symbol},
        "condition": {"value": random.choice(FAKE_TUMOUR_TYPES)},
        "somatic:summary_interpretation": {"value": f"Fake somatic curation of a {record.gene_symbol} variant"},
    }
    evidence.update(_tier_evidence(tier))
    return evidence


def _reviewed_somatic_evidence(previous: dict, tier: str, significance_value: str, curated_on: date) -> dict:
    """ The evidence as it stood after a somatic review - the tier and its AMP level move together """
    evidence = {key: dict(value) for key, value in previous.items()}
    evidence[SpecialEKeys.CLINICAL_SIGNIFICANCE] = {"value": significance_value}
    evidence[SpecialEKeys.CURATION_DATE] = {"value": curated_on.isoformat()}
    for amp_level in SpecialEKeys.AMP_LEVELS_TO_LEVEL:
        evidence.pop(amp_level, None)
    evidence.update(_tier_evidence(tier))

    for key in random.sample(SOMATIC_NARRATIVE_KEYS, random.randrange(1, 3)):
        evidence[key] = {"value": f"Fake {key} revised {curated_on.isoformat()}"}
    return evidence


def _tier_evidence(tier: str) -> dict:
    """ The tier and the AMP level it is reported at, which is what the somatic sort is built from """
    evidence = {SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE: {"value": tier}}
    if amp_levels := TIER_AMP_LEVELS.get(tier):
        evidence[random.choice(amp_levels)] = {"value": [random.choice(AMP_CONTEXTS)]}
    return evidence


def _apply_criterion(evidence: dict, pool: dict[str, str]):
    """ Turns a criterion on, or makes one already on count for more """
    already_met = [key for key in pool if evidence.get(key, {}).get("value") in STRONGER_STRENGTH]
    if already_met and random.random() < 0.35:
        key = random.choice(sorted(already_met))
        evidence[key] = {"value": STRONGER_STRENGTH[evidence[key]["value"]]}
        return
    unmet = [key for key in pool if evidence.get(key, {}).get("value") in (None, CriteriaEvaluation.NOT_MET)]
    if unmet:
        key = random.choice(sorted(unmet))
        evidence[key] = {"value": pool[key]}


def _drop_criterion(evidence: dict, pool: dict[str, str]):
    """ Turns a criterion off, or weakens one that no longer carries as much """
    met = [key for key in pool if evidence.get(key, {}).get("value") not in (None, CriteriaEvaluation.NOT_MET)]
    if not met:
        return
    key = random.choice(sorted(met))
    current = evidence[key]["value"]
    if weaker := WEAKER_STRENGTH.get(current):
        evidence[key] = {"value": weaker if random.random() < 0.5 else CriteriaEvaluation.NOT_MET}
    else:
        evidence[key] = {"value": CriteriaEvaluation.NOT_MET}
