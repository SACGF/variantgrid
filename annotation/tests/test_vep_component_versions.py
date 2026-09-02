import importlib

from django.conf import settings
from django.contrib import admin
from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import VariantAnnotationVersion
from annotation.vep_config import parse_cosmic_version_from_filename, vep_component_version_kwargs
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild

# Package-default GRCh38 config, trimmed to the keys this derivation reads. Spelled out rather than read
# from settings.ANNOTATION so the expected versions don't move with a deployment's own annotation config.
GRCH38_SETTINGS = {
    "reference_fasta": "/data/annotation/fasta/GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
    "vep_config": {
        "sift": True,
        "dbscsnv": "annotation_data/GRCh38/dbscSNV1.1_GRCh38.txt.gz",
        "eve": "annotation_data/GRCh38/eve_merged.vcf.gz",
        "gnomad_sv": "annotation_data/GRCh38/gnomad.v4.0.sv.merged.no_filters.vcf.gz",
        "mastermind": "annotation_data/GRCh38/mastermind_cited_variants_reference-2023.10.02-grch38.vcf.gz",
        "maxentscan": "annotation_data/all_builds/maxentscan",
        "phastcons100way": "annotation_data/GRCh38/hg38.phastCons100way.bw",
        "phastcons30way": "annotation_data/GRCh38/hg38.phastCons30way.bw",
        "phastcons46way": None,
        "phylop100way": "annotation_data/GRCh38/hg38.phyloP100way.bw",
        "phylop30way": "annotation_data/GRCh38/hg38.phyloP30way.bw",
        "phylop46way": None,
        "promoter_ai": "annotation_data/GRCh38/promoterAI_tss500.tsv.bgz",
        "protvar": "annotation_data/all_builds/ProtVar_data.db",
        "repeatmasker": "annotation_data/GRCh38/repeatmasker_hg38.bed.gz",
        "topmed": "annotation_data/GRCh38/TOPMED_GRCh38_20180418.vcf.gz",
        "transcript_blocklist": "annotation_data/GRCh38/blocklist_brca1_new_transcripts.txt",
        "uk10k": "annotation_data/GRCh38/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz",
    },
}

GRCH37_SETTINGS = {
    "reference_fasta": "/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz",
    "vep_config": {
        "sift": True,
        "gnomad_sv": "annotation_data/GRCh37/gnomad_v2.1_sv.sites.grch37.converted.no_filters.vcf.gz",
        "topmed": "annotation_data/GRCh37/TOPMED_GRCh37.vcf.gz",
        "eve": None,
        "promoter_ai": None,
        "transcript_blocklist": None,
    },
}


class VEPComponentVersionKwargsTests(TestCase):
    """ #462 - the VEP command line components the ##VEP= header doesn't report. """

    def test_versions_parsed_from_data_files(self):
        kwargs = vep_component_version_kwargs(GRCH38_SETTINGS)
        self.assertEqual("2023.10.02", kwargs["mastermind"])
        self.assertEqual("1.1", kwargs["dbscsnv"])
        self.assertEqual("20180418", kwargs["topmed"])
        self.assertEqual("20160215", kwargs["uk10k"])
        self.assertEqual("4.0", kwargs["gnomad_sv"])

    def test_gnomad_sv_grch37_naming(self):
        kwargs = vep_component_version_kwargs(GRCH37_SETTINGS)
        self.assertEqual("2.1", kwargs["gnomad_sv"])

    def test_unparseable_version_falls_back_to_basename(self):
        kwargs = vep_component_version_kwargs(GRCH37_SETTINGS)
        self.assertEqual("TOPMED_GRCh37.vcf.gz", kwargs["topmed"])

    def test_unversioned_files_pin_basename(self):
        kwargs = vep_component_version_kwargs(GRCH38_SETTINGS)
        self.assertEqual("repeatmasker_hg38.bed.gz", kwargs["repeat_masker"])
        self.assertEqual("blocklist_brca1_new_transcripts.txt", kwargs["transcript_blocklist"])
        self.assertEqual("ProtVar_data.db", kwargs["protvar"])
        self.assertEqual("maxentscan", kwargs["maxentscan"])

    def test_unconfigured_data_is_none(self):
        kwargs = vep_component_version_kwargs(GRCH37_SETTINGS)
        self.assertIsNone(kwargs["eve"])
        self.assertIsNone(kwargs["promoter_ai"])
        self.assertIsNone(kwargs["transcript_blocklist"])

    def test_conservation_pins_configured_tracks_only(self):
        kwargs = vep_component_version_kwargs(GRCH38_SETTINGS)
        self.assertEqual("hg38.phastCons100way.bw,hg38.phastCons30way.bw,"
                         "hg38.phyloP100way.bw,hg38.phyloP30way.bw", kwargs["conservation"])

    def test_no_filesystem_access(self):
        """ The migration backfill runs wherever the upgrade runs, which may not have the data mounted """
        build_settings = dict(GRCH38_SETTINGS)
        build_settings["vep_config"] = dict(GRCH38_SETTINGS["vep_config"])
        build_settings["vep_config"]["mastermind"] = \
            "annotation_data/nowhere/mastermind_cited_variants_reference-2029.01.01-grch38.vcf.gz"
        kwargs = vep_component_version_kwargs(build_settings)
        self.assertEqual("2029.01.01", kwargs["mastermind"])

    def test_fasta_pins_reference_fasta(self):
        kwargs = vep_component_version_kwargs(GRCH38_SETTINGS)
        self.assertEqual("GCF_000001405.39_GRCh38.p13_genomic.fna.gz", kwargs["fasta"])

    def test_fasta_pins_pre_vep112_override(self):
        """ use_pre_vep112_fasta() sets vep_config "fasta", which is what reaches the command line """
        build_settings = dict(GRCH38_SETTINGS)
        build_settings["vep_config"] = dict(GRCH38_SETTINGS["vep_config"])
        build_settings["vep_config"]["fasta"] = \
            "/data/annotation/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
        kwargs = vep_component_version_kwargs(build_settings)
        self.assertEqual("Homo_sapiens.GRCh38.dna.toplevel.fa.gz", kwargs["fasta"])

    def test_sift_flag_pinned(self):
        self.assertTrue(vep_component_version_kwargs(GRCH38_SETTINGS)["sift_enabled"])
        self.assertFalse(vep_component_version_kwargs({"vep_config": {"sift": False}})["sift_enabled"])

    @override_settings(ANNOTATION_VEP_PICK_ORDER="mane_select,canonical",
                       ANNOTATION_VEP_ARGS=["--mane", "--shift_hgvs", "0"],
                       ANNOTATION_VEP_SV_MAX_SIZE=5_000_000,
                       ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION=0.5)
    def test_command_line_settings_pinned(self):
        kwargs = vep_component_version_kwargs(GRCH38_SETTINGS)
        self.assertEqual("mane_select,canonical", kwargs["pick_order"])
        self.assertEqual("--mane --shift_hgvs 0", kwargs["vep_args"])
        self.assertEqual(5_000_000, kwargs["sv_max_size"])
        self.assertEqual(0.5, kwargs["sv_overlap_min_fraction"])

    @override_settings(ANNOTATION_VEP_PICK_ORDER=None, ANNOTATION_VEP_ARGS=[])
    def test_empty_command_line_settings_are_none(self):
        kwargs = vep_component_version_kwargs(GRCH38_SETTINGS)
        self.assertIsNone(kwargs["vep_args"])
        self.assertIsNone(kwargs["pick_order"])


class ParseCosmicVersionTests(TestCase):
    """ #1673 - the release drives which INFO field carries the sample count, and it's only ever
        written down in the distributed VCF's name. """

    def test_release_parsed_from_distributed_names(self):
        expected = {
            "annotation_data/GRCh37/CosmicCodingMuts_v95_20211101_grch37.normal.vcf.gz": 95,
            "annotation_data/GRCh38/Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz": 99,
            "annotation_data/GRCh37/Cosmic_GenomeScreensMutant_Normal_v101_GRCh37.vcf.gz": 101,
        }
        for path, cosmic_version in expected.items():
            with self.subTest(path=path):
                self.assertEqual(parse_cosmic_version_from_filename(path), cosmic_version)

    def test_unversioned_name_is_none(self):
        self.assertIsNone(
            parse_cosmic_version_from_filename("annotation_data/GRCh37/CosmicCodingMuts.normal.grch37.vcf.gz"))


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class VEPComponentVersionBackfillTests(TestCase):
    """ #462 - migration 0161 stamps the versions the scheduler can still match, so a deployment whose
        annotation config hasn't changed keeps matching its existing VariantAnnotationVersion rather than
        creating a new one and reannotating. """

    MIGRATION = "annotation.migrations.0161_vep_command_line_component_versions"
    PIN_FIELDS = ("mastermind", "dbscsnv", "topmed", "uk10k", "gnomad_sv", "conservation",
                  "repeat_masker", "maxentscan", "protvar", "eve", "promoter_ai", "fasta",
                  "transcript_blocklist", "pick_order", "sift_enabled", "sv_max_size",
                  "sv_overlap_min_fraction", "vep_args")

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _make_unpinned_vav(self, status) -> VariantAnnotationVersion:
        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = status
        kwargs.update(dict.fromkeys(self.PIN_FIELDS))
        return VariantAnnotationVersion.objects.create(**kwargs)

    def _run_backfill(self):
        migration = importlib.import_module(self.MIGRATION)

        class _Apps:
            @staticmethod
            def get_model(_app_label, _model_name):
                return VariantAnnotationVersion

        migration._backfill_vep_components(_Apps(), None)

    def test_active_version_stamped_with_current_settings(self):
        vav = self._make_unpinned_vav(VariantAnnotationVersion.Status.ACTIVE)
        self._run_backfill()
        vav.refresh_from_db()

        expected = vep_component_version_kwargs(settings.ANNOTATION[self.grch37.name])
        for field, value in expected.items():
            self.assertEqual(value, getattr(vav, field), field)

    def test_new_version_stamped(self):
        vav = self._make_unpinned_vav(VariantAnnotationVersion.Status.NEW)
        self._run_backfill()
        vav.refresh_from_db()
        self.assertEqual(settings.ANNOTATION_VEP_SV_MAX_SIZE, vav.sv_max_size)

    def test_historical_version_left_unknown(self):
        vav = self._make_unpinned_vav(VariantAnnotationVersion.Status.HISTORICAL)
        self._run_backfill()
        vav.refresh_from_db()
        for field in self.PIN_FIELDS:
            self.assertIsNone(getattr(vav, field), field)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AdminBlankPinTests(TestCase):
    """ Blanking a nullable TextField in the admin gives '' - which isn't the None the pins are
        derived as for "not configured", so the version would stop matching current VEP """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.admin_user = User.objects.create_superuser(username="admin", email="a@example.com", password="x")

    def test_blanked_pin_still_matches_current_vep(self):
        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs.pop("id")
        kwargs["pick_order"] = None
        vav = VariantAnnotationVersion.objects.create(**kwargs)

        model_admin = admin.site._registry[VariantAnnotationVersion]
        request = RequestFactory().get("/")
        request.user = self.admin_user
        pick_order_field = model_admin.get_form(request, vav, change=True).base_fields["pick_order"]

        vav.pick_order = pick_order_field.clean("")  # blanked out by an admin edit
        vav.save()
        vav.refresh_from_db()
        self.assertIsNone(vav.pick_order)
        _, created = VariantAnnotationVersion.objects.get_or_create(**kwargs)
        self.assertFalse(created, "Edited version still matches current VEP, rather than reannotating")
