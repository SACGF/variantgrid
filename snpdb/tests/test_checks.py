"""
The settings system checks (snpdb/checks.py). Each check is called directly under override_settings; the
suite's own settings must pass them all, which Django verifies before any test runs.
"""
import runpy
from copy import deepcopy

from django.conf import settings
from django.test import SimpleTestCase, override_settings

from snpdb import checks
from variantgrid.settings.components import annotation_settings
from variantgrid.settings.components.annotation_settings import (
    BUILD_GRCH37,
    BUILD_GRCH38,
    pin_annotation_to_columns_version_3,
    use_pre_vep112_fasta,
)

OIDC_BACKENDS = (*settings.AUTHENTICATION_BACKENDS[:1], checks.OIDC_BACKEND, *settings.AUTHENTICATION_BACKENDS[1:])
OIDC_CONFIGURED = {
    "USE_OIDC": True, "LOGIN_URL": checks.OIDC_LOGIN_URL, "AUTHENTICATION_BACKENDS": OIDC_BACKENDS,
    "OIDC_RP_CLIENT_ID": "shariant-test", "OIDC_RP_CLIENT_SECRET": "secret",
    "OIDC_REQUIRED_GROUP": "/variantgrid/shariant_test",
    "LAB_HEAD_MANAGE_MEMBERS": False, "USER_CREATE_ORG_LABS": {}, "REGISTRATION_OPEN": False,
}


# Re-run rather than imported: an env settings file edits the ANNOTATION dict it star-imports in place
PACKAGE_ANNOTATION_SETTINGS = runpy.run_path(annotation_settings.__file__)
PACKAGE_VEP_SETTINGS = {name: PACKAGE_ANNOTATION_SETTINGS[name] for name in (
    "ANNOTATION_VEP_VERSION", "ANNOTATION_VEP_ENSEMBL_GENCODE", "ANNOTATION_VEP_FORK", "ANNOTATION_VEP_BUFFER_SIZE",
    "ANNOTATION_VEP_BATCH_MIN", "ANNOTATION_VEP_BATCH_MAX", "ANNOTATION_MAX_BENIGN_RANKSCORE",
    "ANNOTATION_MIN_PATHOGENIC_RANKSCORE")}


def ids(messages) -> set[str]:
    return {m.id for m in messages}


class OidcCheckTest(SimpleTestCase):

    def test_local_auth_defaults_pass(self):
        self.assertEqual(ids(checks.check_oidc(None)), set())

    @override_settings(AUTHENTICATION_BACKENDS=OIDC_BACKENDS)
    def test_backend_without_use_oidc(self):
        self.assertEqual(ids(checks.check_oidc(None)), {"variantgrid.E005"})

    @override_settings(**OIDC_CONFIGURED)
    def test_configured_oidc_only_lacks_the_app(self):
        """ 'oidc_auth' can't be added to INSTALLED_APPS in a test, so a fully configured OIDC reports just that """
        messages = checks.check_oidc(None)
        self.assertEqual(ids(messages), {"variantgrid.E005"})
        self.assertIn("'oidc_auth' in INSTALLED_APPS", messages[0].msg)

    @override_settings(**{**OIDC_CONFIGURED, "LAB_HEAD_MANAGE_MEMBERS": True,
                          "USER_CREATE_ORG_LABS": {"unknown": "unknown"}, "REGISTRATION_OPEN": True,
                          "OIDC_REQUIRED_GROUP": None})
    def test_local_account_features_under_oidc(self):
        self.assertEqual(ids(checks.check_oidc(None)),
                         {"variantgrid.E002", "variantgrid.E003", "variantgrid.E004", "variantgrid.E005",
                          "variantgrid.E006"})


class FeatureDependencyCheckTest(SimpleTestCase):
    """ Checked against a baseline with every feature off, not the settings of whichever host runs the suite """
    BASELINE = {
        "MME_ENABLED": False, "MME_CONTACT": None,
        "BEACON_ENABLED": False, "BEACON_CONFIG": {**settings.BEACON_CONFIG, "beacon_id": checks.BEACON_PLACEHOLDER_ID},
        "DISCORDANCE_ENABLED": False, "DISCORDANCE_EMAIL": None,
        "CONTACT_US_ENABLED": False, "KEYCLOAK_SYNC_DETAILS": None, "REGISTRATION_OPEN": False,
    }

    def test_defaults_pass(self):
        with override_settings(**self.BASELINE):
            self.assertEqual(ids(checks.check_feature_dependencies(None)), set())

    def test_feature_without_its_dependency(self):
        cases = {
            "variantgrid.E007": {"MME_ENABLED": True, "MME_FROM_EMAIL": None, "SEND_EMAILS": True,
                                 "MME_CONTACT": {"name": "Lab", "href": "mailto:lab@example.org"}},
            "variantgrid.E008": {"BEACON_ENABLED": True},
            "variantgrid.E009": {"DISCORDANCE_ENABLED": True, "DISCORDANCE_EMAIL": None},
            "variantgrid.W001": {"DISCORDANCE_EMAIL": "updates@example.org", "SEND_EMAILS": False},
            "variantgrid.E010": {"CONTACT_US_ENABLED": True, "ADMIN_EMAIL_NOTIFICATION": None,
                                 "SLACK": {"enabled": False}},
            "variantgrid.E011": {"REGISTRATION_OPEN": True, "RECAPTCHA_PUBLIC_KEY": "", "RECAPTCHA_PRIVATE_KEY": "x"},
        }
        for expected, overrides in cases.items():
            with self.subTest(expected), override_settings(**{**self.BASELINE, **overrides}):
                self.assertEqual(ids(checks.check_feature_dependencies(None)), {expected})

    def test_mme_needs_a_contact(self):
        with override_settings(**{**self.BASELINE, "MME_ENABLED": True, "MME_FROM_EMAIL": "mme@example.org",
                                  "SEND_EMAILS": True}):
            self.assertEqual(ids(checks.check_feature_dependencies(None)), {"variantgrid.E007"})

    def test_keycloak_admin_needs_an_accounts_email(self):
        with override_settings(**{**self.BASELINE, "KEYCLOAK_SYNC_DETAILS": {"host": "https://auth.example.org"},
                                  "ACCOUNTS_EMAIL": None}):
            self.assertEqual(ids(checks.check_feature_dependencies(None)), {"variantgrid.E010"})


@override_settings(**PACKAGE_VEP_SETTINGS)
class AnnotationCheckTest(SimpleTestCase):
    """ Checked against the package default (latest) annotation settings, not the host's env file """

    @staticmethod
    def _annotation(**build_overrides):
        annotation = deepcopy(PACKAGE_ANNOTATION_SETTINGS["ANNOTATION"])
        annotation[BUILD_GRCH37]["enabled"] = True
        annotation[BUILD_GRCH38]["enabled"] = True
        for build, overrides in build_overrides.items():
            annotation[build].update(overrides)
        return annotation

    def test_latest_defaults_pass(self):
        with override_settings(ANNOTATION=self._annotation()):
            self.assertEqual(ids(checks.check_annotation(None)), set())

    def test_old_vep_against_latest_annotation(self):
        """ The historical VEP with the package-default (latest) annotation: no chromosome-named fasta, the
            columns_version 5 plugin columns, and GRCh38 configures the VEP 116 plugins """
        with override_settings(ANNOTATION_VEP_VERSION="110", ANNOTATION=self._annotation()):
            self.assertEqual(ids(checks.check_annotation(None)),
                             {"variantgrid.E012", "variantgrid.E013", "variantgrid.W002"})

    def test_pinned_historical_annotation_passes(self):
        """ The env-file recipe for staying on VEP 110 satisfies the check """
        annotation = self._annotation()
        pin_annotation_to_columns_version_3(annotation)
        use_pre_vep112_fasta(annotation)
        with override_settings(ANNOTATION_VEP_VERSION="110", ANNOTATION=annotation):
            self.assertEqual(ids(checks.check_annotation(None)), set())

    def test_invalid_values(self):
        no_builds = self._annotation(**{BUILD_GRCH37: {"enabled": False}, BUILD_GRCH38: {"enabled": False}})
        no_builds["T2T-CHM13v2.0"]["enabled"] = False
        cases = [
            ("no build enabled", {"ANNOTATION": no_builds}, "variantgrid.E014"),
            ("consortium", {"ANNOTATION": self._annotation(**{BUILD_GRCH37: {"annotation_consortium": "refseq"}})},
             "variantgrid.E014"),
            ("gencode", {"ANNOTATION_VEP_ENSEMBL_GENCODE": "all"}, "variantgrid.E014"),
            ("buffer vs fork", {"ANNOTATION_VEP_FORK": 50}, "variantgrid.E015"),
            ("batch range", {"ANNOTATION_VEP_BATCH_MIN": 30_000}, "variantgrid.E016"),
            ("rankscores", {"ANNOTATION_MAX_BENIGN_RANKSCORE": 0.9}, "variantgrid.E017"),
        ]
        for name, overrides, expected in cases:
            overrides.setdefault("ANNOTATION", self._annotation())
            with self.subTest(name), override_settings(**overrides):
                self.assertEqual(ids(checks.check_annotation(None)), {expected})


class OrderedThresholdCheckTest(SimpleTestCase):

    def test_defaults_pass(self):
        self.assertEqual(ids(checks.check_ordered_thresholds(None)), set())

    def test_out_of_order(self):
        cases = [
            ("heartbeat vs lease", {"ANNOTATION_RUN_LEASE_HEARTBEAT_SECONDS": 900}, "variantgrid.E018"),
            ("upload heartbeat", {"UPLOAD_STEP_LEASE_HEARTBEAT_SECONDS": 600}, "variantgrid.E018"),
            ("nearby timeout", {"VARIANT_DETAILS_NEARBY_STATEMENT_TIMEOUT_SECONDS": 601}, "variantgrid.E019"),
            ("slot expiry", {"MAJOR_OPERATION_SLOT_EXPIRE_SECONDS": 60}, "variantgrid.E019"),
            ("ref mismatch", {"VCF_IMPORT_REF_MISMATCH_WARN_FRACTION": 0.6}, "variantgrid.E020"),
            ("genotype quality", {"VARIANT_GRID_GENOTYPE_QUALITY_THRESHOLDS": {"gq": {"bad": 30, "marginal": 20},
                                                                               "pl": {"bad": 50, "marginal": 20}}},
             "variantgrid.E021"),
            ("sv min size", {"ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE": 10}, "variantgrid.W003"),
        ]
        for name, overrides, expected in cases:
            with self.subTest(name), override_settings(**overrides):
                self.assertEqual(ids(checks.check_ordered_thresholds(None)), {expected})


class EnumerationCheckTest(SimpleTestCase):

    def test_defaults_pass(self):
        self.assertEqual(ids(checks.check_enumerations(None)), set())

    def test_value_outside_the_set(self):
        cases = {
            "bucket": {"ALLELE_ORIGIN_NOT_PROVIDED_BUCKET": "X"},
            "sample node levels": {"ANALYSIS_SAMPLE_NODE_LEVELS": ["S", "Q"]},
            "granularity order": {"BEACON_CONFIG": {**settings.BEACON_CONFIG, "default_granularity": "record",
                                                    "max_granularity": "boolean"}},
            "publish level": {"CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL": "everyone"},
        }
        for name, overrides in cases.items():
            with self.subTest(name), override_settings(**overrides):
                self.assertEqual(ids(checks.check_enumerations(None)), {"variantgrid.E022"})


class ReferenceCheckTest(SimpleTestCase):

    def test_defaults_pass(self):
        self.assertEqual(ids(checks.check_dotted_paths(None)), set())
        self.assertEqual(ids(checks.check_url_register(None)), set())
        self.assertEqual(ids(checks.check_gene_grid_columns(None)), set())
        self.assertEqual(ids(checks.check_public_sample_gene_matrix(None)), set())

    @override_settings(MENUS="variantgrid.menus.NO_SUCH_MENUS",
                       FINISH_IMPORT_VCF_STEP_TASKS_CLASSES=["no_such_app.tasks.Task"])
    def test_dotted_path_typo(self):
        messages = checks.check_dotted_paths(None)
        self.assertEqual(ids(messages), {"variantgrid.E023"})
        self.assertEqual(len(messages), 2)

    def test_url_register_typo(self):
        register = {**settings.URLS_NAME_REGISTER, "view_patinet": False}
        with override_settings(URLS_NAME_REGISTER=register, URLS_APP_REGISTER={"analyses": False}):
            messages = checks.check_url_register(None)
        self.assertEqual(ids(messages), {"variantgrid.W004"})
        self.assertIn("view_patinet", messages[0].msg)
        self.assertIn("analyses", messages[1].msg)

    @override_settings(GENE_GRID_ENRICHMENT_KIT_COLUMN_LABELS=["Mean"])
    def test_gene_grid_lists_must_line_up(self):
        self.assertEqual(ids(checks.check_gene_grid_columns(None)), {"variantgrid.E024"})

    def test_public_matrix_needs_a_build(self):
        register = {**settings.URLS_NAME_REGISTER, "public_global_sample_gene_matrix": True}
        annotation = deepcopy(settings.ANNOTATION)
        annotation[BUILD_GRCH37]["enabled"] = annotation[BUILD_GRCH38]["enabled"] = True
        with override_settings(URLS_NAME_REGISTER=register, ANNOTATION=annotation,
                               PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD=None):
            self.assertEqual(ids(checks.check_public_sample_gene_matrix(None)), {"variantgrid.E025"})
        with override_settings(URLS_NAME_REGISTER=register, ANNOTATION=annotation,
                               PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD=BUILD_GRCH37):
            self.assertEqual(ids(checks.check_public_sample_gene_matrix(None)), set())


class ServedConfigurationCheckTest(SimpleTestCase):
    """ Checked against a known-good baseline, not the settings of whichever host runs the suite """
    BASELINE = {
        "DEBUG": False, "UNIT_TEST": True,
        "SECURE_PROXY_SSL_HEADER": None,
        "REST_FRAMEWORK": {**settings.REST_FRAMEWORK, "DEFAULT_SCHEMA_CLASS": checks.SPECTACULAR_SCHEMA_CLASS},
        "MIDDLEWARE": (*(m for m in settings.MIDDLEWARE if m != checks.AXES_MIDDLEWARE), checks.AXES_MIDDLEWARE),
    }

    def test_defaults_pass(self):
        with override_settings(**self.BASELINE):
            self.assertEqual(ids(checks.check_served_configuration(None)), set())

    def test_misconfigurations(self):
        rest_framework_without_schema = {k: v for k, v in self.BASELINE["REST_FRAMEWORK"].items()
                                         if k != "DEFAULT_SCHEMA_CLASS"}
        cases = [
            ("dev switch on a server", {"UNIT_TEST": False, "CACHE_GENERATED_FILES": False}, "variantgrid.W005"),
            ("https without secure cookies", {"SECURE_PROXY_SSL_HEADER": ('HTTP_X_FORWARDED_PROTOCOL', 'https'),
                                              "SESSION_COOKIE_SECURE": False, "CSRF_COOKIE_SECURE": True},
             "variantgrid.W006"),
            ("schema class dropped", {"REST_FRAMEWORK": rest_framework_without_schema}, "variantgrid.E026"),
            ("axes not last", {"MIDDLEWARE": (*self.BASELINE["MIDDLEWARE"], 'eventlog.middleware.PageViewsMiddleware')},
             "variantgrid.E027"),
        ]
        for name, overrides, expected in cases:
            with self.subTest(name), override_settings(**{**self.BASELINE, **overrides}):
                self.assertEqual(ids(checks.check_served_configuration(None)), {expected})
