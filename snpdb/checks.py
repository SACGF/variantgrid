"""
System checks for VariantGrid's own settings: features that contradict each other, thresholds that must be
ordered, and names (dotted paths, URL names) that must resolve. Django runs them on every manage.py command
and at server start; `manage.py deployment_check` reports them alongside its own checks
(`variantgrid/deployment_validation/settings_check.py`), so the upgrader stops on an Error.

A setting the code already validates at the point of use isn't repeated here. Ids: variantgrid.E0xx stop the
command, variantgrid.W0xx are printed. Silence one a deployment does on purpose with SILENCED_SYSTEM_CHECKS.
"""
import sys
from importlib import import_module

from django.conf import settings
from django.core.checks import Error, Tags, Warning, register
from django.urls import get_resolver
from django.utils.module_loading import import_string

from classification.enums.classification_enums import ShareLevel
from genes.enums import AnnotationConsortium
from patients.models_enums import SampleSourceLevel
from snpdb.models.models_enums import AlleleOriginBucket

OIDC_BACKEND = 'oidc_auth.backend.VariantGridOIDCAuthenticationBackend'
OIDC_LOGIN_URL = '/oidc_login/'
AXES_MIDDLEWARE = 'axes.middleware.AxesMiddleware'
SPECTACULAR_SCHEMA_CLASS = 'drf_spectacular.openapi.AutoSchema'
BEACON_PLACEHOLDER_ID = "org.variantgrid.beacon"
VEP_FASTA_RENAME_FIXED_VERSION = 112  # below this VEP renamed contigs against an NCBI fasta, dropping plugin output
VEP_PLUGIN_COLUMNS_VERSION = 5  # columns_version that added the VEP 116 plugins (#1638)
VEP_PLUGIN_MIN_VERSION = 116
VEP_116_PLUGINS = ("eve", "popeve", "promoter_ai")
GENCODE_VALUES = ("primary", "basic", None)
SV_OVERLAP_METHODS = ("greatest_overlap", "lowest_af", "exact_or_lowest_af")
MME_GENE_ID_PREFERENCES = ("ensembl", "symbol")
BEACON_ENVIRONMENTS = ("prod", "test", "dev", "staging")
BEACON_GRANULARITIES = ("boolean", "count", "record")
# Setting -> the dotted paths it holds. Each must import: a typo is otherwise an ImportError on the first use
DOTTED_PATH_SETTINGS = {
    "MENUS": "string",
    "CLINGEN_ALLELE_REGISTRY_API_CLASS": "string",
    "TRANSCRIPT_SEQUENCE_FETCHER_CLASS": "string",
    "PATIENT_SPECIMEN_REFERENCE_NORMALIZER": "string",
    "FINISH_IMPORT_VCF_STEP_TASKS_CLASSES": "list",
    "IMPORT_TASK_FACTORY_IMPORTS": "modules",
}


@register(Tags.security)
def check_secure_cookies_not_used_with_runserver(app_configs, **kwargs):
    """runserver only speaks HTTP, and browsers never store a Secure cookie from an HTTP response, so
    logging in silently fails - you get sent back to the login page with no error."""
    errors = []
    if "runserver" in sys.argv:
        secure_cookie_settings = [name for name in ("SESSION_COOKIE_SECURE", "CSRF_COOKIE_SECURE")
                                  if getattr(settings, name)]
        if secure_cookie_settings:
            errors.append(
                Error(
                    f"{', '.join(secure_cookie_settings)} set while running runserver - nobody will be able to log in",
                    hint="Development machines copy variantgrid/settings/env_developers/_settings_template.py, "
                         "which leaves out the https_settings import. If you are serving HTTPS locally, silence "
                         "this with SILENCED_SYSTEM_CHECKS.",
                    id="variantgrid.E001",
                )
            )
    return errors


@register(Tags.security)
def check_oidc(app_configs, **kwargs):
    """Under OIDC, Keycloak owns the account and rewrites the user's groups from its claims on every login
    (oidc_auth/backend.py), so anything that creates accounts or assigns groups locally is undone."""
    messages = []
    backend_installed = OIDC_BACKEND in settings.AUTHENTICATION_BACKENDS
    if not settings.USE_OIDC:
        if backend_installed:
            messages.append(Error(f"{OIDC_BACKEND} is in AUTHENTICATION_BACKENDS but USE_OIDC is False",
                                  hint="Set USE_OIDC = True, or remove the backend.", id="variantgrid.E005"))
        return messages

    if settings.LAB_HEAD_MANAGE_MEMBERS:
        messages.append(Error("LAB_HEAD_MANAGE_MEMBERS is True under OIDC - a lab head's change to membership is "
                              "undone when the user next logs in (Keycloak groups replace the user's groups)",
                              hint="Set LAB_HEAD_MANAGE_MEMBERS = False, as shariantcommon.py does.",
                              id="variantgrid.E002"))
    if settings.USER_CREATE_ORG_LABS:
        messages.append(Error("USER_CREATE_ORG_LABS is set under OIDC - the lab group it adds at user creation is "
                              "removed on the user's next login",
                              hint="Set USER_CREATE_ORG_LABS = {} and put users in labs through Keycloak.",
                              id="variantgrid.E003"))
    if getattr(settings, "REGISTRATION_OPEN", False):
        messages.append(Error("REGISTRATION_OPEN is True under OIDC - accounts live in Keycloak, a locally registered "
                              "user can never log in",
                              hint="Remove REGISTRATION_OPEN.", id="variantgrid.E004"))

    missing = []
    if 'oidc_auth' not in settings.INSTALLED_APPS:
        missing.append("'oidc_auth' in INSTALLED_APPS")
    if not backend_installed:
        missing.append(f"{OIDC_BACKEND} in AUTHENTICATION_BACKENDS")
    if settings.LOGIN_URL != OIDC_LOGIN_URL:
        missing.append(f"LOGIN_URL = '{OIDC_LOGIN_URL}'")
    for name in ("OIDC_RP_CLIENT_ID", "OIDC_RP_CLIENT_SECRET"):
        if not getattr(settings, name, None):
            missing.append(name)
    if missing:
        messages.append(Error(f"USE_OIDC is True but OIDC is only partly configured, missing: {', '.join(missing)}",
                              hint="See the Keycloak section of variantgrid/settings/env/shariantcommon.py.",
                              id="variantgrid.E005"))
    if not settings.OIDC_REQUIRED_GROUP:
        messages.append(Error("USE_OIDC is True without OIDC_REQUIRED_GROUP - any user in the Keycloak realm can log "
                              "in to this environment",
                              hint="Set OIDC_REQUIRED_GROUP to this environment's group, eg '/variantgrid/shariant_test'.",
                              id="variantgrid.E006"))
    return messages


@register(Tags.compatibility)
def check_feature_dependencies(app_configs, **kwargs):
    """A feature switched on without the setting it sends through or identifies itself with."""
    messages = []
    if settings.MME_ENABLED:
        # EmailLog.send_mail is quietly a no-op when from_email is falsy or SEND_EMAILS is False, and MME
        # obliges us to notify depositors of matches
        if not (settings.MME_FROM_EMAIL and settings.SEND_EMAILS):
            messages.append(Error("MME_ENABLED requires MME_FROM_EMAIL and SEND_EMAILS - MatchMaker Exchange must be "
                                  "able to notify depositors of matches", id="variantgrid.E007"))
        contact = settings.MME_CONTACT or {}
        if not (contact.get("name") and contact.get("href")):
            messages.append(Error("MME_ENABLED requires MME_CONTACT with at least 'name' and 'href' - it is sent as the "
                                  "contact block on every outbound patient",
                                  hint="See MME_CONTACT in variantgrid/settings/env/vgaws.py.", id="variantgrid.E007"))
    if settings.BEACON_ENABLED:
        config = settings.BEACON_CONFIG
        if config.get("beacon_id") == BEACON_PLACEHOLDER_ID or not config.get("organization", {}).get("name"):
            messages.append(Error("BEACON_ENABLED with the placeholder BEACON_CONFIG - beacon_id and organization name "
                                  "identify this deployment to the federation",
                                  hint="Override BEACON_CONFIG in the env settings file as vgaws.py does.",
                                  id="variantgrid.E008"))
    if settings.DISCORDANCE_ENABLED and not settings.DISCORDANCE_EMAIL:
        messages.append(Warning("DISCORDANCE_ENABLED without DISCORDANCE_EMAIL - discordance notifications are never sent",
                              id="variantgrid.E009"))
    if settings.DISCORDANCE_EMAIL and not settings.SEND_EMAILS:
        messages.append(Warning("DISCORDANCE_EMAIL is set but SEND_EMAILS is False - discordance emails are logged, "
                                "not sent", hint="Expected on a test deployment; silence with SILENCED_SYSTEM_CHECKS.",
                                id="variantgrid.W001"))
    if settings.CONTACT_US_ENABLED and not (settings.ADMIN_EMAIL_NOTIFICATION or settings.SLACK.get("enabled")):
        messages.append(Error("CONTACT_US_ENABLED but neither ADMIN_EMAIL_NOTIFICATION nor SLACK is set - a contact form "
                              "post reaches nobody (AdminNotificationBuilder)", id="variantgrid.E010"))
    if settings.KEYCLOAK_SYNC_DETAILS and not settings.ACCOUNTS_EMAIL:
        messages.append(Error("KEYCLOAK_SYNC_DETAILS without ACCOUNTS_EMAIL - the Keycloak admin page sends the welcome "
                              "email from it", id="variantgrid.E010"))
    if getattr(settings, "REGISTRATION_OPEN", False) and not (settings.RECAPTCHA_PUBLIC_KEY and
                                                               settings.RECAPTCHA_PRIVATE_KEY):
        messages.append(Error("REGISTRATION_OPEN without RECAPTCHA_PUBLIC_KEY and RECAPTCHA_PRIVATE_KEY - the "
                              "registration form is the reCAPTCHA one",
                              hint="RECAPTCHA_PUBLIC_KEY = get_secret('RECAPTCHA.public_key') etc, as vgaws.py does.",
                              id="variantgrid.E011"))
    return messages


def _enabled_builds() -> dict:
    return {name: config for name, config in settings.ANNOTATION.items() if config.get("enabled")}


@register(Tags.compatibility)
def check_annotation(app_configs, **kwargs):
    """ANNOTATION / VEP settings that have to agree with the installed VEP version and each other."""
    messages = []
    vep_version = int(settings.ANNOTATION_VEP_VERSION)
    builds = _enabled_builds()
    if not builds:
        messages.append(Error("No genome build is enabled in ANNOTATION", id="variantgrid.E014"))
    consortia = {c.label for c in AnnotationConsortium}
    for build_name, config in builds.items():
        vep_config = config["vep_config"]
        if vep_version < VEP_FASTA_RENAME_FIXED_VERSION and not vep_config.get("fasta"):
            messages.append(Error(f"ANNOTATION_VEP_VERSION {vep_version} needs a chromosome-named fasta for {build_name} "
                                  f"- below {VEP_FASTA_RENAME_FIXED_VERSION} VEP renames contigs against an NCBI fasta "
                                  f"and silently drops plugin annotations",
                                  hint="Call use_pre_vep112_fasta(ANNOTATION) in the env settings file.",
                                  id="variantgrid.E012"))
        if vep_version < VEP_PLUGIN_MIN_VERSION:
            if configured := [p for p in VEP_116_PLUGINS if vep_config.get(p)]:
                messages.append(Warning(f"{build_name} configures {', '.join(configured)} but they need VEP "
                                        f">= {VEP_PLUGIN_MIN_VERSION} (ANNOTATION_VEP_VERSION is {vep_version}) so are "
                                        f"skipped", id="variantgrid.W002"))
            if config["columns_version"] >= VEP_PLUGIN_COLUMNS_VERSION:
                messages.append(Error(f"{build_name} columns_version {config['columns_version']} has the VEP "
                                      f"{VEP_PLUGIN_MIN_VERSION} plugin columns but ANNOTATION_VEP_VERSION is "
                                      f"{vep_version}",
                                      hint="Install VEP 116 or call pin_annotation_to_columns_version_4(ANNOTATION).",
                                      id="variantgrid.E013"))
        if config["annotation_consortium"] not in consortia:
            messages.append(Error(f"{build_name} annotation_consortium {config['annotation_consortium']!r} is not one "
                                  f"of {sorted(consortia)}", id="variantgrid.E014"))
    if settings.ANNOTATION_VEP_ENSEMBL_GENCODE not in GENCODE_VALUES:
        messages.append(Error(f"ANNOTATION_VEP_ENSEMBL_GENCODE {settings.ANNOTATION_VEP_ENSEMBL_GENCODE!r} is not one "
                              f"of {GENCODE_VALUES}", id="variantgrid.E014"))

    fork = settings.ANNOTATION_VEP_FORK
    for pipeline, buffer_size in settings.ANNOTATION_VEP_BUFFER_SIZE.items():
        if buffer_size <= 2 * fork:
            messages.append(Error(f"ANNOTATION_VEP_BUFFER_SIZE[{pipeline!r}] = {buffer_size} must be well above "
                                  f"2 x ANNOTATION_VEP_FORK ({fork}) - VEP divides the buffer between forks",
                                  id="variantgrid.E015"))
    batch_max = settings.ANNOTATION_VEP_BATCH_MAX
    if batch_max is not None and settings.ANNOTATION_VEP_BATCH_MIN > batch_max:
        messages.append(Error(f"ANNOTATION_VEP_BATCH_MIN ({settings.ANNOTATION_VEP_BATCH_MIN}) is above "
                              f"ANNOTATION_VEP_BATCH_MAX ({batch_max})", id="variantgrid.E016"))
    if settings.ANNOTATION_MAX_BENIGN_RANKSCORE >= settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE:
        messages.append(Error(f"ANNOTATION_MAX_BENIGN_RANKSCORE ({settings.ANNOTATION_MAX_BENIGN_RANKSCORE}) must be "
                              f"below ANNOTATION_MIN_PATHOGENIC_RANKSCORE "
                              f"({settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE})", id="variantgrid.E017"))
    return messages


@register(Tags.compatibility)
def check_ordered_thresholds(app_configs, **kwargs):
    """Pairs of numbers whose order the code relies on. A heartbeat longer than its lease means every live
    run is reclaimed as dead; a per-operation timeout above the global one never fires."""
    messages = []

    def below(lower_name, upper_name, message_id):
        lower, upper = getattr(settings, lower_name), getattr(settings, upper_name)
        if lower >= upper:
            messages.append(Error(f"{lower_name} ({lower}) must be below {upper_name} ({upper})", id=message_id))

    below("ANNOTATION_RUN_LEASE_HEARTBEAT_SECONDS", "ANNOTATION_RUN_LEASE_SECONDS", "variantgrid.E018")
    below("ANNOTATION_COUNT_LEASE_SECONDS", "ANNOTATION_RUN_LEASE_SECONDS", "variantgrid.E018")
    below("UPLOAD_STEP_LEASE_HEARTBEAT_SECONDS", "UPLOAD_STEP_LEASE_SECONDS", "variantgrid.E018")

    for name in ("MAJOR_OPERATION_STATEMENT_TIMEOUT_SECONDS", "VARIANT_DETAILS_NEARBY_STATEMENT_TIMEOUT_SECONDS"):
        if getattr(settings, name) > settings.DATABASE_STATEMENT_TIMEOUT_SECONDS:
            messages.append(Error(f"{name} ({getattr(settings, name)}) is above DATABASE_STATEMENT_TIMEOUT_SECONDS "
                                  f"({settings.DATABASE_STATEMENT_TIMEOUT_SECONDS}) so never applies",
                                  id="variantgrid.E019"))
    if settings.MAJOR_OPERATION_SLOT_EXPIRE_SECONDS < settings.MAJOR_OPERATION_STATEMENT_TIMEOUT_SECONDS:
        messages.append(Error(f"MAJOR_OPERATION_SLOT_EXPIRE_SECONDS ({settings.MAJOR_OPERATION_SLOT_EXPIRE_SECONDS}) "
                              f"is below MAJOR_OPERATION_STATEMENT_TIMEOUT_SECONDS "
                              f"({settings.MAJOR_OPERATION_STATEMENT_TIMEOUT_SECONDS}) - a slot can expire while its "
                              f"query is still running", id="variantgrid.E019"))

    if settings.VCF_IMPORT_REF_MISMATCH_WARN_FRACTION > settings.VCF_IMPORT_REF_MISMATCH_FAIL_FRACTION:
        messages.append(Error(f"VCF_IMPORT_REF_MISMATCH_WARN_FRACTION ({settings.VCF_IMPORT_REF_MISMATCH_WARN_FRACTION}) "
                              f"is above VCF_IMPORT_REF_MISMATCH_FAIL_FRACTION "
                              f"({settings.VCF_IMPORT_REF_MISMATCH_FAIL_FRACTION})", id="variantgrid.E020"))

    thresholds = settings.VARIANT_GRID_GENOTYPE_QUALITY_THRESHOLDS
    gq, pl = thresholds["gq"], thresholds["pl"]
    if gq["bad"] >= gq["marginal"]:
        messages.append(Error(f"VARIANT_GRID_GENOTYPE_QUALITY_THRESHOLDS gq: bad ({gq['bad']}) must be below marginal "
                              f"({gq['marginal']}) - GQ is a minimum", id="variantgrid.E021"))
    if pl["bad"] <= pl["marginal"]:
        messages.append(Error(f"VARIANT_GRID_GENOTYPE_QUALITY_THRESHOLDS pl: bad ({pl['bad']}) must be above marginal "
                              f"({pl['marginal']}) - PL reads the other way", id="variantgrid.E021"))

    if settings.VARIANT_SYMBOLIC_ALT_ENABLED and \
            settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE < settings.VARIANT_SYMBOLIC_ALT_SIZE:
        messages.append(Warning(f"ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE ({settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE}) "
                                f"is below VARIANT_SYMBOLIC_ALT_SIZE ({settings.VARIANT_SYMBOLIC_ALT_SIZE}) - nothing "
                                f"symbolic is that small, so the gap does nothing", id="variantgrid.W003"))
    return messages


@register(Tags.compatibility)
def check_enumerations(app_configs, **kwargs):
    """Settings whose value is one of a fixed set - a typo is otherwise a KeyError deep in a pipeline."""
    messages = []
    beacon = settings.BEACON_CONFIG
    checks = [
        ("ALLELE_ORIGIN_NOT_PROVIDED_BUCKET", settings.ALLELE_ORIGIN_NOT_PROVIDED_BUCKET, AlleleOriginBucket.values),
        ("ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD", settings.ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD,
         SV_OVERLAP_METHODS),
        ("MME_GENE_ID_PREFERENCE", settings.MME_GENE_ID_PREFERENCE, MME_GENE_ID_PREFERENCES),
        ("BEACON_CONFIG['environment']", beacon.get("environment"), BEACON_ENVIRONMENTS),
        ("BEACON_CONFIG['default_granularity']", beacon.get("default_granularity"), BEACON_GRANULARITIES),
        ("BEACON_CONFIG['max_granularity']", beacon.get("max_granularity"), BEACON_GRANULARITIES),
        ("CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL", settings.CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL,
         [level.value for level in ShareLevel]),
    ]
    for name, value, allowed in checks:
        if value not in allowed:
            messages.append(Error(f"{name} = {value!r} is not one of {list(allowed)}", id="variantgrid.E022"))

    if bad_levels := set(settings.ANALYSIS_SAMPLE_NODE_LEVELS) - set(SampleSourceLevel.values):
        messages.append(Error(f"ANALYSIS_SAMPLE_NODE_LEVELS has {sorted(bad_levels)}, not in {SampleSourceLevel.values}",
                              id="variantgrid.E022"))
    default_granularity, max_granularity = beacon.get("default_granularity"), beacon.get("max_granularity")
    if default_granularity in BEACON_GRANULARITIES and max_granularity in BEACON_GRANULARITIES and \
            BEACON_GRANULARITIES.index(default_granularity) > BEACON_GRANULARITIES.index(max_granularity):
        messages.append(Error(f"BEACON_CONFIG default_granularity {default_granularity!r} is above max_granularity "
                              f"{max_granularity!r}", id="variantgrid.E022"))
    return messages


@register(Tags.compatibility)
def check_dotted_paths(app_configs, **kwargs):
    """Every dotted path a setting holds imports."""
    messages = []
    for setting_name, kind in DOTTED_PATH_SETTINGS.items():
        value = getattr(settings, setting_name)
        paths = [value] if kind == "string" else list(value)
        for dotted_path in paths:
            if not dotted_path:
                continue
            try:
                if kind == "modules":
                    import_module(dotted_path)
                else:
                    import_string(dotted_path)
            except ImportError as e:
                messages.append(Error(f"{setting_name}: {dotted_path!r} does not import ({e})", id="variantgrid.E023"))
    if settings.USE_OIDC and (logout_method := getattr(settings, "OIDC_OP_LOGOUT_URL_METHOD", None)):
        try:
            import_string(logout_method)
        except ImportError as e:
            messages.append(Error(f"OIDC_OP_LOGOUT_URL_METHOD: {logout_method!r} does not import ({e})",
                                  id="variantgrid.E023"))
    return messages


@register(Tags.urls)
def check_url_register(app_configs, **kwargs):
    """A URLS_NAME_REGISTER key that isn't a URL name silently fails to show or hide the page. A Warning, as the
    name may belong to an app this deployment turns off in URLS_APP_REGISTER."""
    messages = []
    url_names = {name for name in get_resolver().reverse_dict if isinstance(name, str)}
    url_names.update(settings.URLS_NAME_REGISTER_SETTINGS)
    register_names = set(settings.URLS_NAME_REGISTER) | set(settings.URLS_NAME_REGISTER_SETTINGS)
    if unknown := sorted(register_names - url_names):
        messages.append(Warning(f"URLS_NAME_REGISTER names that are not a URL name: {', '.join(unknown)}",
                                hint="A typo, a renamed URL, or an app turned off in URLS_APP_REGISTER.",
                                id="variantgrid.W004"))
    installed = set(settings.INSTALLED_APPS) | {app.split('.')[0] for app in settings.INSTALLED_APPS}
    if unknown_apps := sorted(set(settings.URLS_APP_REGISTER) - installed):
        messages.append(Warning(f"URLS_APP_REGISTER names that are not an installed app: {', '.join(unknown_apps)}",
                                id="variantgrid.W004"))
    return messages


@register(Tags.compatibility)
def check_gene_grid_columns(app_configs, **kwargs):
    """The gene grid zips these four lists."""
    lengths = {name: len(getattr(settings, name)) for name in ("GENE_GRID_ENRICHMENT_KIT_COLUMNS",
                                                                "GENE_GRID_ENRICHMENT_KIT_COLUMN_TOOL_TIPS",
                                                                "GENE_GRID_ENRICHMENT_KIT_COLUMN_LABELS",
                                                                "GENE_GRID_ENRICHMENT_KIT_COLUMN_LABEL_TOOL_TIPS")}
    if len(set(lengths.values())) > 1:
        return [Error(f"GENE_GRID_ENRICHMENT_KIT_COLUMN* lists must be the same length: {lengths}",
                      id="variantgrid.E024")]
    return []


@register(Tags.compatibility)
def check_public_sample_gene_matrix(app_configs, **kwargs):
    """The public (no login) matrix has no user to take a build from."""
    if settings.URLS_NAME_REGISTER["public_global_sample_gene_matrix"] and \
            len(_enabled_builds()) > 1 and not settings.PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD:
        return [Error("PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD must be set when more than one genome build is enabled",
                      hint="Or turn off 'public_global_sample_gene_matrix' in URLS_NAME_REGISTER.",
                      id="variantgrid.E025")]
    return []


@register(Tags.security)
def check_served_configuration(app_configs, **kwargs):
    """Settings that only make sense one way on a served deployment."""
    messages = []
    if not settings.DEBUG and not settings.UNIT_TEST:
        for name in ("CACHE_GENERATED_FILES", "VARIANT_ANNOTATION_VALIDATE"):
            if not getattr(settings, name):
                messages.append(Warning(f"{name} is False with DEBUG off - it is a development-only switch",
                                        id="variantgrid.W005"))
    if getattr(settings, "SECURE_PROXY_SSL_HEADER", None):
        insecure = [name for name in ("SESSION_COOKIE_SECURE", "CSRF_COOKIE_SECURE") if not getattr(settings, name)]
        if insecure:
            messages.append(Warning(f"SECURE_PROXY_SSL_HEADER is set (served over HTTPS) but {', '.join(insecure)} "
                                    f"is False", hint="Import variantgrid.settings.components.https_settings.",
                                    id="variantgrid.W006"))
    if 'drf_spectacular' in settings.INSTALLED_APPS and \
            settings.REST_FRAMEWORK.get("DEFAULT_SCHEMA_CLASS") != SPECTACULAR_SCHEMA_CLASS:
        messages.append(Error(f"REST_FRAMEWORK['DEFAULT_SCHEMA_CLASS'] must be {SPECTACULAR_SCHEMA_CLASS!r} - "
                              f"drf-spectacular asserts on any other, so /api/schema returns 500",
                              hint="An env file overriding REST_FRAMEWORK should start from {**REST_FRAMEWORK, ...}.",
                              id="variantgrid.E026"))
    if AXES_MIDDLEWARE in settings.MIDDLEWARE and settings.MIDDLEWARE[-1] != AXES_MIDDLEWARE:
        messages.append(Error(f"{AXES_MIDDLEWARE} must be the last entry in MIDDLEWARE - it renders the lockout response "
                              f"after every other middleware has run",
                              hint="An env file adding middleware should insert before it, as shariantcommon.py does.",
                              id="variantgrid.E027"))
    return messages
