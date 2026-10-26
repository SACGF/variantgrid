# noinspection PyUnresolvedReferences
from classification.signals.classification_health_checks import *   # so we load the receivers
from classification.signals.classification_hooks_assign_owner import *   # so we load the receivers
from classification.signals.classification_hooks_discordance_status import *  # so we load the receivers
from classification.signals.classification_hooks_share_flags import *  # so we load the receivers
from classification.signals.classification_hooks_significant_change import *  # so we load the receivers
from classification.signals.classification_hooks_variants_classification_changes import *  # so we load the receivers
from classification.signals.classification_search import *  # so we load the receivers
from classification.signals.discordance_report_review_detail import *  # so we load the receivers
from classification.signals.classification_hooks_import_notifications import *  # so we load the receivers
from classification.signals.discordance_report_triage_changes import *  # so we load the receivers
from classification.signals.discordance_report_search import *  # so we load the receivers
from classification.signals.classification_liftover import *  # so we load the receivers
from classification.signals.classification_hooks_grouping import *  # so we load the receivers
from classification.signals.classification_hooks_grouping_search_terms import *  # so we load the receivers
from classification.signals.classification_import_integration_status import *  # so we load the receivers
from classification.signals.clinvar_export_integration_status import *  # so we load the receivers
from classification.signals.classification_hooks_allele_merge import *  # so we load the receivers
from classification.signals.classification_hooks_overlaps import *
from classification.signals.clinvar_record_to_contribution import *
from classification.signals.overlap_review_formatting import *
from classification.signals.classification_hooks_citations import *
from classification.signals.classification_flag_extra_info import *  # so we load the receivers
from classification.signals.classification_settings_hooks import *  # so we load the receivers
from classification.signals.classification_hooks_p_hgvs import *  # so we load the receivers
