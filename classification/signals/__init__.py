# noinspection PyUnresolvedReferences
from classification.signals.classification_health_checks import *   # so we load the receivers
from classification.signals.classification_hooks_assign_owner import *   # so we load the receivers
from classification.signals.classification_hooks_discordance_notifications import *
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
from classification.signals.classification_hooks_pending_flags import *  # so we load the receivers
