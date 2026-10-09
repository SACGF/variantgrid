// @ts-check
// classification/templates/classification/classification_dashboard.html
const dashboardData = readJsonData("classification-dashboard-data");

// for showing issues
function classificationReviews(data) {
    // data.flags = JSON.stringify(['classification_pending_changes', 'classification_internal_review', 'classification_suggestion']);
    data.pending = "true";
    data.labs = dashboardData.lab_ids_str;
}

function classificationSigChanges(data) {
    data.flags = JSON.stringify(['classification_significance_change']);
    data.labs = dashboardData.lab_ids_str;
}

function classificationMatchingVariant(data) {
    data.flags = JSON.stringify(['classification_matching_variant', 'classification_matching_variant_warning', 'classification_transcript_version_change']);
    data.labs = dashboardData.lab_ids_str;
}

function classificationUnshared(data) {
    data.flags = JSON.stringify(['classification_unshared']);
    data.labs = dashboardData.lab_ids_str;
}
function classificationWithdrawn(data) {
    data.flags = JSON.stringify(['classification_withdrawn']);
    data.labs = dashboardData.lab_ids_str;
}
function classificationExcludeClinVar(data) {
    data.flags = JSON.stringify(['classification_not_public']);
    data.labs = dashboardData.lab_ids_str;
}
