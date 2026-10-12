// @ts-check
// snpdb/templates/snpdb/data/watch_manual_variant_entry.html
function populatePage(data) {
    console.log(data);
    $("#variant-status").text(data["first_variant_annotation_status"]);
    $("#created-date").empty().append(createTimestampDom(data["created"], true));
}

function displayError(message) {
    $("#glowing-cloud").remove();
    $("#variant-status").addClass("error").text(message);
}

function handleMVECData(mvecId, data) {
    // if we can redirect - do so
    // otherwise populate screen
    if (data["is_ready"]) {
        window.location.href = Urls.view_variant(data["first_variant_id"]);
    } else {
        if (data["import_status"] == 'E') {
            displayError("Variant import failed.");
        } else {
            populatePage(data);
            setTimeout(() => pollServer(mvecId), 2000);
        }
    }
}

function pollServer(mvecId) {
    $.ajax({
        url: Urls.api_manual_variant_entry_collection(mvecId),
        success: (data) => handleMVECData(mvecId, data),
        error: function() {
            displayError("Error contacting server - please reload page.");
        }
    });
}


