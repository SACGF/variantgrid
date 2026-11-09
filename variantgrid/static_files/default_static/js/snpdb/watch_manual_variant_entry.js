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

function handleSampleData(mvecId, data) {
    // The sample import succeeds once its variants are annotated
    const sample = data["sample"];
    if (data["import_status"] === 'E' || (sample && sample["import_status"] === 'E')) {
        displayError("Sample import failed.");
    } else if (sample && sample["import_status"] === 'S') {
        window.location.href = Urls.view_sample(sample["id"]);
    } else {
        populatePage(data);
        if (sample) {
            $("#variant-status").text(`Sample ${sample["import_status_display"]}`);
        }
        setTimeout(() => pollServer(mvecId, true), 2000);
    }
}

function handleMVECData(mvecId, asSample, data) {
    if (asSample) {
        handleSampleData(mvecId, data);
    } else if (data["is_ready"]) {
        // if we can redirect - do so
        window.location.href = Urls.view_variant(data["first_variant_id"]);
    } else {
        if (data["import_status"] == 'E') {
            displayError("Variant import failed.");
        } else {
            populatePage(data);
            setTimeout(() => pollServer(mvecId, false), 2000);
        }
    }
}

function pollServer(mvecId, asSample) {
    $.ajax({
        url: Urls.api_manual_variant_entry_collection(mvecId),
        success: (data) => handleMVECData(mvecId, asSample, data),
        error: function() {
            displayError("Error contacting server - please reload page.");
        }
    });
}


