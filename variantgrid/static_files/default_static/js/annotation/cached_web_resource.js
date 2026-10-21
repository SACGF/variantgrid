// @ts-check
// annotation/templates/annotation/cached_web_resource.html
const CACHED_WEB_RESOURCE_POLL_MS = 5000;

function pollCachedWebResource(cwrName) {
    const container = $(`.cached-web-resource-container[cached-web-resource="${cwrName}"]`);
    $.get(Urls.view_cached_web_resource(cwrName), function(html) {
        const newContainer = $($.parseHTML(html)).filter(".cached-web-resource-container");
        container.replaceWith(newContainer);
        // 'C' is the brief window between the resource being re-created and its task starting
        const importStatus = newContainer.attr("import-status");
        if (importStatus === "C" || importStatus === "I") {
            setTimeout(() => pollCachedWebResource(cwrName), CACHED_WEB_RESOURCE_POLL_MS);
        }
    });
}

function initCachedWebResources() {
    $(document).on("click", ".load-cached-web-resource", function() {
        const cwrName = $(this).attr("cached-web-resource");
        const container = $(this).closest(".cached-web-resource-container");

        $.ajax({
            type: "POST",
            url: Urls.load_cached_web_resource(cwrName),
            success: function() {
                container.html('<i class="fa fa-spinner"></i>');
                setTimeout(() => pollCachedWebResource(cwrName), CACHED_WEB_RESOURCE_POLL_MS);
            }
        });
    });

    $('.cached-web-resource-container[import-status="I"]').each(function() {
        const cwrName = $(this).attr("cached-web-resource");
        setTimeout(() => pollCachedWebResource(cwrName), CACHED_WEB_RESOURCE_POLL_MS);
    });
}
