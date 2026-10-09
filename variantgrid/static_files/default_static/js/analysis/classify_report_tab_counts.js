// @ts-check
// analysis/templates/analysis/classify_report_tab_counts.html
// The Classify & Report tab's own label carries the counts, so the page says whether there is anything
// to do without the tab being opened. Fetched after render - working out which taggings are this case's
// walks every analysis its samples are in, which is too much for page load
$(document).ready(function() {
    const pageData = readJsonData("classify-report-tab-counts-data");
    const link = $("a.nav-link[data-href='" + pageData.classify_report_tab_url + "']");
    if (!link.length) {
        return;
    }

    function countBadge(iconClass, count, title) {
        return $("<span>", {class: "ml-2", title: title})
            .append($("<i>", {class: iconClass + " mr-1"}), count);
    }

    function loadCounts() {
        $.getJSON(pageData.classify_report_summary_url, function(data) {
            const counts = $("<span>", {class: "classify-report-counts small"});
            if (data.outstanding) {
                counts.append(countBadge("fa-solid fa-tags", data.outstanding,
                                         "Tags awaiting classification"));
            }
            if (data.classifications) {
                counts.append(countBadge("fa-solid fa-file-medical", data.classifications,
                                         "Classifications on this case"));
            }
            $(".classify-report-counts", link).remove();
            link.append(counts);
        });
    }

    // Classifying inside the tab reloads it - the label has to follow, or it keeps the counts you arrived with
    $(document).on("classifyReportChanged", loadCounts);
    loadCounts();
});
