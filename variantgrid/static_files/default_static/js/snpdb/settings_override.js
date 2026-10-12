// @ts-check
// snpdb/templates/snpdb/tags/settings_override.html
/* The collection a select points at is a page of its own - link straight to it, and to the listing
   page when nothing is chosen. Used at user, lab and organization level (this tag renders all three) */
function setupCollectionCrossLink(fieldName, urlFunc, listingUrl, listingText) {
    const select = $(`#id_${fieldName}`);
    if (!select.length) {
        return;
    }
    const crossLink = $("<a class='hover-link' target='_blank' rel='noopener'></a>");
    const listingLink = $(`<a class='hover-link' target='_blank' rel='noopener' href='${listingUrl}'>${listingText}</a>`);
    select.after(listingLink);
    select.after(crossLink);

    function updateCrossLink() {
        const pk = select.val();
        crossLink.text(select.find("option:selected").text());
        setCrossLink(crossLink, urlFunc, pk);
        listingLink.toggle(!pk);
    }

    select.change(updateCrossLink);
    updateCrossLink();
}

$(document).ready(function() {
    setupCollectionCrossLink("columns", Urls.view_custom_columns, Urls.custom_columns(), "Manage Custom Columns");
    setupCollectionCrossLink("tag_config", Urls.view_tag_config_collection, Urls.tag_settings(), "Manage Tag Config");
});
