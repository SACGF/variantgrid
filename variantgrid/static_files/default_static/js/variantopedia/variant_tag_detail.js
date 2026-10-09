// @ts-check
// variantopedia/templates/variantopedia/variant_tag_detail.html
function removeTagCallback() {
    // when you click the detail - it closes for some reason, so can't really update it.
    // instead we'll just update the overall grid of counts
    // $('#variant-tag-detail-datatable').DataTable().ajax.reload();
    $('#variant-tags-datatable').DataTable().ajax.reload();
}

function removeTag(variant_tag_id) {
    const tagDetailData = readJsonData("variant-tag-detail-data");
    let data = 'variant_id=' + tagDetailData.variant_id;
    data += '&tag_id=' + tagDetailData.tag_id;
    data += '&variant_tag_id=' + variant_tag_id;
    data += '&op=' + 'del';
    $.ajax({
        type: "POST",
        data: data,
        url: Urls.set_variant_tag('V'),
        success : removeTagCallback,
    });
}

function tagDetailRenderer(data, type, row) {
    let cellValue = null;
    if (row.can_write) {
        const deleteLink = $("<a/>").attr({
            href: `javascript:removeTag(${row.id})`,
        });
        const button = $('<span style="display: inline-block;" class="click-to-delete-button" title="" data-original-title="Remove tag" data-p="1"></span>');
        deleteLink.append(button);
        cellValue = deleteLink.prop('outerHTML');
    }
    return cellValue;
}
