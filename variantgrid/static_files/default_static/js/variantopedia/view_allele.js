// @ts-check
// variantopedia/templates/variantopedia/view_allele.html
/* global classificationGroupingRedraw */ // classification/classification_groupings.js
let viewAlleleId = null;  // set by initViewAllele
let alleleOriginFilter = null;
let testingContextFilter = null;

function makeFilterBox(text) {
    $('#classificationFilterBox').empty();
    const box = $('<a>', {class: 'filter-item', text: text, onclick:"clearFilterBox()"});
    $('#classificationFilterBox').addClass('mr-3').append(box);
}
function clearFilterBox() {
    alleleOriginFilter = null;
    testingContextFilter = null;
    updateDiff();
    $('#classificationFilterBox').empty().removeClass('mr-3');
    classificationGroupingRedraw();
}

function filterAlleleOrigin(alleleOrigin, label) {
    testingContextFilter = null;
    alleleOriginFilter = alleleOrigin;
    updateDiff();
    makeFilterBox(label);
    classificationGroupingRedraw();
}

function filterTestingContext(testingContext, label) {
    testingContextFilter = testingContext;
    alleleOriginFilter = null;
    updateDiff();
    makeFilterBox(label);
    classificationGroupingRedraw();
}

function classificationGroupingFilterExtra(data) {
    data.allele_id = viewAlleleId;
    data.allele_origin = alleleOriginFilter;
    data.testing_context = testingContextFilter;
}

function updateDiff() {
    let href = Urls.classification_diff() + "?allele=" + viewAlleleId + "&latest=true";
    if (alleleOriginFilter != null) {
        href += `&allele_origin=${alleleOriginFilter}`;
    }
    if (testingContextFilter != null) {
        href += `&testing_context=${testingContextFilter}`;
    }
    $('#showDiffLink').attr('href', href);
}

function initViewAllele(alleleId) {
    viewAlleleId = alleleId;
    $(document).ready(updateDiff);
}
