// @ts-check
// variantopedia/templates/variantopedia/search.html - the search box and the results' category links

function initSearch() {
    const bigSearch = $('#big_search');
    bigSearch.select();
    bigSearch.focus(() => {
        bigSearch.select();
    });
}

function updateSearchMode() {
    const feelingLucky = $("#feeling-lucky").is(':checked');
    $("#search-mode").val(feelingLucky ? "lucky" : "preview");
}

function scrollToResults(category) {
    const scrollMe = isPhoneMode() ? $('body') : $('.main-content');
    const scrollToOffset = $(`#results-${category}`).position().top;

    scrollMe.animate({scrollTop: scrollToOffset - 30}, 100, 'swing', () => {
        $(`.result-${category}`).animate({'background-color': '#cceecc'}, 150).animate({'background-color': 'inherit'}, 100);
    });
    return false;
}

function searchAgain() {
    const scrollMe = isPhoneMode() ? $('body') : $('.main-content');
    scrollMe.animate({scrollTop: 0}, 100);
    $('#big_search').focus();
}
