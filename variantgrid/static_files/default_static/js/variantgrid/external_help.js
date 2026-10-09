// @ts-check
// variantgrid/templates/default_templates/external_help.html
/* global populateSideMenu */ // help/help_menu.js
const DEFAULT_PAGE = 'about.html';

function getQueryParams() {
    let qs = document.location.search;
    qs = qs.split('+').join(' ');

    const params = {};
    let tokens;
    const re = /[?&]?([^=]+)=([^&]*)/g;

    while (tokens = re.exec(qs)) {
        params[decodeURIComponent(tokens[1])] = decodeURIComponent(tokens[2]);
    }

    return params;
}

// There are 2 copies of this function, here and internal page.
function helpPage(pageName) {
    const static_iframe = $("iframe#static-iframe");
    const static_iframe_url = "/static/help/" + pageName;
    static_iframe.attr("src", static_iframe_url);
} 


$(document).ready(function() {
    const menu = $('ul#help-side-menu');
    const baseUrl = window.location.pathname + "?page=";
    populateSideMenu(baseUrl, menu);

    const args = getQueryParams();
    const pageName = args["page"] || DEFAULT_PAGE;
    if (pageName) {
        helpPage(pageName);
    }
}); 
