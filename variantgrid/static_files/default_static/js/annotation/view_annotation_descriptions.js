// @ts-check
// annotation/templates/annotation/view_annotation_descriptions.html
// Each example cell is drawn by the grid's own client renderer from the column definition and
// fictional row the view built - the same path DataTableDefinition takes, minus DataTables
// @see snpdb/grid_columns/composite_examples.py
const EXAMPLE_EXTRA = readJsonData("view-annotation-descriptions-data").example_extra;

function setupExampleCell(table) {
    const col = table.data("columnJson");
    const row = table.data("rowJson");
    const th = table.find("thead th");
    th.text(col.label);
    if (col.headerTitle) {
        th.attr("title", col.headerTitle);
    }
    if (col.width) {
        th.css("width", col.width);
    }

    let html = row[col.data];
    if (col.render) {
        html = eval(col.render)(row[col.data], "display", row, {extra: EXAMPLE_EXTRA, kwargs: col.renderKwargs || null});
        if (html instanceof jQuery) {
            html = html.prop("outerHTML");
        }
    }
    const td = table.find("tbody td");
    td.html(html == null ? "" : String(html));
    // The cell is a picture of a fictional variant - nothing in it goes anywhere
    td.find("a").removeAttr("href").removeAttr("target").removeAttr("orig_href");

    if (!col.sortMenu || !col.sortMenu.length) {
        return;
    }
    const menu = $('<span>', {class: 'dt-sort-menu'});
    const toggle = $('<a>', {class: 'dt-sort-menu-toggle', href: 'javascript:void(0)',
                             title: 'Sort this column by', text: '\u25be'});
    const items = $('<div>', {class: 'dropdown-menu dt-sort-menu-items'});
    for (const entry of col.sortMenu) {
        $('<a>', {class: 'dropdown-item', href: 'javascript:void(0)', text: entry.label})
            .on('click', function(event) {
                event.stopPropagation();
                FloatingPanel.hide();
            }).appendTo(items);
    }
    toggle.on('click', function(event) {
        event.stopPropagation();
        if (FloatingPanel.isShowing(items)) {
            FloatingPanel.hide();
        } else {
            FloatingPanel.show(items, this, {alignRight: true});
        }
    });
    menu.append(toggle);
    th.append(menu);
}

// Each row carries the columns version range its definition is written in - a row without one
// is written in every version. Sections left with nothing to show for the chosen version go too
function showColumnsVersion(version) {
    const showAll = version === "all";
    $(".annotation-row").each(function() {
        const row = $(this);
        const minVersion = parseInt(row.attr("data-min-cv"));
        const maxVersion = parseInt(row.attr("data-max-cv"));
        row.toggle(showAll || ((isNaN(minVersion) || version >= minVersion) &&
                               (isNaN(maxVersion) || version <= maxVersion)));
    });
    $(".annotation-columns-version-note").toggle(showAll);
    $(".composite-section, .annotation-level-card").each(function() {
        $(this).toggle($(this).find(".annotation-row:visible").length > 0);
    });
    $(".columns-version-toggle .btn").removeClass("active");
    $(".columns-version-toggle .btn").filter(function() {
        return $(this).data("columnsVersion") === version;
    }).addClass("active");
}

$(document).ready(() => {
    $("table.composite-example").each(function() {
        setupExampleCell($(this));
    });
    $(".columns-version-toggle .btn").on("click", function() {
        showColumnsVersion($(this).data("columnsVersion"));
    });
    showColumnsVersion(readJsonData("view-annotation-descriptions-data").latest_columns_version);
});
