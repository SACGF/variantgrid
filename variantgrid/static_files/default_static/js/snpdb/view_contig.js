// @ts-check
// snpdb/templates/snpdb/genomics/view_contig.html
// An ideogram of the contig's cytobands. The selected band (?band=, also a whole arm or a range like p13-p12)
// is shared by the ideogram (ideogram.js), the cytoband table and the genes and variants grid tabs, and changes
// without a reload

let contigPage = null;  // contig_page_data from snpdb.views.views.view_contig: contig and cytobands are
                        // snpdb.serializers.ContigSerializer / CytobandSerializer output
let selectedBand = null;  // {band, start, end} - 1-based inclusive - or null for the whole contig
let contigIdeogram = null;

function initViewContig(contigPageData) {
    contigPage = contigPageData;
    if (contigPage.selectedInterval) {
        selectedBand = Object.assign({band: contigPage.band}, contigPage.selectedInterval);
    }
    $(document).ready(setupViewContig);
}

// data-datatable-data callback for the genes grid
function contigGenesFilter(data) {
    data.contig_id = contigPage.contig.id;
    if (selectedBand) {
        data.region_start = selectedBand.start;
        data.region_end = selectedBand.end;
    }
}

function setupViewContig() {
    const container = document.getElementById("contig-ideogram");
    if (container) {
        contigIdeogram = new Ideogram(container, {
            contig: contigPage.contig,
            cytobands: contigPage.cytobands,
            onBandClick: (cytoband) => selectBand({band: cytoband.name, start: cytoband.start, end: cytoband.end}),
        });
        contigIdeogram.draw();
    }

    $("#cytoband-table").on("click", "a.cytoband-link", function(event) {
        if (event.ctrlKey || event.metaKey || event.shiftKey || event.button !== 0) {
            return;  // a new tab still gets the link
        }
        event.preventDefault();
        selectBand(cytobandForRow($(this).closest("tr")));
    });
    $("#selected-band-clear").click(() => selectBand(null));
    // The grids are built while their tab is hidden, so the shown one needs its column widths measured
    $("#contig-grids a[data-toggle='tab']").on("shown.bs.tab", () => {
        $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
    });

    history.replaceState({selectedBand: selectedBand}, "");
    window.addEventListener("popstate", (event) => {
        if (event.state && "selectedBand" in event.state) {
            showSelectedBand(event.state.selectedBand);
        }
    });
    showSelectedBand(selectedBand, false);
}

function cytobandForRow(row) {
    return {band: row.data("band"), start: Number(row.data("start")), end: Number(row.data("end"))};
}

function selectBand(band) {
    const url = new URL(window.location.href);
    if (band) {
        url.searchParams.set("band", band.band);
    } else {
        url.searchParams.delete("band");
    }
    history.pushState({selectedBand: band}, "", url.toString());
    showSelectedBand(band);
}

function showSelectedBand(band, refilterGrids = true) {
    selectedBand = band;
    if (contigIdeogram) {
        contigIdeogram.setSelection(band);
    }

    $("#cytoband-table tbody tr").each(function() {
        const row = $(this);
        const inBand = Boolean(band) && band.start <= Number(row.data("start")) && Number(row.data("end")) <= band.end;
        row.toggleClass("table-active", inBand);
    });

    $("#selected-band-summary").toggle(Boolean(band));
    $("#ideogram-click-hint").toggle(!band);
    $("#unknown-band-card").toggle(!band && Boolean(new URL(window.location.href).searchParams.get("band")));
    if (band) {
        $("#selected-band-name").text(contigPage.contig.name + band.band);
    }

    if (refilterGrids) {
        const genesGrid = $("#contig-genes-grid");
        if ($.fn.DataTable.isDataTable(genesGrid)) {
            genesGrid.DataTable().ajax.reload();
        }
        if (typeof setAllVariantsLockedFilters === "function" && contigPage.numVariantTypes !== undefined) {
            setAllVariantsLockedFilters(contigVariantsLockedFilters(), contigPage.numVariantTypes);
        }
    }
}

function contigVariantsLockedFilters() {
    return {
        contig_ids: [contigPage.contig.id],
        position_range: selectedBand ? [selectedBand.start, selectedBand.end] : null,
    };
}
