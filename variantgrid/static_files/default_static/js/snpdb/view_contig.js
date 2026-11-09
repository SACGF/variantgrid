// @ts-check
// snpdb/templates/snpdb/genomics/view_contig.html
// An ideogram of the contig's cytobands. The selected band (?band=, also a whole arm or a range like p13-p12)
// is shared by the ideogram, the cytoband table and the genes and variants grids, and changes without a reload

const IDEOGRAM_HEIGHT = 24;
const IDEOGRAM_LABEL_HEIGHT = 14;
const IDEOGRAM_MIN_LABEL_WIDTH = 34;

let contigPage = null;  // contig_page_data from snpdb.views.views.view_contig
let selectedBand = null;  // {band, start, end} - 1-based inclusive - or null for the whole contig

function initViewContig(contigPageData) {
    contigPage = contigPageData;
    if (contigPage.selectedInterval) {
        selectedBand = Object.assign({band: contigPage.band}, contigPage.selectedInterval);
    }
    $(document).ready(setupViewContig);
}

// data-datatable-data callback for the genes grid
function contigGenesFilter(data) {
    data.contig_id = contigPage.contigId;
    if (selectedBand) {
        data.region_start = selectedBand.start;
        data.region_end = selectedBand.end;
    }
}

function setupViewContig() {
    drawIdeogram();
    let resizeTimer = null;
    $(window).on("resize", () => {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(drawIdeogram, 100);
    });

    $("#cytoband-table").on("click", "a.cytoband-link", function(event) {
        if (event.ctrlKey || event.metaKey || event.shiftKey || event.button !== 0) {
            return;  // a new tab still gets the link
        }
        event.preventDefault();
        selectBand(cytobandForRow($(this).closest("tr")));
    });
    $("#selected-band-clear").click(() => selectBand(null));

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
    drawSelection();

    $("#cytoband-table tbody tr").each(function() {
        const row = $(this);
        const inBand = Boolean(band) && band.start <= Number(row.data("start")) && Number(row.data("end")) <= band.end;
        row.toggleClass("table-active", inBand);
    });

    $("#selected-band-card").toggle(Boolean(band));
    $("#unknown-band-card").toggle(!band && Boolean(new URL(window.location.href).searchParams.get("band")));
    if (band) {
        $("#selected-band-name").text(contigPage.contigName + band.band);
        $("#selected-band-coordinates").text(contigPage.contigLabel + ":" + band.start.toLocaleString() + "-" + band.end.toLocaleString());
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
        contig_ids: [contigPage.contigId],
        position_range: selectedBand ? [selectedBand.start, selectedBand.end] : null,
    };
}

function stainFill(stain) {
    // gneg is white, gposN a grey N% of the way to black
    const gpos = /^gpos(\d*)$/.exec(stain);
    if (gpos) {
        const percent = gpos[1] ? Number(gpos[1]) : 100;
        const level = Math.round(255 * (1 - percent / 100));
        return `rgb(${level},${level},${level})`;
    }
    if (stain === "acen") {
        return "#d92f27";
    }
    if (stain === "gvar" || stain === "stalk") {
        return "url(#ideogram-hatch)";
    }
    return "#ffffff";
}

function svgElement(name, attributes) {
    const element = document.createElementNS(SVG_NS, name);
    for (const [key, value] of Object.entries(attributes)) {
        element.setAttribute(key, String(value));
    }
    return element;
}

function ideogramScale(width) {
    return width / contigPage.contigLength;
}

function drawIdeogram() {
    const container = document.getElementById("contig-ideogram");
    if (!container || !contigPage.cytobands.length) {
        return;
    }
    const width = container.clientWidth;
    const height = IDEOGRAM_HEIGHT + IDEOGRAM_LABEL_HEIGHT + 4;
    const scale = ideogramScale(width);
    const top = 2;
    const middle = top + IDEOGRAM_HEIGHT / 2;

    const svg = svgElement("svg", {width: width, height: height, class: "ideogram"});
    const defs = svgElement("defs", {});
    const hatch = svgElement("pattern", {id: "ideogram-hatch", width: 4, height: 4, patternUnits: "userSpaceOnUse",
                                         patternTransform: "rotate(45)"});
    hatch.appendChild(svgElement("rect", {width: 4, height: 4, fill: "#ffffff"}));
    hatch.appendChild(svgElement("line", {x1: 0, y1: 0, x2: 0, y2: 4, stroke: "#888888", "stroke-width": 2}));
    defs.appendChild(hatch);
    svg.appendChild(defs);

    for (const cytoband of contigPage.cytobands) {
        // Cytoband coordinates are 1-based inclusive
        const x = (cytoband.start - 1) * scale;
        const bandWidth = Math.max((cytoband.end - cytoband.start + 1) * scale, 1);
        let shape;
        if (cytoband.stain === "acen") {
            // Pinch the centromere - each arm's acen band narrows towards the other arm
            const points = cytoband.name.startsWith("p")
                ? [[x, top], [x + bandWidth, middle], [x, top + IDEOGRAM_HEIGHT]]
                : [[x + bandWidth, top], [x, middle], [x + bandWidth, top + IDEOGRAM_HEIGHT]];
            shape = svgElement("polygon", {points: points.map((p) => p.join(",")).join(" ")});
        } else {
            shape = svgElement("rect", {x: x, y: top, width: bandWidth, height: IDEOGRAM_HEIGHT, stroke: "#333333",
                                        "stroke-width": 0.5});
        }
        shape.setAttribute("fill", stainFill(cytoband.stain));
        shape.setAttribute("class", "ideogram-band");
        shape.setAttribute("data-band", cytoband.name);
        const title = svgElement("title", {});
        title.textContent = `${contigPage.contigName}${cytoband.name} (${cytoband.stain}) ${cytoband.start.toLocaleString()}-${cytoband.end.toLocaleString()}`;
        shape.appendChild(title);
        shape.addEventListener("click", () => selectBand({band: cytoband.name, start: cytoband.start, end: cytoband.end}));
        svg.appendChild(shape);

        if (bandWidth >= IDEOGRAM_MIN_LABEL_WIDTH) {
            const label = svgElement("text", {x: x + bandWidth / 2, y: top + IDEOGRAM_HEIGHT + IDEOGRAM_LABEL_HEIGHT,
                                              "text-anchor": "middle", class: "ideogram-label"});
            label.textContent = cytoband.name;
            svg.appendChild(label);
        }
    }

    svg.appendChild(svgElement("rect", {id: "ideogram-selection", y: top - 1, height: IDEOGRAM_HEIGHT + 2,
                                        fill: "none", stroke: "#ff0000", "stroke-width": 2,
                                        "pointer-events": "none"}));
    container.replaceChildren(svg);
    drawSelection();
}

function drawSelection() {
    const selection = document.getElementById("ideogram-selection");
    const container = document.getElementById("contig-ideogram");
    if (!selection || !container) {
        return;
    }
    if (!selectedBand) {
        selection.setAttribute("visibility", "hidden");
        return;
    }
    const scale = ideogramScale(container.clientWidth);
    selection.setAttribute("x", String((selectedBand.start - 1) * scale));
    selection.setAttribute("width", String(Math.max((selectedBand.end - selectedBand.start + 1) * scale, 2)));
    selection.setAttribute("visibility", "visible");
}
