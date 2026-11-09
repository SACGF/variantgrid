// @ts-check
// snpdb/templates/snpdb/genomics/view_contig.html, genes/templates/genes/view_gene_symbol.html
// A horizontal ideogram of one contig, drawn as inline SVG from its file-level cytobands
// (snpdb.serializers.CytobandSerializer output - start/end 1-based inclusive), with an optional red outline around a
// selected range and an optional arrow marking a feature (a gene) too small to see at chromosome scale

const IDEOGRAM_HEIGHT = 24;
const IDEOGRAM_LABEL_HEIGHT = 14;
const IDEOGRAM_MIN_LABEL_WIDTH = 34;
const IDEOGRAM_MARKER_HEIGHT = 10;

class Ideogram {
    /**
     * @param {HTMLElement} container - drawn at its width, and redrawn when the window resizes
     * @param {{contig: {length: number}, cytobands: Array<Object>,
     *          onBandClick?: function(Object): void,
     *          marker?: {start: number, end: number, label: string}}} options
     */
    constructor(container, options) {
        this.container = container;
        this.contigLength = options.contig.length;
        this.cytobands = options.cytobands;
        this.onBandClick = options.onBandClick || null;
        this.marker = options.marker || null;
        this.selection = null;  // {start, end} or null
        this.selectionRect = null;

        let resizeTimer = null;
        $(window).on("resize", () => {
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(() => this.draw(), 100);
        });
    }

    /** @param {?{start: number, end: number}} range - 1-based inclusive, null to clear */
    setSelection(range) {
        this.selection = range;
        this._drawSelection();
    }

    draw() {
        if (!this.cytobands.length) {
            return;
        }
        const width = this.container.clientWidth;
        const scale = width / this.contigLength;
        const top = 2 + (this.marker ? IDEOGRAM_MARKER_HEIGHT : 0);
        const middle = top + IDEOGRAM_HEIGHT / 2;
        const height = top + IDEOGRAM_HEIGHT + IDEOGRAM_LABEL_HEIGHT + 2;

        const svg = svgElement("svg", {width: width, height: height, class: "ideogram"});
        const defs = svgElement("defs", {});
        const hatch = svgElement("pattern", {id: "ideogram-hatch", width: 4, height: 4, patternUnits: "userSpaceOnUse",
                                             patternTransform: "rotate(45)"});
        hatch.appendChild(svgElement("rect", {width: 4, height: 4, fill: "#ffffff"}));
        hatch.appendChild(svgElement("line", {x1: 0, y1: 0, x2: 0, y2: 4, stroke: "#888888", "stroke-width": 2}));
        defs.appendChild(hatch);
        svg.appendChild(defs);

        for (const cytoband of this.cytobands) {
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
            title.textContent = `${cytoband.label} (${cytoband.stain}) ${cytoband.start.toLocaleString()}-${cytoband.end.toLocaleString()}`;
            shape.appendChild(title);
            if (this.onBandClick) {
                const onBandClick = this.onBandClick;
                shape.addEventListener("click", () => onBandClick(cytoband));
            }
            svg.appendChild(shape);

            if (bandWidth >= IDEOGRAM_MIN_LABEL_WIDTH) {
                const label = svgElement("text", {x: x + bandWidth / 2, y: top + IDEOGRAM_HEIGHT + IDEOGRAM_LABEL_HEIGHT,
                                                  "text-anchor": "middle", class: "ideogram-label"});
                label.textContent = cytoband.name;
                svg.appendChild(label);
            }
        }

        this.selectionRect = svgElement("rect", {y: top - 1, height: IDEOGRAM_HEIGHT + 2, fill: "none",
                                                 stroke: "#ff0000", "stroke-width": 2, "pointer-events": "none"});
        svg.appendChild(this.selectionRect);
        if (this.marker) {
            svg.appendChild(this._markerElement(scale, top));
        }
        this.container.replaceChildren(svg);
        this._drawSelection();
    }

    _markerElement(scale, top) {
        const group = svgElement("g", {class: "ideogram-marker", "pointer-events": "none"});
        const title = svgElement("title", {});
        title.textContent = this.marker.label;
        group.appendChild(title);
        const start = (this.marker.start - 1) * scale;
        const spanWidth = Math.max((this.marker.end - this.marker.start + 1) * scale, 2);
        group.appendChild(svgElement("rect", {x: start, y: top, width: spanWidth, height: IDEOGRAM_HEIGHT,
                                              fill: "#ff0000", "fill-opacity": 0.6}));
        const x = start + spanWidth / 2;
        const half = IDEOGRAM_MARKER_HEIGHT / 2;
        const arrow = [[x - half, top - IDEOGRAM_MARKER_HEIGHT], [x + half, top - IDEOGRAM_MARKER_HEIGHT], [x, top]];
        group.appendChild(svgElement("polygon", {points: arrow.map((p) => p.join(",")).join(" "), fill: "#ff0000"}));
        return group;
    }

    _drawSelection() {
        if (!this.selectionRect) {
            return;
        }
        if (!this.selection) {
            this.selectionRect.setAttribute("visibility", "hidden");
            return;
        }
        const scale = this.container.clientWidth / this.contigLength;
        this.selectionRect.setAttribute("x", String((this.selection.start - 1) * scale));
        this.selectionRect.setAttribute("width", String(Math.max((this.selection.end - this.selection.start + 1) * scale, 2)));
        this.selectionRect.setAttribute("visibility", "visible");
    }
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
