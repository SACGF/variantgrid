// SVG pedigree node rendering for TrioNode and QuadNode.
//
// Shapes follow standard pedigree notation:
//   rect   = male
//   circle = female
//   filled = affected, white = unaffected
//
// Sex convention matches samplenode.js: sex === 'F' → circle, anything else → rect.
// Father/mother sex is definitional (always rect/circle respectively).
// Proband is always rendered filled (they are the affected individual under study).

const PEDIGREE_SHAPE_SIZE    = 14;
const PEDIGREE_AFFECTED_FILL  = "#222";
const PEDIGREE_UNAFFECTED_FILL = "#fff";
const PEDIGREE_STROKE       = "#333";
const PEDIGREE_STROKE_WIDTH  = 1.5;

function pedigreeRect(svg, cx, cy, size, filled) {
    const half = size / 2;
    svg.append("svg:rect")
       .attr("x", cx - half).attr("y", cy - half)
       .attr("width", size).attr("height", size)
       .attr("stroke", PEDIGREE_STROKE).attr("stroke-width", PEDIGREE_STROKE_WIDTH)
       .attr("fill", filled ? PEDIGREE_AFFECTED_FILL : PEDIGREE_UNAFFECTED_FILL);
}

function pedigreeCircle(svg, cx, cy, size, filled) {
    svg.append("svg:circle")
       .attr("cx", cx).attr("cy", cy).attr("r", size / 2)
       .attr("stroke", PEDIGREE_STROKE).attr("stroke-width", PEDIGREE_STROKE_WIDTH)
       .attr("fill", filled ? PEDIGREE_AFFECTED_FILL : PEDIGREE_UNAFFECTED_FILL);
}

function pedigreeShape(svg, cx, cy, size, sex, filled) {
    if (sex === 'F') {
        pedigreeCircle(svg, cx, cy, size, filled);
    } else {
        pedigreeRect(svg, cx, cy, size, filled);
    }
}

function pedigreeLine(svg, x1, y1, x2, y2) {
    svg.append("svg:line")
       .attr("x1", x1).attr("y1", y1).attr("x2", x2).attr("y2", y2)
       .attr("stroke", PEDIGREE_STROKE).attr("stroke-width", PEDIGREE_STROKE_WIDTH);
}


// ── TrioNode ──────────────────────────────────────────────────────────────────
//
//   ■ ———————— ●       father (rect) left, mother (circle) right
//        |
//        ■/●           proband (sex-aware), always filled

const TRIO_W = 64, TRIO_H = 50;

function createTrioNode() {
    const div = $('<div/>').addClass("window default-node-container");

    const svgEl = $('<div class="trio-pedigree-svg"/>').css({position: 'absolute', top: 0, left: 0});
    div.append(svgEl);

    const nodeOverlay = $('<div/>').addClass("node-overlay").css({
        width: '100%', height: '100%', position: 'absolute', top: 0, left: 0, 'z-index': 30,
    });
    nodeOverlay.append("<span class='node-name'></span>");
    $("<div />", {class: "node-color-overlay"}).appendTo(nodeOverlay);
    div.append(nodeOverlay);

    div[0].updateState = function(args) {
        $('svg', this).remove();
        const svg = d3.select($('.trio-pedigree-svg', this)[0])
            .append("svg:svg").attr("width", TRIO_W).attr("height", TRIO_H);

        const fatherX = 16, motherX = 48, parentsY = 16;
        const probandX = 32, probandY = 40;
        const midX = (fatherX + motherX) / 2;

        pedigreeLine(svg, fatherX, parentsY, motherX, parentsY);
        pedigreeLine(svg, midX, parentsY, probandX, probandY);

        pedigreeRect(svg, fatherX, parentsY, PEDIGREE_SHAPE_SIZE, args.father_affected);
        pedigreeCircle(svg, motherX, parentsY, PEDIGREE_SHAPE_SIZE, args.mother_affected);
        pedigreeShape(svg, probandX, probandY, PEDIGREE_SHAPE_SIZE, args.proband_sex, true);
    };

    return div;
}


// ── QuadNode ──────────────────────────────────────────────────────────────────
//
//   ■ ————————————— ●   father (rect) left, mother (circle) right
//            |
//      ■/● ——|—— ■/●   sibling (left, sex-aware, affected flag),
//                       proband (right, sex-aware, always filled)

const QUAD_W = 80, QUAD_H = 55;

function createQuadNode() {
    const div = $('<div/>').addClass("window default-node-container");

    const svgEl = $('<div class="quad-pedigree-svg"/>').css({position: 'absolute', top: 0, left: 0});
    div.append(svgEl);

    const nodeOverlay = $('<div/>').addClass("node-overlay").css({
        width: '100%', height: '100%', position: 'absolute', top: 0, left: 0, 'z-index': 30,
    });
    nodeOverlay.append("<span class='node-name'></span>");
    $("<div />", {class: "node-color-overlay"}).appendTo(nodeOverlay);
    div.append(nodeOverlay);

    div[0].updateState = function(args) {
        $('svg', this).remove();
        const svg = d3.select($('.quad-pedigree-svg', this)[0])
            .append("svg:svg").attr("width", QUAD_W).attr("height", QUAD_H);

        const fatherX = 16, motherX = 64, parentsY = 14;
        const siblingX = 24, probandX = 56, childrenY = 42;
        const dropX = (fatherX + motherX) / 2;
        const barY = (parentsY + childrenY) / 2;

        pedigreeLine(svg, fatherX, parentsY, motherX, parentsY);
        pedigreeLine(svg, dropX, parentsY, dropX, barY);
        pedigreeLine(svg, siblingX, barY, probandX, barY);
        pedigreeLine(svg, siblingX, barY, siblingX, childrenY);
        pedigreeLine(svg, probandX, barY, probandX, childrenY);

        pedigreeRect(svg, fatherX, parentsY, PEDIGREE_SHAPE_SIZE, args.father_affected);
        pedigreeCircle(svg, motherX, parentsY, PEDIGREE_SHAPE_SIZE, args.mother_affected);
        pedigreeShape(svg, siblingX, childrenY, PEDIGREE_SHAPE_SIZE, args.sibling_sex, args.sibling_affected);
        pedigreeShape(svg, probandX, childrenY, PEDIGREE_SHAPE_SIZE, args.proband_sex, true);
    };

    return div;
}
