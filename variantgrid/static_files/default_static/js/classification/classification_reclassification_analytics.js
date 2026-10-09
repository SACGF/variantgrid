// @ts-check
// classification/templates/classification/classification_reclassification_analytics.html
const reclassificationCharts = readJsonData("classification-reclassification-analytics-data").chart_data;
const AXIS_COLOUR = "#52514e";
const GRID_COLOUR = "#e6e5e1";

function chartLayout(extra) {
    return Object.assign({
        height: 420,
        margin: {t: 20, r: 20, b: 60, l: 70},
        paper_bgcolor: "#fcfcfb",
        plot_bgcolor: "#fcfcfb",
        font: {color: AXIS_COLOUR},
        xaxis: {gridcolor: GRID_COLOUR, zeroline: false},
        yaxis: {gridcolor: GRID_COLOUR, zeroline: false}
    }, extra || {});
}

const plotConfig = {responsive: true, displayModeBar: false};

const timeToReclassification = reclassificationCharts.time_to_reclassification;
if (timeToReclassification.series.length) {
    const edges = timeToReclassification.edges;
    const centres = edges.slice(0, -1).map((edge, index) => (edge + edges[index + 1]) / 2);
    const barWidth = (edges[1] - edges[0]) * 0.94;
    const traces = timeToReclassification.series.map(series => ({
        x: centres,
        y: series.counts,
        width: barWidth,
        name: series.label,
        type: 'bar',
        marker: {color: series.colour, line: {color: "#fcfcfb", width: 1}},
        customdata: series.counts.map((count, index) =>
            [Math.round(Math.pow(10, edges[index]) - 1), Math.round(Math.pow(10, edges[index + 1]) - 1)]),
        hovertemplate: '%{y} records ' + series.label +
                       '<br>%{customdata[0]} to %{customdata[1]} days<extra></extra>'
    }));
    Plotly.newPlot('time-to-reclassification', traces, chartLayout({
        barmode: 'stack',
        bargap: 0,
        legend: {traceorder: 'normal'},
        xaxis: {title: {text: 'Days at the starting significance'}, range: [0, edges[edges.length - 1]],
                tickmode: 'array', tickvals: timeToReclassification.tick_values,
                ticktext: timeToReclassification.tick_labels,
                gridcolor: GRID_COLOUR, zeroline: false},
        yaxis: {title: {text: 'Reclassifications'}, gridcolor: GRID_COLOUR, zeroline: false}
    }), plotConfig);
}

const SANKEY_HEIGHT = 520, SANKEY_MARGIN = 25, SANKEY_PAD = 16;
// plotly places a node around the centre of node.y, and skips any node whose y is a falsy zero
const SANKEY_MIN_Y = 1e-4;

/** Stacks a column of nodes top to bottom in the order given, sized by its share of the flow,
 *  so a significance sits at the same depth on both sides and a downgrade always reads downwards. */
function sankeyColumnCentres(totals, grandTotal, usable, plotHeight) {
    let cursor = 0;
    return totals.map(total => {
        const height = (total / grandTotal) * usable;
        const centre = (cursor + height / 2) / plotHeight;
        cursor += height + SANKEY_PAD;
        return Math.max(centre, SANKEY_MIN_Y);
    });
}

const flow = reclassificationCharts.flow;
if (flow.values.length) {
    const columnSize = flow.node_totals.length / 2;
    const grandTotal = flow.values.reduce((running, value) => running + value, 0);
    const plotHeight = SANKEY_HEIGHT - SANKEY_MARGIN * 2;
    // a significance with nothing flowing through it is unlinked, so plotly heaps it in with the
    // first column - leave room for that or the bands it does draw run past the bottom
    const unlinked = flow.node_totals.filter(total => !total).length;
    const usable = plotHeight - SANKEY_PAD * (columnSize + unlinked - 1);
    const centres = sankeyColumnCentres(flow.node_totals.slice(0, columnSize), grandTotal, usable, plotHeight)
        .concat(sankeyColumnCentres(flow.node_totals.slice(columnSize), grandTotal, usable, plotHeight));
    Plotly.newPlot('significance-flow', [{
        type: 'sankey',
        orientation: 'h',
        arrangement: 'fixed',
        node: {
            label: flow.labels,
            color: flow.node_colours,
            x: flow.node_totals.map((total, index) => index < columnSize ? 0.001 : 0.999),
            y: centres,
            pad: SANKEY_PAD,
            thickness: 16,
            line: {color: "#fcfcfb", width: 2}
        },
        link: {
            source: flow.sources,
            target: flow.targets,
            value: flow.values,
            color: flow.link_colours
        }
    }], chartLayout({
        height: SANKEY_HEIGHT,
        margin: {t: SANKEY_MARGIN, r: SANKEY_MARGIN, b: SANKEY_MARGIN, l: SANKEY_MARGIN}
    }), plotConfig);
}

const survival = reclassificationCharts.survival;
if (survival.intervals.length) {
    const traces = survival.curves.map(curve => ({
        x: survival.intervals,
        y: curve.survival.map(remaining => remaining * 100),
        name: curve.label,
        type: 'scatter',
        mode: 'lines+markers',
        line: {color: curve.colour, width: 2, shape: 'hv'},
        marker: {color: curve.colour, size: 6},
        hovertemplate: '%{y:.1f}% ' + curve.label.toLowerCase() + ' after %{x} years<extra></extra>'
    }));
    Plotly.newPlot('survival', traces, chartLayout({
        hovermode: 'x unified',
        legend: {traceorder: 'normal'},
        xaxis: {title: {text: 'Years since the cohort was taken'}, tickmode: 'linear', dtick: 0.5,
                gridcolor: GRID_COLOUR, zeroline: false},
        yaxis: {title: {text: '% of the cohort still waiting'}, range: [0, 100],
                gridcolor: GRID_COLOUR, zeroline: false}
    }), plotConfig);
}

const activity = reclassificationCharts.activity;
if (activity.years.length) {
    // the two rates share a left axis, the ratio between them needs its own
    const traces = [
        {
            x: activity.years, y: activity.reviewed, name: 'Records re-evaluated', type: 'bar',
            marker: {color: activity.review_colour, opacity: 0.35},
            hovertemplate: '%{y} records re-evaluated in %{x}<extra></extra>', yaxis: 'y3'
        },
        {
            x: activity.years, y: activity.reclassified, name: 'Records reclassified', type: 'bar',
            marker: {color: activity.change_colour, opacity: 0.35},
            hovertemplate: '%{y} records reclassified in %{x}<extra></extra>', yaxis: 'y3'
        },
        {
            x: activity.years, y: activity.reviewed_percents, name: 'Re-evaluation rate',
            type: 'scatter', mode: 'lines+markers', line: {color: activity.review_colour, width: 2},
            customdata: activity.populations,
            hovertemplate: '%{y}% of %{customdata} records re-evaluated<extra></extra>'
        },
        {
            x: activity.years, y: activity.reclassified_percents, name: 'Reclassification rate',
            type: 'scatter', mode: 'lines+markers',
            line: {color: activity.change_colour, width: 2, dash: 'dot'},
            customdata: activity.populations,
            hovertemplate: '%{y}% of %{customdata} records reclassified<extra></extra>'
        },
        {
            x: activity.years, y: activity.reviews_per_change, name: 'Records reviewed per change',
            type: 'scatter', mode: 'lines', line: {color: AXIS_COLOUR, width: 1, dash: 'dash'},
            yaxis: 'y2', hovertemplate: '%{y} reviewed per change<extra></extra>'
        }
    ];
    Plotly.newPlot('curation-activity', traces, chartLayout({
        height: 480,
        barmode: 'group',
        hovermode: 'x unified',
        legend: {orientation: 'h', y: -0.2},
        margin: {t: 20, r: 70, b: 80, l: 70},
        xaxis: {title: {text: 'Year'}, tickmode: 'linear', gridcolor: GRID_COLOUR, zeroline: false},
        yaxis: {title: {text: '% of the catalogue'}, gridcolor: GRID_COLOUR, zeroline: false},
        yaxis2: {title: {text: 'Reviewed per change'}, overlaying: 'y', side: 'right', showgrid: false,
                 zeroline: false, rangemode: 'tozero'},
        yaxis3: {overlaying: 'y', side: 'right', showgrid: false, zeroline: false, visible: false,
                 rangemode: 'tozero'}
    }), plotConfig);

    Plotly.newPlot('serial-activity', [
        {
            x: activity.serial_labels, y: activity.serial_reviewed, name: 'Re-evaluated', type: 'bar',
            marker: {color: activity.review_colour},
            hovertemplate: '%{y} records re-evaluated in %{x}<extra></extra>'
        },
        {
            x: activity.serial_labels, y: activity.serial_reclassified, name: 'Reclassified',
            type: 'bar', marker: {color: activity.change_colour},
            hovertemplate: '%{y} records reclassified in %{x}<extra></extra>'
        }
    ], chartLayout({
        barmode: 'group',
        xaxis: {title: {text: 'Distinct years the record was touched in'}, gridcolor: GRID_COLOUR, zeroline: false},
        yaxis: {title: {text: 'Records'}, type: 'log', gridcolor: GRID_COLOUR, zeroline: false}
    }), plotConfig);
}

const geneBurden = reclassificationCharts.gene_burden;
if (geneBurden.symbols.length) {
    Plotly.newPlot('gene-burden', [{
        x: geneBurden.vus_counts.slice().reverse(),
        y: geneBurden.symbols.slice().reverse(),
        type: 'bar',
        orientation: 'h',
        marker: {color: geneBurden.colour},
        customdata: geneBurden.percents.slice().reverse(),
        hovertemplate: '%{y}: %{x} VUS records, %{customdata}% of the gene<extra></extra>'
    }], chartLayout({
        height: 700,
        margin: {t: 20, r: 20, b: 60, l: 120},
        xaxis: {title: {text: 'Records currently at VUS'}, gridcolor: GRID_COLOUR, zeroline: false},
        yaxis: {gridcolor: GRID_COLOUR, zeroline: false, automargin: true}
    }), plotConfig);
}

const labs = reclassificationCharts.labs;
if (labs.labs.length) {
    $('#lab-league').dataTable({
        paginate: false,
        searching: false,
        info: false,
        order: [[2, 'desc']]
    });

    Plotly.newPlot('lab-rates', [
        {
            x: labs.labs, y: labs.reviewed_percents, name: 'Re-evaluation rate', type: 'bar',
            marker: {color: labs.review_colour}, customdata: labs.held,
            hovertemplate: '%{x}: %{y}% of %{customdata} records re-evaluated<extra></extra>'
        },
        {
            x: labs.labs, y: labs.reclassified_percents, name: 'Reclassification rate', type: 'bar',
            marker: {color: labs.change_colour}, customdata: labs.held,
            hovertemplate: '%{x}: %{y}% of %{customdata} records reclassified<extra></extra>'
        }
    ], chartLayout({
        height: 480,
        barmode: 'group',
        margin: {t: 20, r: 20, b: 160, l: 70},
        xaxis: {gridcolor: GRID_COLOUR, zeroline: false, automargin: true, tickangle: -40},
        yaxis: {title: {text: '% of records held'}, gridcolor: GRID_COLOUR, zeroline: false}
    }), plotConfig);
}

const points = reclassificationCharts.points;
if (points.labels.length) {
    const traces = points.labels.map((label, index) => ({
        x: points.deltas[index],
        name: label,
        type: 'box',
        orientation: 'h',
        boxpoints: false,
        marker: {color: points.colours[index]},
        line: {color: points.colours[index]},
        hovertemplate: label + ': %{x} points<extra></extra>'
    }));
    Plotly.newPlot('points-transitions', traces, chartLayout({
        height: Math.max(320, 44 * traces.length + 120),
        showlegend: false,
        margin: {t: 20, r: 20, b: 60, l: 120},
        xaxis: {title: {text: 'ACMG points travelled, positive towards pathogenic'},
                gridcolor: GRID_COLOUR, zeroline: true, zerolinecolor: AXIS_COLOUR},
        yaxis: {gridcolor: GRID_COLOUR, zeroline: false, automargin: true}
    }), plotConfig);
}

/** Unapplied and weakened to the left, applied, strengthened and plain changes to the right, keys
 *  running down in the order the view put them in. The keys that folded into the other row open up
 *  when it's clicked. */
function plotEvidenceMovement(selector, evidence) {
    if (!evidence.collapsed.labels.length) {
        return;
    }
    const chart = document.getElementById(selector);
    const toggle = document.getElementById(selector + '-toggle');
    let expanded = false;

    function draw() {
        const series = expanded ? evidence.expanded : evidence.collapsed;
        const labels = series.labels.slice().reverse();
        const bar = (values, name, colour, verb) => ({
            x: values.slice().reverse(), y: labels, name: name, type: 'bar', orientation: 'h',
            marker: {color: colour},
            hovertemplate: '%{y} ' + verb + ' in %{x} reclassifications<extra></extra>'
        });
        const traces = [
            bar(series.unapplied, 'Unapplied', evidence.unapplied_colour, 'unapplied'),
            bar(series.weakened, 'Weakened', evidence.weakened_colour, 'weakened'),
            bar(series.applied, 'Applied', evidence.applied_colour, 'applied'),
            bar(series.strengthened, 'Strengthened', evidence.strengthened_colour, 'strengthened'),
            bar(series.changed, 'Changed', evidence.changed_colour, 'changed')
        ];
        Plotly.newPlot(selector, traces, chartLayout({
            height: Math.max(320, 26 * labels.length + 120),
            barmode: 'relative',
            margin: {t: 20, r: 20, b: 60, l: 220},
            legend: {traceorder: 'normal'},
            xaxis: {title: {text: 'Unapplied or weakened \u2190  reclassifications  \u2192 applied, strengthened or changed'},
                    tickformat: 'd', gridcolor: GRID_COLOUR, zeroline: true, zerolinecolor: AXIS_COLOUR},
            yaxis: {gridcolor: GRID_COLOUR, zeroline: false, automargin: true}
        }), plotConfig);

        if (toggle && evidence.folded_count) {
            toggle.textContent = expanded ? 'Gather the smaller keys back into other'
                                          : 'Show all ' + evidence.folded_count + ' keys in other';
        }

        const keys = series.keys.slice().reverse();
        chart.on('plotly_click', event => {
            if (!expanded && keys[event.points[0].pointNumber] === evidence.other_key) {
                expanded = true;
                draw();
            }
        });
    }

    if (toggle && evidence.folded_count) {
        toggle.classList.remove('d-none');
        toggle.addEventListener('click', event => {
            event.preventDefault();
            expanded = !expanded;
            draw();
        });
    }
    draw();
}

plotEvidenceMovement('evidence-towards-pathogenic', reclassificationCharts.evidence_towards_pathogenic);
plotEvidenceMovement('evidence-towards-benign', reclassificationCharts.evidence_towards_benign);
