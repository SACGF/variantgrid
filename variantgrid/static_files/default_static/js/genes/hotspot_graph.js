// @ts-check
// genes/templates/genes/hotspot_graph.html
/* Loaded by ajax, possibly several graphs per page - each one calls this with its own uuid and values
   hotspotData: molecularConsequenceColors, numCodons, domains, variantData, transcriptUrls, title, yTitle,
                transcriptVersionId */
function setupHotspotGraph(uuid, hotspotData) {
    const GNOMAD_PERCENT = [0.1, 1, 5, 100];
    const molecularConsequenceColors = hotspotData.molecularConsequenceColors;
    const numCodons = hotspotData.numCodons;
    const domains = hotspotData.domains;
    const variantData = hotspotData.variantData;
    const transcriptUrls = hotspotData.transcriptUrls;
    const hotSpotId = "hotspot-graph-" + uuid;
    const hotSpotDiv = $("#" + hotSpotId);

    function drawTranscriptModel(selector, gnomADMaxPercent) {
        const DEFAULT_COLORS = [
            '#d62728',  // brick red
            '#2ca02c',  // cooked asparagus green
            '#1f77b4',  // muted blue
            '#ff7f0e',  // safety orange
            '#e377c2',  // raspberry yogurt pink
            '#bcbd22',  // curry yellow-green
            '#17becf',  // blue-teal
            '#9467bd',  // muted purple
            '#8c564b'   // chestnut brown
        ];

        let titleText = escapeHtml(hotspotData.title);
        if (gnomADMaxPercent != 100) {
            titleText += "\ngnomAD max " + gnomADMaxPercent + "%";
        }
        const gnomadAFMax = gnomADMaxPercent / 100.0;
        const consequences = Object.keys(molecularConsequenceColors).sort();
        let yMax = 5; // show minimum so gene diagram doesn't get too fat
        const variantsByConsequence = {};
        for (let i = 0; i < consequences.length; ++i) {
            variantsByConsequence[consequences[i]] = {x: [], y: [], text: []};
        }

        for (let i = 0; i < variantData.length; ++i) {
            const vd = variantData[i];
            const text = vd[0];
            const x = vd[1];
            const consequence = vd[2];
            const gnomadAf = vd[3];
            const numSamples = vd[4];
            const vbc = variantsByConsequence[consequence];

            if (gnomadAf && gnomadAf > gnomadAFMax) {
                continue;
            }

            vbc["x"].push(x);
            vbc["y"].push(numSamples);
            vbc["text"].push(text);
        }

        const BAR_WIDTH = 1.5;
        const data = [];
        // variantData is unique per allele, but there may be different alleles that have the same protein position
        // and consequence on the graph. The bar is by default stacked so this will push it up but we need to put
        // the lollypop on top after working out how high it is.
        const mergedXData = {};

        for (let i = 0; i < consequences.length; ++i) {
            const consequence = consequences[i];
            const v = variantsByConsequence[consequence];
            const color = molecularConsequenceColors[consequence];

            const barTrace = {
                x: v.x,
                y: v.y,
                width: Array.from({length: v.y.length}).map(x => BAR_WIDTH),
                hoverinfo: 'skip',
                name: consequence,
                type: 'bar',
                marker: {
                    color: color
                }
            };
            data.push(barTrace);

            for (let j = 0; j < v.x.length; ++j) {
                const x = v.x[j];
                const mergedData = mergedXData[x] || {count: 0, text: {}, last_consequence: null};
                mergedData.count += v.y[j];
                mergedData.text[v.text[j]] = 1;
                mergedData.last_consequence = consequence;
                mergedXData[x] = mergedData;
            }
        }

        // Lollypop head
        const consequenceLollypops = {};
        for (const [ x, mergedData ] of Object.entries(mergedXData)) {
            const lpData = consequenceLollypops[mergedData.last_consequence] || {x: [], y: [], text: []};
            lpData.x.push(x);
            lpData.y.push(mergedData.count);
            lpData.text.push(Object.keys(mergedData.text).join(" "));
            consequenceLollypops[mergedData.last_consequence] = lpData;
        }

        for (const [consequence, lpData] of Object.entries(consequenceLollypops)) {
            const color = molecularConsequenceColors[consequence];
            const scatterTrace = {
                x: lpData.x,
                y: lpData.y,
                text: lpData.text,
                name: consequence,
                showlegend: false,
                mode: 'markers',
                type: 'scatter',
                marker: {
                    color: color,
                    size: BAR_WIDTH * 5
                }
            };
            data.push(scatterTrace);
            yMax = Math.max(yMax, ...lpData.y);
        }


        const geneThickness = yMax / 5;
        const domainThickness = geneThickness; // by default no extra thickness
        const geneYTop = 0;
        const geneYBottom = geneYTop - geneThickness;
        // make domain lie entirely below, as it was sometimes obscuring classifications on top
        const domainYTop = geneYTop;
        const domainYBottom = geneYBottom - (domainThickness - geneThickness);

        // Need to give domains consistent colors
        let colorIndex = 0;
        const domainColors = {};

        const shapes = [
            {
                type: "rect",
                x0: 0,
                y0: geneYTop,
                x1: numCodons + 1, // Ends "on" that AA so need to go to next one to cover it.
                y1: geneYBottom,
                fillcolor: "#b0b0b0",
            }
        ];

        /* Fonts are specified in font-size units so we have to scale them ourselves
           and abbreviate or not show if they overflow the protein domain shape */
        const MIN_CHARS = 2;
        const graphSize = $("#" + selector).width();
        const fontScale = graphSize / (numCodons + 1);
        const fontSize = 10;
        const fontWidth = fontSize / fontScale;

        const annotations = [];
        for (let i = 0; i < domains.length; ++i) {
            const d = domains[i];
            const domainName = d[0];
            const domainStart = d[2];
            const domainEnd = d[3] + 1; // Ends "on" that AA so need to go to next one to cover it.
            let domainColor = domainColors[domainName];
            if (typeof (domainColor) == 'undefined') {
                domainColor = DEFAULT_COLORS[colorIndex];
                colorIndex++;
                domainColors[domainName] = domainColor;
            }
            const domainShape = {
                type: "rect",
                x0: domainStart,
                y0: domainYTop,
                x1: domainEnd,
                y1: domainYBottom,
                fillcolor: domainColor,
            };
            shapes.push(domainShape);
            const domainWidth = domainEnd - domainStart;
            const domainMaxChars = domainWidth / fontWidth;
            let domainText = "";
            if (domainMaxChars >= MIN_CHARS) {
                domainText = domainName;
                if (domainText.length > domainMaxChars) {
                    domainText = domainText.substring(0, domainMaxChars - 1) + ".";
                }
            }
            const annotation = {
                showarrow: false,
                x: domainStart + (domainEnd - domainStart) / 2,
                y: domainYBottom + domainThickness / 2,
                text: "<b>" + domainText + "</b>",
                font: {
                    color: "white",
                    size: fontSize
                },
                xanchor: "center",
            };
            annotations.push(annotation);
        }
        const layout = {
            barmode: 'stack',
            title: {
                text: titleText
            },
            showlegend: true,
            xaxis: {
                title: {text: "Amino acid"},
                range: [0, numCodons + 1],
                showgrid: false
            },
            yaxis: {
                title: {text: escapeHtml(hotspotData.yTitle)},
                range: [domainYBottom, yMax + 2],
                nticks: 5,
                showgrid: false
            },
            shapes: shapes,
            annotations: annotations,
        };

        const config = {
            showLink: true,
            plotlyServerURL: "https://chart-studio.plotly.com"
        };

        Plotly.newPlot(selector, data, layout, config);

        const myPlot = document.getElementById(selector);
        myPlot.on('plotly_click', function (data) {
            const hotspot_graph_click_func = hotSpotDiv.attr("hotspot_graph_click_func");
            if (typeof (hotspot_graph_click_func) != 'undefined') {
                const barClicked = data.points[0];
                const i = barClicked.data.x.indexOf(String(barClicked.x));
                const text = "Hotspot click " + barClicked.data.text[i];
                const fn = window[hotspot_graph_click_func];
                if (typeof fn === 'function') {
                    fn(String(hotspotData.transcriptVersionId), text, barClicked.x);
                }
            }
        });
    }

    function clickTranscript() {
        const accession = $(this).attr("accession");
        const url = transcriptUrls[accession];
        const container = hotSpotDiv.parent();
        const hotspot_graph_click_func = hotSpotDiv.attr("hotspot_graph_click_func");
        container.empty();
        container.load(url, function() {
            // put click handler back
            $(".hotspot-graph", this).attr("hotspot_graph_click_func", hotspot_graph_click_func);
        });
    }

    $(document).ready(function () {
        const slider = $("#hotspot-graph-af-slider-" + uuid);
        function drawMyTranscriptModel() {
            let gnomADMaxPercent = 100;
            if (slider.length) {
                gnomADMaxPercent = GNOMAD_PERCENT[Number(slider.val())];
            }
            drawTranscriptModel(hotSpotId, gnomADMaxPercent);
        }

        const gnomadMax = GNOMAD_PERCENT.length - 1;
        slider.attr({min: 0, max: gnomadMax, step: 1}).val(gnomadMax).on("change", drawMyTranscriptModel);
        // A label per step, positioned along the track
        for (let i = 0; i <= gnomadMax; i++) {
            const el = $('<label>' + GNOMAD_PERCENT[i] + '%</label>').css('left', (i / gnomadMax * 100) + '%');
            slider.parent().append(el);
        }

        drawMyTranscriptModel();

        $("a.hotspot-load-transcript-link", "#hotspot-transcripts-" + uuid).on('click', clickTranscript);
    });
}
