/**
 * Variant details "Samples" section.
 *
 * The server renders a shell, this fills it in from the variant_sample_genotypes API. Requests are
 * per variant (not per build) - builds sharing a contig resolve to the same variant, and each row
 * carries its own sample's build.
 */
const VariantSampleInformation = (function () {
    // Ordered as the zygosity filter checkboxes are drawn
    const ZYGOSITIES = [
        {code: 'R', label: 'REF'},
        {code: 'E', label: 'HET'},
        {code: 'O', label: 'HOM_ALT'},
        {code: 'U', label: 'Unknown Zygosity'},
    ];
    const ZYGOSITY_LABELS = Object.fromEntries(ZYGOSITIES.map(z => [z.code, z.label]));
    const DEFAULT_VISIBLE_ZYGOSITIES = ['E', 'O'];
    // Legacy data (missing values in FreeBayes before PythonKnownVariantsImporter v12) has
    // -2147483648 for empty values @see https://github.com/SACGF/variantgrid/issues/59
    const MISSING_FORMAT_VALUES = new Set([-1, '-1', -2147483648, '-2147483648', null, undefined]);
    // Matches the server rendered patient graphs in patients/templatetags/patient_graph_tags.py
    const HPO_COLOR = "#FF8C00";
    const OMIM_COLOR = "#22529e";
    const MONDO_COLOR = "#99CD83";
    const MAX_GRAPH_TERMS = 20;
    const GRAPH_WIDTH = 512;
    const GRAPH_HEIGHT = 384;

    let config = null;
    let dataTable = null;
    let allRows = [];
    let zygosityCounts = {};
    let checkedZygosities = null;  // null until we've seen a response and know what to default to

    function genotypesUrl(variantId) {
        return config.genotypesUrlForVariant.replace(/0$/, variantId);
    }

    function ontologyTermUrl(termType, term) {
        return Urls.ontology_term_text(termType, term);
    }

    function renderSampleLink(data, type, row) {
        if (type !== 'display') {
            return data;
        }
        return $('<a>', {href: row.sample_url, text: data}).prop('outerHTML');
    }

    function renderZygosity(data, type) {
        if (type === 'sort') {
            return data;
        }
        return ZYGOSITY_LABELS[data] || data;
    }

    function filterText(row) {
        const parts = [];
        if (row.filters) {
            parts.push(row.filters);
        }
        if (row.sample_filters) {
            parts.push("FT: " + row.sample_filters);
        }
        return parts.join(' ');
    }

    function filterSpan(value, text) {
        const cssClass = value === 'PASS' ? 'text-muted' : 'text-danger font-weight-bold';
        return $('<span>', {class: cssClass, text: text}).prop('outerHTML');
    }

    function renderFilters(data, type, row) {
        if (type !== 'display') {
            return filterText(row);
        }
        const spans = [];
        if (row.filters) {
            spans.push(filterSpan(row.filters, row.filters));
        }
        if (row.sample_filters) {
            spans.push(filterSpan(row.sample_filters, "FT: " + row.sample_filters));
        }
        return spans.join(' ');
    }

    function ontologyRenderer(termType) {
        return function (data, type) {
            if (type !== 'display') {
                return data ? data.split('|').join(', ') : '';
            }
            return getOntologyTermLinks(termType, data, ontologyTermUrl);
        };
    }

    function renderFormatField(data) {
        return MISSING_FORMAT_VALUES.has(data) ? '.' : data;
    }

    function columns() {
        return [
            {data: 'sample_name', title: 'Sample Name', render: renderSampleLink},
            {data: 'zygosity', title: 'Zygosity', render: renderZygosity},
            {data: 'filters', title: 'FILTER', render: renderFilters},
            {data: 'enrichment_kit', title: 'EnrichmentKit', defaultContent: ''},
            {data: 'patient_hpo', title: 'Patient HPO', className: 'no-word-wrap', render: ontologyRenderer('HP')},
            {data: 'patient_omim', title: 'Patient OMIM', className: 'no-word-wrap', render: ontologyRenderer('OMIM')},
            {data: 'patient_mondo', title: 'Patient MONDO', className: 'no-word-wrap', render: ontologyRenderer('MONDO')},
            {data: 'vcf', title: 'VCF', defaultContent: ''},
            {data: 'allele_frequency', title: 'Allele Frequency', render: renderFormatField},
            {data: 'allele_depth', title: 'Allele Depth', render: renderFormatField},
            {data: 'read_depth', title: 'Read Depth', render: renderFormatField},
            {data: 'phred_likelihood', title: 'Phred Likelihood', render: renderFormatField},
        ];
    }

    function csvFilename() {
        const slug = (config.gHgvs || config.variantLabel).replace(/[^a-zA-Z0-9]+/g, '_').replace(/^_|_$/g, '');
        const dateStr = new Date().toISOString().slice(0, 10).replace(/-/g, '');
        return 'variant_genotypes_' + slug + '_' + dateStr;
    }

    function createDataTable() {
        dataTable = $("#genotype-grid").DataTable({
            data: [],
            columns: columns(),
            pageLength: 20,
            order: [],
            dom: 'Bfrtip',
            buttons: [{
                extend: 'csv',
                text: 'CSV',
                exportOptions: {orthogonal: 'export'},
                filename: csvFilename,
            }],
        });

        $.fn.dataTable.ext.search.push(function (settings, searchData, dataIndex, rowData) {
            if (settings.nTable.id !== 'genotype-grid' || checkedZygosities === null) {
                return true;
            }
            return checkedZygosities.has(rowData.zygosity);
        });
    }

    function drawZygosityFilters() {
        if (checkedZygosities === null) {
            // Show variants by default, but don't hide everything if that's all there is
            const hasDefaults = DEFAULT_VISIBLE_ZYGOSITIES.some(code => zygosityCounts[code]);
            const visible = hasDefaults ? DEFAULT_VISIBLE_ZYGOSITIES : ZYGOSITIES.map(z => z.code);
            checkedZygosities = new Set(visible);
        }

        const container = $("#zygosity-filters").empty();
        for (const zygosity of ZYGOSITIES) {
            const count = zygosityCounts[zygosity.code];
            if (!count) {
                continue;
            }
            const checkbox = $('<input>', {
                type: 'checkbox',
                class: 'genotype-filter',
                value: zygosity.code,
            }).prop('checked', checkedZygosities.has(zygosity.code));
            container.append($('<label>', {text: `${zygosity.label} (${count}):`}), checkbox,
                $('<span>', {class: 'separator', text: '|'}));
        }

        container.find(".genotype-filter").change(function () {
            checkedZygosities = new Set(container.find(".genotype-filter:checked")
                .map(function () { return $(this).val(); }).get());
            dataTable.draw();
        });
    }

    function drawLocusCounts(locusCounts) {
        const showUnknown = locusCounts.some(row => row["Unknown"]);
        const headers = ['Allele', 'Total', 'Hom Ref', 'Het', 'Hom'];
        const fields = ['total', 'REF', 'HET', 'HOM_ALT'];
        if (showUnknown) {
            headers.push('Unknown');
            fields.push('Unknown');
        }

        const table = $('<table>', {class: 'table'});
        const headerRow = $('<tr>');
        headers.forEach(h => headerRow.append($('<th>', {text: h})));
        table.append($('<thead>').append(headerRow));

        const tbody = $('<tbody>');
        for (const row of locusCounts) {
            const link = $('<a>', {class: 'hover-link', href: row.url, text: row.variant});
            const description = row.description ? ` (${row.description})` : '';
            const tr = $('<tr>').append($('<td>').append(link).append(document.createTextNode(description)));
            fields.forEach(f => tr.append($('<td>', {text: row[f]})));
            tbody.append(tr);
        }
        table.append(tbody);
        $("#locus-counts").empty().append(table);
    }

    function drawMultiallelic(multiallelic) {
        const container = $("#other-loci-variants-by-multiallelic").empty();
        if (!multiallelic.length) {
            return;
        }
        container.append($('<h4>', {text: 'Variants from multi-allelic VCF records'}));
        container.append($('<p>', {
            text: 'Multi-allelic VCF records are split then separately normalised. This may cause records ' +
                'that were on a line together in a VCF to now have different loci.'
        }));
        for (const entry of multiallelic) {
            container.append($('<b>', {text: entry.multiallelic}));
            const list = $('<ul>');
            for (const variant of entry.variants) {
                list.append($('<li>').append($('<a>', {href: variant.url, text: variant.label})));
            }
            container.append(list);
        }
    }

    function drawSummary(data) {
        const builds = data.genome_builds.join(', ');
        const samples = data.samples;
        let summary = `Searching ${samples.total} samples for ${builds}`;
        summary += samples.visible < samples.total ? ` (you can see ${samples.visible}).` : '.';

        const container = $("#samples-summary").empty();
        if (!data.observations.total) {
            container.append($('<p>', {text: `No results found (searched ${samples.total} samples for ${builds}).`}));
            return;
        }
        container.append($('<p>', {text: summary}));

        if (data.observations.invisible) {
            const hidden = $('<div>').append(
                $('<img>', {src: config.unknownSampleIconUrl}),
                $('<b>', {text: `You cannot see ${data.observations.invisible} samples`}),
                $('<p>', {
                    text: `You only have permission to view ${data.observations.visible} out of ` +
                        `${data.observations.total} total observations in the database.`
                }));
            container.append(hidden);
        }
    }

    function plotAlleleFrequencies() {
        const x = [];
        const alleleFrequencies = [];
        const labels = [];
        for (const row of allRows) {
            x.push((Math.random() - 0.5) * 0.1);  // jitter
            alleleFrequencies.push(row.allele_frequency);
            let sampleDescription = `Sample: (${row.sample})`;
            if (row.sample_name) {
                sampleDescription += " " + row.sample_name;
            }
            sampleDescription += " DP=" + row.read_depth;
            labels.push(sampleDescription);
        }

        const scatterData = {
            name: "Allele Frequency",
            type: 'scatter',
            mode: 'markers',
            x: x,
            y: alleleFrequencies,
            text: labels,
        };
        const layout = defaultLayout(config.variantLabel + " Allele Frequency");
        layout.xaxis = Object.assign(layout.xaxis || {}, {zeroline: false, showgrid: false});
        layout.yaxis = Object.assign(layout.yaxis || {}, {range: [0, 1.05], showgrid: false, zeroline: false});
        Plotly.newPlot('sample-allele-frequency-scatter', [scatterData], layout);

        // Need to clamp 100 so that histo fits in 10 bins
        const clampedAf = alleleFrequencies.map(af => af >= 1 ? 0.99 : af);
        const histogramData = {
            name: "Allele Frequency",
            type: 'histogram',
            x: clampedAf,
            autobinx: false,
            histnorm: "count",
            xbins: {start: 0, size: .05, end: 1},
        };
        const histogramLayout = defaultLayout(config.variantLabel + " Allele Frequency Histogram");
        histogramLayout.xaxis = {autotick: true, range: [0, 1]};
        histogramLayout.yaxis = {showline: false, zeroline: false, showticklabels: false};
        Plotly.newPlot('sample-allele-frequency-histogram', [histogramData], histogramLayout);
    }

    function plotPatientTermCounts(selector, title, color, rowField) {
        // The graph counts patients, while the grid is per sample - a patient can have several
        const patientsByTerm = {};
        for (const row of allRows) {
            if (!row.patient || !row[rowField]) {
                continue;
            }
            for (const term of row[rowField].split('|')) {
                (patientsByTerm[term] = patientsByTerm[term] || new Set()).add(row.patient);
            }
        }

        const counts = Object.entries(patientsByTerm).map(([term, patients]) => [term, patients.size]);
        if (!counts.length) {
            return false;
        }
        counts.sort((a, b) => b[1] - a[1]);
        const top = counts.slice(0, MAX_GRAPH_TERMS).reverse();  // Plotly draws h-bars bottom up
        plotHBarArrays(selector, title + "<br>Patient Count", top.map(c => c[1]), top.map(c => c[0]),
            GRAPH_WIDTH, GRAPH_HEIGHT, color, {l: 200});
        return true;
    }

    function drawPhenotypeGraphs() {
        const hpo = plotPatientTermCounts('sample-hpo-graph', "Human Phenotype Ontology", HPO_COLOR, 'patient_hpo');
        const omim = plotPatientTermCounts('sample-omim-graph', "OMIM", OMIM_COLOR, 'patient_omim');
        const mondo = plotPatientTermCounts('sample-mondo-graph', "MONDO", MONDO_COLOR, 'patient_mondo');
        $("#variant-sample-phenotype-graphs").toggle(hpo || omim || mondo);
    }

    function addResponse(data) {
        config.gHgvs = config.gHgvs || data.g_hgvs;
        drawSummary(data);
        drawLocusCounts(data.locus_counts);
        drawMultiallelic(data.multiallelic);

        for (const [zygosity, count] of Object.entries(data.zygosity_counts)) {
            zygosityCounts[zygosity] = (zygosityCounts[zygosity] || 0) + count;
        }

        if (data.rows.length) {
            allRows = allRows.concat(data.rows);
            $("#genotype-grid-container").show();
            drawZygosityFilters();
            dataTable.rows.add(data.rows).draw();
            $("#allele-frequency-graphs-container").show();
            plotAlleleFrequencies();
            drawPhenotypeGraphs();
        }
    }

    function init(cfg) {
        config = cfg;
        createDataTable();

        // getOntologyTermLinks makes these as the grid renders, so delegate from the document
        $(document).on('click', '.ontology-terms-container .collapsed-term', expandCollapsedOntologyTerm);

        $.getJSON(genotypesUrl(config.variantId), addResponse).fail(function () {
            $("#samples-summary").html($('<p>', {class: 'text-danger', text: 'Error retrieving samples'}));
        }).always(function () {
            $("#samples-loading").hide();
        });
    }

    return {
        init: init,
    };
})();
