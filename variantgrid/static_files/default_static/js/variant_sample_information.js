/**
 * Variant details "Samples" section.
 *
 * The server renders a shell, this fills it in from the variant_sample_genotypes API. The same mutation
 * is a separate variant in each genome build, so the server hands us one variant per build (via the
 * allele) and we request each - responses accumulate into the one grid, keyed by the row's own build.
 *
 * The allele may be created (and linked to the other builds' variants) by the variant details page while
 * this is loading, so we also listen for ALLELE_VARIANTS_EVENT and pick up any variants that turned up.
 */
// var not const - the analysis' variant details tabs load this script again for each variant
var VariantSampleInformation = (function () {
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
    const AF_COLOR = "#1f77b4";
    // Graphs share a line, sizing themselves to their flex column - see .graph-row
    const RESPONSIVE = {responsive: true};
    const MAX_GRAPH_TERMS = 20;
    // The expandable row detail - patient/specimen/extraction identifiers that would make an
    // already wide grid wider (private#2837). The server leaves out anything it has no value for
    const DETAIL_FIELDS = [
        {key: 'patient_code', label: 'Patient', urlKey: 'patient_url'},
        {key: 'specimen', label: 'Specimen', urlKey: 'specimen_url'},
        {key: 'specimen_tissue', label: 'Tissue'},
        {key: 'specimen_tissue_status', label: 'Tissue Status'},
        {key: 'specimen_collection_date', label: 'Collected'},
        {key: 'extraction', label: 'Extraction', urlKey: 'extraction_url'},
        {key: 'nucleic_acid', label: 'Nucleic Acid'},
    ];
    const GRAPH_HEIGHT = 384;
    const ALLELE_VARIANTS_EVENT = "allele-variants-loaded";
    // Namespaced so re-loading the section (another variant details tab) replaces its document handlers
    const EVENT_NAMESPACE = ".variantSampleInformation";

    let config = null;
    let dataTable = null;
    const responses = new Map();  // variant_id -> payload, replaced when that variant is re-loaded uncapped
    let requestedVariantIds = new Set();
    let pendingRequests = 0;
    let allRows = [];
    let checkedZygosities = null;  // null until we've seen a response and know what to default to
    let genomeBuildChecked = {};  // build -> bool, a build starts checked when its first rows arrive
    let graphedFilter = null;  // Redrawing 10k plotly points on every page change is a waste

    function genotypesUrl(variantId, limit) {
        const url = config.genotypesUrlForVariant.replace(/0$/, variantId);
        return limit === undefined ? url : url + "?limit=" + limit;
    }

    function loadVariant(variantId, limit) {
        return $.getJSON(genotypesUrl(variantId, limit)).done(addResponse);
    }

    function showError(message) {
        $("#samples-messages").empty().append(createMessage("error", message));
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
        // PASS is the boring case so let it recede, anything else reads as normal body text
        const cssClass = value === 'PASS' ? 'text-muted' : '';
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

    function classificationPill(classification) {
        // Same shape as the server rendered clinical_significance_values tag - CSS is in global.scss
        const pill = $('<span>', {
            class: 'd-inline-block c-pill ' + classification.css_class,
            title: classification.title,
            text: classification.label,
        });
        if (classification.pending_from) {
            pill.append($('<i>', {
                class: 'fa-solid fa-clock ml-1',
                title: `This is a pending change, the classification is currently ${classification.pending_from}`,
            }));
        }
        return $('<a>', {href: classification.url, class: 'mr-1', title: classification.lab})
            .append(pill).prop('outerHTML');
    }

    function renderClassifications(data, type) {
        if (!data || !data.length) {
            return '';
        }
        if (type !== 'display') {
            return data.map(classification => classification.label).join(', ');
        }
        return data.map(classificationPill).join('');
    }

    function renderFormatField(data) {
        return MISSING_FORMAT_VALUES.has(data) ? '.' : data;
    }

    function hasDetails(row) {
        return DETAIL_FIELDS.some(field => row[field.key]);
    }

    function renderDetailsControl(data, type, row) {
        if (type !== 'display' || !hasDetails(row)) {
            return '';
        }
        return $('<i>', {class: 'fas fa-chevron-right', title: 'Patient / specimen details'}).prop('outerHTML');
    }

    function detailsList(row) {
        const list = $('<dl>', {class: 'row sample-details'});
        for (const field of DETAIL_FIELDS) {
            const value = row[field.key];
            if (!value) {
                continue;
            }
            const url = field.urlKey && row[field.urlKey];
            const dd = $('<dd>', {class: 'col-sm-9'});
            dd.append(url ? $('<a>', {href: url, text: value}) : document.createTextNode(value));
            list.append($('<dt>', {class: 'col-sm-3', text: field.label}), dd);
        }
        return list;
    }

    function toggleDetails() {
        const tr = $(this).closest('tr');
        const row = dataTable.row(tr);
        if (!row.data() || !hasDetails(row.data())) {
            return;
        }
        if (row.child.isShown()) {
            row.child.hide();
        } else {
            row.child(detailsList(row.data())).show();
        }
        $(this).find('i').toggleClass('fa-chevron-right fa-chevron-down');
    }

    function columns() {
        return [
            {data: null, name: 'details', title: '', className: 'details-control',
                orderable: false, searchable: false, render: renderDetailsControl},
            {data: 'sample_name', title: 'Sample Name', render: renderSampleLink},
            {data: 'genome_build', name: 'genome_build', title: 'Genome Build', visible: false},
            {data: 'zygosity', title: 'Zygosity', render: renderZygosity},
            // Hidden until a row actually has one - see addResponse
            {data: 'classifications', name: 'classifications', title: 'Classifications',
                visible: false, render: renderClassifications},
            {data: 'allele_frequency', title: 'Allele Frequency', render: renderFormatField},
            {data: 'allele_depth', title: 'Allele Depth', render: renderFormatField},
            {data: 'read_depth', title: 'Read Depth', render: renderFormatField},
            {data: 'phred_likelihood', title: 'Phred Likelihood', render: renderFormatField},
            {data: 'filters', title: 'FILTER', render: renderFilters},
            {data: 'enrichment_kit', title: 'EnrichmentKit', defaultContent: ''},
            {data: 'patient_hpo', title: 'Patient HPO', className: 'no-word-wrap', render: ontologyRenderer('HP')},
            {data: 'patient_omim', title: 'Patient OMIM', className: 'no-word-wrap', render: ontologyRenderer('OMIM')},
            {data: 'patient_mondo', title: 'Patient MONDO', className: 'no-word-wrap', render: ontologyRenderer('MONDO')},
            {data: 'vcf', title: 'VCF', defaultContent: ''},
            // Drawn in the row detail rather than taking up width, but still columns so that the
            // search box and the CSV reach them - that's how you spot rows sharing a patient
            {data: 'patient_code', title: 'Patient', visible: false, defaultContent: ''},
            {data: 'specimen', title: 'Specimen', visible: false, defaultContent: ''},
            {data: 'extraction', title: 'Extraction', visible: false, defaultContent: ''},
        ];
    }

    function csvFilename() {
        const slug = (config.gHgvs || config.variantLabel).replace(/[^a-zA-Z0-9]+/g, '_').replace(/^_|_$/g, '');
        const dateStr = new Date().toISOString().slice(0, 10).replace(/-/g, '');
        return 'variant_genotypes_' + slug + '_' + dateStr;
    }

    function truncatedVariantIds() {
        return [...responses.values()].filter(data => data.truncated).map(data => data.variant_id);
    }

    function exportCsv(e, dt, node, csvConfig) {
        const csvAction = $.fn.dataTable.ext.buttons.csvHtml5.action;
        const variantIds = truncatedVariantIds();
        if (!variantIds.length) {
            csvAction.call(this, e, dt, node, csvConfig);
            return;
        }

        // A download must not silently be a truncated view of the data
        const button = $(node);
        const buttonHtml = button.html();
        button.addClass('disabled').empty()
            .append($('<i>', {class: 'fas fa-spinner fa-spin mr-1'}), document.createTextNode('Loading all rows...'));
        $.when(...variantIds.map(variantId => loadVariant(variantId, 0)))
            .done(() => csvAction.call(this, e, dt, node, csvConfig))
            .fail(() => showError("Error loading all rows - CSV export cancelled"))
            .always(() => button.removeClass('disabled').html(buttonHtml));
    }

    function createDataTable() {
        dataTable = $("#genotype-grid").DataTable({
            data: [],
            columns: columns(),
            pageLength: 20,
            order: [],
            scrollX: true,  // Wide grid scrolls itself rather than the variant details page
            dom: 'Bfrtip',
            buttons: [{
                extend: 'csv',
                text: 'CSV',
                // Everything but the expand control, so the detail columns are in the download too
                exportOptions: {orthogonal: 'export', columns: ':not(.details-control)'},
                filename: csvFilename,
                action: exportCsv,
            }],
        });

        $("#genotype-grid tbody").on('click', 'td.details-control', toggleDetails);

        $.fn.dataTable.ext.search.push(function (settings, searchData, dataIndex, rowData) {
            if (settings.nTable.id !== 'genotype-grid') {
                return true;
            }
            if (genomeBuildChecked[rowData.genome_build] === false) {
                return false;
            }
            return checkedZygosities === null || checkedZygosities.has(rowData.zygosity);
        });

        // Graphs summarise what the grid is showing, so follow the zygosity/build/search filters
        dataTable.on('draw', drawGraphs);
    }

    /** The rows the grid is currently showing - across all pages, but after zygosity/build/search filtering */
    function visibleRows() {
        return dataTable.rows({search: 'applied'}).data().toArray();
    }

    function zygosityCounts() {
        const counts = {};
        for (const data of responses.values()) {
            for (const [zygosity, count] of Object.entries(data.zygosity_counts)) {
                counts[zygosity] = (counts[zygosity] || 0) + count;
            }
        }
        return counts;
    }

    function genomeBuildCounts() {
        const counts = {};
        for (const row of allRows) {
            counts[row.genome_build] = (counts[row.genome_build] || 0) + 1;
        }
        return counts;
    }

    /** Draws "label (count): [x] |" checkboxes, handing the checked values to onChange */
    function drawFilterCheckboxes(container, entries, isChecked, onChange) {
        container.empty();
        for (const entry of entries) {
            const checkbox = $('<input>', {
                type: 'checkbox',
                class: 'genotype-filter',
                value: entry.value,
            }).prop('checked', isChecked(entry.value));
            container.append($('<label>', {text: `${entry.label} (${entry.count}):`}), checkbox,
                $('<span>', {class: 'separator', text: '|'}));
        }

        container.find(".genotype-filter").change(function () {
            onChange(container.find(".genotype-filter:checked")
                .map(function () { return $(this).val(); }).get());
            dataTable.draw();
        });
    }

    function drawZygosityFilters() {
        const counts = zygosityCounts();
        if (checkedZygosities === null) {
            // Show variants by default, but don't hide everything if that's all there is
            const hasDefaults = DEFAULT_VISIBLE_ZYGOSITIES.some(code => counts[code]);
            const visible = hasDefaults ? DEFAULT_VISIBLE_ZYGOSITIES : ZYGOSITIES.map(z => z.code);
            checkedZygosities = new Set(visible);
        }

        const entries = ZYGOSITIES.filter(z => counts[z.code])
            .map(z => ({value: z.code, label: z.label, count: counts[z.code]}));
        drawFilterCheckboxes($("#zygosity-filters"), entries,
            code => checkedZygosities.has(code),
            values => checkedZygosities = new Set(values));
    }

    function drawGenomeBuildFilters() {
        const counts = genomeBuildCounts();
        const builds = Object.keys(counts).sort();
        for (const build of builds) {
            if (!(build in genomeBuildChecked)) {
                genomeBuildChecked[build] = true;  // Builds load one at a time, so this can be after a filter change
            }
        }

        // A single build is the whole grid, so there'd be nothing to filter
        $("#genome-build-filter").toggle(builds.length > 1);
        const entries = builds.map(build => ({value: build, label: build, count: counts[build]}));
        drawFilterCheckboxes($("#genome-build-filters"), entries,
            build => genomeBuildChecked[build],
            function (values) {
                builds.forEach(build => genomeBuildChecked[build] = values.includes(build));
            });
    }

    function locusCountsTable(locusCounts) {
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
        return table.append(tbody);
    }

    function drawLocusCounts() {
        // Each build has its own locus, so they get their own table rather than being merged
        const container = $("#locus-counts").empty();
        for (const data of responses.values()) {
            if (responses.size > 1) {
                container.append($('<h5>', {text: data.genome_builds.join(', ')}));
            }
            container.append(locusCountsTable(data.locus_counts));
        }
    }

    function drawMultiallelic() {
        const container = $("#other-loci-variants-by-multiallelic").empty();
        const entries = [...responses.values()].flatMap(data => data.multiallelic);
        if (!entries.length) {
            return;
        }
        container.append($('<h4>', {text: 'Variants from multi-allelic VCF records'}));
        container.append($('<p>', {
            text: 'Multi-allelic VCF records are split then separately normalised. This may cause records ' +
                'that were on a line together in a VCF to now have different loci.'
        }));
        for (const entry of entries) {
            container.append($('<b>', {text: entry.multiallelic}));
            const list = $('<ul>');
            for (const variant of entry.variants) {
                list.append($('<li>').append($('<a>', {href: variant.url, text: variant.label})));
            }
            container.append(list);
        }
    }

    function mergedSummary() {
        const genomeBuilds = [];
        const countedBuilds = new Set();  // Builds sharing a contig come back in one response - count them once
        const samples = {total: 0, visible: 0};
        const observations = {total: 0, visible: 0, invisible: 0};
        let partial = false;

        for (const data of responses.values()) {
            const buildKey = data.genome_builds.join(',');
            if (!countedBuilds.has(buildKey)) {
                countedBuilds.add(buildKey);
                data.genome_builds.forEach(b => genomeBuilds.push(b));
                samples.total += data.samples.total;
                samples.visible += data.samples.visible;
            }
            if (data.partial) {
                // Stopped scanning, so its total is the precomputed database-wide one and it never
                // worked out how many of those this user can see
                partial = true;
                observations.total += data.observations.total || 0;
            } else {
                observations.total += data.observations.total;
                observations.visible += data.observations.visible;
                observations.invisible += data.observations.invisible;
            }
        }
        genomeBuilds.sort();
        return {genomeBuilds: genomeBuilds, samples: samples, observations: observations, partial: partial};
    }

    function loadAllRows(button) {
        button.prop('disabled', true).text('Loading...');
        $.when(...truncatedVariantIds().map(variantId => loadVariant(variantId, 0)))
            .fail(() => showError("Error loading all rows"));
    }

    function truncatedRows(observations, partial) {
        const button = $('<button>', {
            type: 'button',
            class: 'btn btn-outline-secondary btn-sm',
            text: 'Load all rows',
        }).click(function () {
            loadAllRows($(this));
        });

        const shown = allRows.length.toLocaleString();
        // A partial scan has nothing but the precomputed counts to give a total, and they may not be there
        const text = observations.total
            ? `Showing ${shown} of ${(partial ? '~' : '') + observations.total.toLocaleString()} rows. `
            : `Showing ${shown} rows, there are more. `;
        const block = $('<p>').append($('<span>', {text: text}), button);
        if (partial) {
            block.append($('<span>', {
                class: 'text-muted',
                text: ' Counts are the precomputed totals for the whole database - loading all rows counts' +
                    ' the ones you have permission to see.',
            }));
        }
        return block;
    }

    function drawSummary() {
        const {genomeBuilds, samples, observations, partial} = mergedSummary();
        const builds = genomeBuilds.join(', ');
        let summary = `Searching ${samples.total} samples for ${builds}`;
        summary += samples.visible < samples.total ? ` (you can see ${samples.visible}).` : '.';

        const container = $("#samples-summary").empty();
        if (!observations.total && !allRows.length) {
            container.append($('<p>', {text: `No results found (searched ${samples.total} samples for ${builds}).`}));
            return;
        }
        container.append($('<p>', {text: summary}));

        if (observations.invisible) {
            const hidden = $('<div>').append(
                $('<img>', {src: config.unknownSampleIconUrl}),
                $('<b>', {text: `You cannot see ${observations.invisible} samples`}),
                $('<p>', {
                    text: `You only have permission to view ${observations.visible} out of ` +
                        `${observations.total} total observations in the database.`
                }));
            container.append(hidden);
        }

        if (truncatedVariantIds().length) {
            container.append(truncatedRows(observations, partial));
        }
    }

    function plotAlleleFrequencies(rows) {
        const x = [];
        const alleleFrequencies = [];
        const labels = [];
        for (const row of rows) {
            if (row.allele_frequency_unit === null || row.allele_frequency_unit === undefined) {
                continue;
            }
            x.push((Math.random() - 0.5) * 0.1);  // jitter
            alleleFrequencies.push(row.allele_frequency_unit);
            let sampleDescription = `Sample: (${row.sample})`;
            if (row.sample_name) {
                sampleDescription += " " + row.sample_name;
            }
            sampleDescription += " DP=" + row.read_depth;
            labels.push(sampleDescription);
        }

        $("#allele-frequency-graphs-container").toggle(alleleFrequencies.length > 0);
        if (!alleleFrequencies.length) {
            return;
        }

        const scatterData = {
            name: "Allele Frequency",
            type: 'scatter',
            mode: 'markers',
            marker: {color: AF_COLOR},
            x: x,
            y: alleleFrequencies,
            text: labels,
        };
        const layout = defaultLayout("Allele Frequency", null, GRAPH_HEIGHT);
        layout.xaxis = Object.assign(layout.xaxis || {}, {zeroline: false, showgrid: false, showticklabels: false});
        layout.yaxis = Object.assign(layout.yaxis || {}, {range: [0, 1.05], showgrid: false, zeroline: false});
        Plotly.newPlot('sample-allele-frequency-scatter', [scatterData], layout, RESPONSIVE);

        // Need to clamp 100 so that histo fits in 10 bins
        const clampedAf = alleleFrequencies.map(af => af >= 1 ? 0.99 : af);
        const histogramData = {
            name: "Allele Frequency",
            type: 'histogram',
            orientation: 'h',
            marker: {color: AF_COLOR},
            y: clampedAf,
            autobiny: false,
            histnorm: "count",
            ybins: {start: 0, size: .05, end: 1},
        };
        const histogramLayout = defaultLayout("Allele Frequency Histogram", null, GRAPH_HEIGHT);
        // Bins run up the same scale as the scatter beside it, so it reads as that plot's summary
        histogramLayout.xaxis = {tickmode: 'auto', showgrid: false, zeroline: false};
        histogramLayout.yaxis = {range: [0, 1.05], showgrid: false, zeroline: false, showticklabels: false};
        histogramLayout.margin = {l: 20};
        Plotly.newPlot('sample-allele-frequency-histogram', [histogramData], histogramLayout, RESPONSIVE);
    }

    function plotPatientTermCounts(rows, selector, title, color, rowField) {
        // The graph counts patients, while the grid is per sample - a patient can have several
        const patientsByTerm = {};
        for (const row of rows) {
            if (!row.patient || !row[rowField]) {
                continue;
            }
            for (const term of row[rowField].split('|')) {
                (patientsByTerm[term] = patientsByTerm[term] || new Set()).add(row.patient);
            }
        }

        const counts = Object.entries(patientsByTerm).map(([term, patients]) => [term, patients.size]);
        // An empty graph would still take up a column of the flex row
        $("#" + selector).toggle(counts.length > 0);
        if (!counts.length) {
            Plotly.purge(selector);
            return false;
        }
        counts.sort((a, b) => b[1] - a[1]);
        const top = counts.slice(0, MAX_GRAPH_TERMS).reverse();  // Plotly draws h-bars bottom up
        plotHBarArrays(selector, title + "<br>Patient Count", top.map(c => c[1]), top.map(c => c[0]),
            null, GRAPH_HEIGHT, color, {l: 200});
        return true;
    }

    function drawPhenotypeGraphs(rows) {
        // The graphs size themselves to their container, so it has to be visible before we plot into it
        const container = $("#variant-sample-phenotype-graphs").show();
        const hpo = plotPatientTermCounts(rows, 'sample-hpo-graph', "Human Phenotype Ontology", HPO_COLOR, 'patient_hpo');
        const omim = plotPatientTermCounts(rows, 'sample-omim-graph', "OMIM", OMIM_COLOR, 'patient_omim');
        const mondo = plotPatientTermCounts(rows, 'sample-mondo-graph', "MONDO", MONDO_COLOR, 'patient_mondo');
        container.toggle(hpo || omim || mondo);
    }

    function drawGraphs() {
        const zygosities = checkedZygosities === null ? '' : [...checkedZygosities].sort().join(',');
        const builds = Object.keys(genomeBuildChecked).filter(b => genomeBuildChecked[b]).sort().join(',');
        const filter = [allRows.length, dataTable.search(), zygosities, builds].join('|');
        if (filter === graphedFilter) {
            return;  // Paging/sorting doesn't change what the graphs summarise
        }
        graphedFilter = filter;

        const rows = visibleRows();
        plotAlleleFrequencies(rows);
        drawPhenotypeGraphs(rows);
    }

    function addResponse(data) {
        // A truncated response is replaced when its rows are loaded in full, so the grid and everything
        // drawn from more than one response is rebuilt rather than appended to
        responses.set(data.variant_id, data);
        config.gHgvs = config.gHgvs || data.g_hgvs;
        allRows = [...responses.values()].flatMap(response => response.rows);

        drawSummary();
        drawLocusCounts();
        drawMultiallelic();

        dataTable.clear();
        if (allRows.length) {
            $("#genotype-grid-container").show();
            drawZygosityFilters();
            drawGenomeBuildFilters();
            dataTable.column('genome_build:name').visible(Object.keys(genomeBuildCounts()).length > 1, false);
            const hasClassifications = allRows.some(row => row.classifications && row.classifications.length);
            dataTable.column('classifications:name').visible(hasClassifications, false);
            dataTable.column('details:name').visible(allRows.some(hasDetails), false);
            dataTable.rows.add(allRows);
            graphedFilter = null;  // The rows changed, so the graphs have to follow
            dataTable.draw();  // 'draw' redraws the graphs
            dataTable.columns.adjust();  // Grid was hidden when created, so scrollX header/body need re-sizing
        }
    }

    function requestComplete() {
        if (--pendingRequests === 0) {
            $("#samples-loading").hide();
            if (!responses.size) {
                $("#samples-summary").html($('<p>', {class: 'text-danger', text: 'Error retrieving samples'}));
            }
        }
    }

    function requestVariants(variantIds) {
        for (const variantId of variantIds || []) {
            if (requestedVariantIds.has(variantId)) {
                continue;
            }
            requestedVariantIds.add(variantId);
            pendingRequests++;
            $("#samples-loading").show();
            loadVariant(variantId).always(requestComplete);
        }
    }

    function init(cfg) {
        config = cfg;
        createDataTable();

        $(document).off(EVENT_NAMESPACE);

        // getOntologyTermLinks makes these as the grid renders, so delegate from the document
        $(document).on('click' + EVENT_NAMESPACE, '.ontology-terms-container .collapsed-term',
                       expandCollapsedOntologyTerm);

        // The allele can be created (linking the other builds' variants) while this section is loading
        $(document).on(ALLELE_VARIANTS_EVENT + EVENT_NAMESPACE, (event, variantIds) => requestVariants(variantIds));

        requestVariants(config.variantIds);
        requestVariants(window.alleleVariantIds);  // Allele resolved before we got here, so we missed the event
    }

    return {
        init: init,
    };
})();
