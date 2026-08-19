/**
 * Tag stats page - @see https://github.com/SACGF/variantgrid/issues/1751
 *
 * Every card lazy loads its own JSON endpoint so the page paints straight away. Responses are cached
 * server side, so each card shows when its numbers were calculated.
 */

const TAG_STATS_COLORS = ['#4c78a8', '#f58518', '#54a24b', '#e45756', '#72b7b2',
                          '#ff9da6', '#9d755d', '#b279a2', '#eeca3b', '#bab0ac'];

function tagStatsTimeAgo(isoString) {
    const seconds = Math.max(0, (new Date() - new Date(isoString)) / 1000);
    const units = [["day", 86400], ["hour", 3600], ["minute", 60]];
    for (const [name, unitSeconds] of units) {
        const amount = Math.floor(seconds / unitSeconds);
        if (amount >= 1) {
            return amount + " " + name + (amount === 1 ? "" : "s") + " ago";
        }
    }
    return "just now";
}

function tagStatsLayout(title, extra) {
    return Object.assign({
        title: title,
        height: 400,
        margin: {t: 50, b: 100},
        legend: {orientation: "h"},
        colorway: TAG_STATS_COLORS,
    }, extra || {});
}

/** Fetch a card's JSON then hand it to render(data, $content) */
function loadTagStatsCard(cardId, url, render) {
    const $card = $("#" + cardId);
    const $content = $card.find(".card-content");
    const $spinner = $card.find(".card-spinner");

    $spinner.show();
    $content.hide();
    $.getJSON(url).done((data) => {
        $spinner.hide();
        $card.find(".card-calculated").text(tagStatsTimeAgo(data.calculated));
        render(data, $content);
        $content.show();
    }).fail(() => {
        $spinner.hide();
        $content.html('<div class="alert alert-danger">Could not load this card</div>').show();
    });
}

/** Tags that failed liftover have no allele, so they sit outside every allele based count */
function noAlleleWarning(count) {
    if (!count) {
        return "";
    }
    return `<p class="text-muted small mt-2">${count.toLocaleString()} tag
            ${count === 1 ? "event is" : "events are"} not included - liftover failed, so they have no allele.</p>`;
}

function renderTagStatsHeadline(data, $content) {
    const stats = [
        ["Tag events", data.tag_events, "Every time someone applied a tag"],
        ["Distinct (allele, tag)", data.distinct_allele_tags, "Ignoring who tagged it, and when"],
        ["Tagged alleles", data.distinct_alleles, "Alleles carrying at least one tag"],
        ["Analyses", data.analyses, "Analyses tags were applied in"],
        ["Taggers this year", data.taggers_this_year, "Users who tagged in the last 365 days"],
    ];
    const cells = stats.map(([label, value, help]) =>
        `<div class="col text-center" title="${help}">
            <div class="h3 mb-0">${value.toLocaleString()}</div>
            <div class="text-muted small">${label}</div>
         </div>`);
    $content.html('<div class="row">' + cells.join("") + '</div>' + noAlleleWarning(data.no_allele));
}

function renderTagStatsOverTime(data, $content) {
    $content.html('<div id="tags-over-time-graph"></div><div id="tags-proportion-graph"></div>');

    const traces = data.series.map((s) => ({x: data.months, y: s.counts, name: s.name, type: 'bar'}));
    Plotly.newPlot('tags-over-time-graph', traces,
                   tagStatsLayout("Tag events per month", {barmode: 'stack'}), {responsive: true});

    // Deployments have a long tail of one-off tags - collapse anything under 1% into "other"
    const total = Object.values(data.totals).reduce((a, b) => a + b, 0);
    const labels = [];
    const values = [];
    let other = 0;
    for (const [label, count] of Object.entries(data.totals)) {
        if (count < total * 0.01) {
            other += count;
        } else {
            labels.push(label);
            values.push(count);
        }
    }
    if (other) {
        labels.push("other");
        values.push(other);
    }
    Plotly.newPlot('tags-proportion-graph', [{labels: labels, values: values, type: 'pie'}],
                   tagStatsLayout("Share of all tag events"), {responsive: true});
}

function renderTagStatsUser(data, $content) {
    const topTags = data.top_tags.map((t) => `<li>${t.tag} &times; ${t.count.toLocaleString()}</li>`).join("");
    $content.html(
        `<p><strong>${data.tag_events.toLocaleString()}</strong> tag events on
         <strong>${data.distinct_alleles.toLocaleString()}</strong> alleles</p>
         <ul>${topTags || "<li>No tags yet</li>"}</ul>
         ${noAlleleWarning(data.no_allele)}
         <div id="user-tags-sparkline"></div>`);

    Plotly.newPlot('user-tags-sparkline',
                   [{x: data.recent_days, y: data.recent_counts, type: 'bar'}],
                   tagStatsLayout("Tag events in the last 30 days", {height: 250}), {responsive: true});
}

function renderTagStatsByLab(data, $content) {
    const userTable = (users) => {
        if (!users.length) {
            return '<p class="text-muted">No tag events.</p>';
        }
        const rows = users.map((u) =>
            `<tr><td>${u.user}</td><td>${u.labs.join(", ")}</td><td>${u.count.toLocaleString()}</td></tr>`).join("");
        return `<table class="table table-sm"><thead><tr><th>User</th><th>Labs</th><th>Tag events</th></tr></thead>
                <tbody>${rows}</tbody></table>`;
    };
    const labRows = data.labs.map((l) =>
        `<tr><td>${l.lab}</td><td>${l.count.toLocaleString()}</td></tr>`).join("");

    $content.html(
        `<h5>Top users - all time</h5>` + userTable(data.top_users) +
        `<h5>Top users - last 30 days</h5>` + userTable(data.top_users_recent) +
        `<h5>By lab</h5>
         <p class="text-muted small">A user in more than one lab counts into each, so lab totals can add up to
         more than the number of tag events.</p>
         <table class="table table-sm"><thead><tr><th>Lab</th><th>Tag events</th></tr></thead>
         <tbody>${labRows}</tbody></table>`);
}

function renderTagStatsGenes(data, $content) {
    $content.html('<div id="tag-genes-graph"></div>');
    if (!data.genes.length) {
        $content.html('<p>No tagged variants matched.</p>');
        return;
    }
    const traces = data.series.map((s) => ({y: data.genes, x: s.counts, name: s.name,
                                            type: 'bar', orientation: 'h'}));
    Plotly.newPlot('tag-genes-graph', traces,
                   tagStatsLayout("Tag events per gene", {barmode: 'stack',
                                                          height: 100 + 40 * data.genes.length,
                                                          yaxis: {automargin: true}}),
                   {responsive: true});
}

function renderTagStatsCoOccurrence(data, $content, tagsUrlFunc) {
    const pairs = data.top_pairs.map((p) => {
        const url = tagsUrlFunc(p.tags);
        return `<li><a href="${url}">${p.tags.join(" + ")}</a>: ${p.alleles.toLocaleString()} alleles</li>`;
    }).join("");

    $content.html('<div id="tag-co-occurrence-graph"></div>' +
                  '<h5>Most common pairs</h5><ul class="columns-2">' + pairs + '</ul>');

    Plotly.newPlot('tag-co-occurrence-graph',
                   [{z: data.matrix, x: data.tags, y: data.tags, type: 'heatmap', colorscale: 'YlOrRd'}],
                   tagStatsLayout("Alleles carrying both tags",
                                  {height: 600, xaxis: {automargin: true}, yaxis: {automargin: true}}),
                   {responsive: true});
}

function renderTagStatsReTagged(data, $content, alleleUrlFunc) {
    if (!data.alleles.length) {
        $content.html(`<p>Nothing tagged '${data.tag}'.</p>`);
        return;
    }
    const rows = data.alleles.map((a) =>
        `<tr><td><a href="${alleleUrlFunc(a.allele_id)}">${a.variant}</a></td>
         <td>${a.count.toLocaleString()}</td></tr>`).join("");
    const percent = Math.round(100 * data.top_events / data.tag_events);

    $content.html(
        `<p>These ${data.alleles.length} alleles account for
         <strong>${data.top_events.toLocaleString()}</strong> of the
         ${data.tag_events.toLocaleString()} '${data.tag}' tag events (${percent}%), spread over
         ${data.distinct_alleles.toLocaleString()} distinct alleles.</p>
         <p class="text-muted">A Tags node set to "Exclude" filters these out of an analysis, so known
         artefacts don't need re-confirming case by case.</p>
         <table class="table table-sm"><thead><tr><th>Variant</th><th>Tag events</th></tr></thead>
         <tbody>${rows}</tbody></table>`);
}

function renderTagStatsTagGenesOverTime(data, $content) {
    $content.html('<div id="tag-genes-over-time-graph"></div>');
    const traces = data.series.map((s) => ({x: data.quarters, y: s.counts, name: s.name,
                                            type: 'scatter', stackgroup: 'one'}));
    Plotly.newPlot('tag-genes-over-time-graph', traces,
                   tagStatsLayout(`'${data.tag}' tag events per gene, by quarter`), {responsive: true});
}
