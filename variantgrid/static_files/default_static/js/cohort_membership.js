/*
 * Cohort membership editor (issue #933)
 *
 * The whole diff (adds, removals, order) is staged client side and committed with one save.
 * CohortGenotype packs every sample into fixed width arrays, so any membership change rebuilds the
 * collection regardless of how many samples moved - batching means one rebuild rather than one per click.
 */

const COHORT_TASK_POLL_FREQUENCY = 1000;

class CohortMembershipEditor {
    constructor(config, sampleRows) {
        this.config = config;
        this.samplesById = {};
        for (const row of sampleRows) {
            this.samplesById[row.id] = row;
        }
        this.candidateIds = [];
        for (const row of (config.candidate_samples || [])) {
            this.samplesById[row.id] = row;
            this.candidateIds.push(row.id);
        }
        this.order = sampleRows.map((row) => row.id);
        this.baseline = this.order.slice();
        this.removed = new Set();
        this.version = config.version;
        this.status = null;

        this.tableContainer = $("#cohort-membership-table");
        this.saveBar = $("#cohort-membership-save-bar");
        this.sampleSelect = $("#id_sample", "#cohort-membership");
        this.vcfSelect = $("#id_vcf", "#cohort-membership");
    }

    init() {
        const that = this;
        this.sampleSelect.change(function () {
            const sampleId = $(this).val();
            if (sampleId) {
                that.stageSamples({sample_ids: JSON.stringify([Number(sampleId)])});
                clearAutocompleteChoice($(this));
            }
        });
        this.vcfSelect.change(function () {
            const vcfId = $(this).val();
            if (vcfId) {
                that.stageSamples({vcf_id: vcfId});
                clearAutocompleteChoice($(this));
            }
        });

        if (this.config.import_status === "I") {
            this.status = "building";
            if (this.config.running_celery_task) {
                this.pollTask(this.config.running_celery_task);
            }
        } else if (this.config.import_status === "E") {
            this.status = "failed";
        } else if (this.config.has_genotype_data) {
            this.status = "ready";
        } else if (this.order.length) {
            this.status = "unbuilt";
        }
        this.render();
    }

    /* Members after pending removals are applied - what would be saved */
    effectiveMembers() {
        return this.order.filter((id) => !this.removed.has(id));
    }

    diff() {
        const members = this.effectiveMembers();
        const memberSet = new Set(members);
        const baselineSet = new Set(this.baseline);
        const added = members.filter((id) => !baselineSet.has(id));
        const removed = this.baseline.filter((id) => !memberSet.has(id));
        const orderChanged = !added.length && !removed.length &&
            JSON.stringify(members) !== JSON.stringify(this.baseline);
        return {
            members: members,
            added: added,
            removed: removed,
            orderChanged: orderChanged,
            dirty: Boolean(added.length || removed.length || orderChanged),
        };
    }

    /* Whether saving this membership takes the free sub cohort path or a genotype rebuild */
    saveMode(members) {
        if (!members.length) {
            return {key: "empty", html: "<b>Cohort has no samples</b> - add at least one before saving."};
        }
        const vcfIds = new Set(members.map((id) => this.samplesById[id].vcf_id));
        if (vcfIds.size === 1) {
            const vcfName = this.samplesById[members[0]].vcf_name;
            return {
                key: "sub",
                html: `All samples come from <b>${escapeHtml(vcfName)}</b> - saved as a <b>sub cohort</b>, ready instantly.`,
            };
        }
        return {
            key: "rebuild",
            html: `Samples come from <b>${vcfIds.size} VCFs</b> - saving rebuilds the genotype data in the background.`,
        };
    }

    /* Duo/Trio/Quad cascade off CohortSample, so a removal can take one with it - say so before saving */
    familyGroupsLostBy(removedSampleIds) {
        const familyGroups = this.config.family_groups || {};
        const labels = [];
        for (const sampleId of removedSampleIds) {
            for (const label of (familyGroups[sampleId] || [])) {
                if (!labels.includes(label)) {
                    labels.push(label);
                }
            }
        }
        return labels;
    }

    stageSamples(params) {
        const that = this;
        $.getJSON(this.config.sample_rows_url, params, function (data) {
            for (const row of data.samples) {
                that.samplesById[row.id] = row;
                if (that.removed.has(row.id)) {
                    that.removed.delete(row.id);
                } else if (!that.order.includes(row.id)) {
                    that.order.push(row.id);
                }
            }
            that.status = null;
            that.render();
        });
    }

    renderTable() {
        const baselineSet = new Set(this.baseline);
        const table = $(`<table class="table membership-table">
            <thead><tr><th></th><th>Sample</th><th>Source VCF</th><th>Het / Hom</th><th>Sex</th><th></th></tr></thead>
            <tbody></tbody></table>`);
        const tbody = $("tbody", table);

        for (const sampleId of this.order) {
            const sample = this.samplesById[sampleId];
            const isNew = !baselineSet.has(sampleId);
            const isRemoved = this.removed.has(sampleId);
            const rowClass = isRemoved ? "pending-remove" : (isNew ? "pending-add" : "");
            const newChip = isNew ? ' <span class="badge badge-success chip-new">new</span>' : "";
            const button = isRemoved
                ? `<button type="button" class="btn btn-sm btn-outline-success undo-sample" data-sample-id="${sampleId}">Undo</button>`
                : `<button type="button" class="btn btn-sm btn-outline-danger remove-sample" data-sample-id="${sampleId}">Remove</button>`;
            tbody.append(`<tr draggable="true" data-sample-id="${sampleId}" class="${rowClass}">
                <td class="drag-handle" title="Drag to reorder">&#8942;&#8942;</td>
                <td class="sample-name"><a href="${sample.url}">${escapeHtml(sample.name)}</a>${newChip}</td>
                <td class="vcf-name">${escapeHtml(sample.vcf_name)}</td>
                <td class="num">${sample.het_hom || "N/A"}</td>
                <td>${escapeHtml(sample.sex)}</td>
                <td class="text-right">${this.config.has_write_permission ? button : ""}</td>
            </tr>`);
        }
        if (!this.order.length) {
            tbody.append('<tr><td colspan="6" class="text-muted">No samples - add some below.</td></tr>');
        }

        // Sub cohort: the rest of the parent VCF's samples, ready to tick back in
        const memberSet = new Set(this.order);
        const available = this.candidateIds.filter((id) => !memberSet.has(id));
        if (available.length) {
            tbody.append(`<tr class="candidate-heading"><td colspan="6">Other samples in this VCF</td></tr>`);
            for (const sampleId of available) {
                const sample = this.samplesById[sampleId];
                const button = `<button type="button" class="btn btn-sm btn-outline-primary add-sample" data-sample-id="${sampleId}">Add</button>`;
                tbody.append(`<tr class="candidate-row" data-candidate-id="${sampleId}">
                    <td></td>
                    <td class="sample-name"><a href="${sample.url}">${escapeHtml(sample.name)}</a></td>
                    <td class="vcf-name">${escapeHtml(sample.vcf_name)}</td>
                    <td class="num">${sample.het_hom || "N/A"}</td>
                    <td>${escapeHtml(sample.sex)}</td>
                    <td class="text-right">${this.config.has_write_permission ? button : ""}</td>
                </tr>`);
            }
        }

        this.tableContainer.empty().append(table);
        if (this.config.has_write_permission) {
            this.bindRowHandlers(tbody);
        }
    }

    bindRowHandlers(tbody) {
        const that = this;
        $(".remove-sample", tbody).click(function () {
            that.removed.add(Number($(this).data("sample-id")));
            that.status = null;
            that.render();
        });
        $(".undo-sample", tbody).click(function () {
            that.removed.delete(Number($(this).data("sample-id")));
            that.status = null;
            that.render();
        });
        $(".add-sample", tbody).click(function () {
            that.order.push(Number($(this).data("sample-id")));
            that.status = null;
            that.render();
        });

        let dragRow = null;
        $("tr[draggable]", tbody).each(function () {
            const row = this;
            row.addEventListener("dragstart", (e) => {
                dragRow = row;
                $(row).addClass("dragging");
                e.dataTransfer.effectAllowed = "move";
                e.dataTransfer.setData("text/plain", $(row).data("sample-id"));
            });
            row.addEventListener("dragover", (e) => {
                e.preventDefault();
                if (!dragRow || dragRow === row) {
                    return;
                }
                const rect = row.getBoundingClientRect();
                const after = e.clientY > rect.top + rect.height / 2;
                tbody[0].insertBefore(dragRow, after ? row.nextSibling : row);
            });
            row.addEventListener("dragend", () => {
                $(dragRow).removeClass("dragging");
                dragRow = null;
                that.order = $("tr[data-sample-id]", tbody).map(function () {
                    return Number($(this).data("sample-id"));
                }).get();
                that.status = null;
                that.render();
            });
        });
    }

    renderSaveBar() {
        const diff = this.diff();
        const mode = this.saveMode(diff.members);

        let pill;
        if (diff.dirty) {
            const parts = [];
            if (diff.added.length) {
                parts.push(`<span class="diff-added">+${diff.added.length} added</span>`);
            }
            if (diff.removed.length) {
                parts.push(`<span class="diff-removed">&minus;${diff.removed.length} removed</span>`);
            }
            if (diff.orderChanged) {
                parts.push('<span class="diff-order">order changed</span>');
            }
            pill = parts.join(" &middot; ");
        } else {
            pill = '<span class="text-muted">No pending changes</span>';
        }

        const versionText = diff.dirty
            ? `v${this.version} &rarr; v${this.version + 1} on save (one bump)`
            : `v${this.version}`;
        const canSave = this.config.has_write_permission && diff.dirty &&
            mode.key !== "empty" && this.status !== "building";
        const canDiscard = diff.dirty && this.status !== "building";

        const lostFamilyGroups = this.familyGroupsLostBy(diff.removed);
        const familyGroupWarning = lostFamilyGroups.length
            ? `<div class="family-group-warning"><i class="fas fa-exclamation-triangle"></i>
               Saving also deletes ${escapeHtml(lostFamilyGroups.join(", "))} - built on the samples you removed.</div>`
            : "";

        const bar = $(`<div class="save-bar">
            <div class="save-bar-top">
                <span class="diff-pill">${pill}</span>
                <span class="version-tag text-muted">${versionText}</span>
                <span class="save-bar-actions">
                    <button type="button" class="btn btn-outline-secondary" id="discard-membership-changes">Discard</button>
                    <button type="button" class="btn btn-primary" id="save-cohort-membership">Save cohort</button>
                </span>
            </div>
            <div class="mode-line mode-${mode.key}">${mode.html}</div>
            ${familyGroupWarning}
            <div class="save-bar-status"></div>
        </div>`);

        $("#save-cohort-membership", bar).prop("disabled", !canSave);
        $("#discard-membership-changes", bar).prop("disabled", !canDiscard);
        $(".save-bar-status", bar).html(this.statusHtml());

        this.saveBar.empty().append(bar);

        const that = this;
        $("#save-cohort-membership", bar).click(() => that.save());
        $("#discard-membership-changes", bar).click(() => that.discard());
        $("#retry-cohort-genotype", bar).click(() => that.launchGenotypeTask(that.config.create_genotype_url, {}));
    }

    statusHtml() {
        switch (this.status) {
            case "building":
                return '<i class="fas fa-spinner fa-spin"></i> Building genotype data - this cohort can be used in analyses when it finishes.';
            case "saved":
                return '<i class="fas fa-check text-success"></i> Saved - cohort is ready for analyses.';
            case "ready":
                return '<i class="fas fa-check text-success"></i> Genotype data is built - cohort is ready for analyses.';
            case "unbuilt":
                return "Genotype data has not been built for this membership yet - press <b>Save cohort</b> to build it.";
            case "failed":
                return `<span class="text-danger">Building genotype data failed.</span>
                        <button type="button" class="btn btn-sm btn-outline-secondary" id="retry-cohort-genotype">Try again</button>`;
            default:
                return "";
        }
    }

    render() {
        this.renderTable();
        this.renderSaveBar();
    }

    discard() {
        this.order = this.baseline.slice();
        this.removed.clear();
        this.status = this.config.has_genotype_data ? "ready" : "unbuilt";
        this.render();
    }

    save() {
        const members = this.effectiveMembers();
        this.launchGenotypeTask(this.config.save_url, {sample_ids: JSON.stringify(members)}, members);
    }

    launchGenotypeTask(url, data, savedMembers) {
        const that = this;
        this.status = "building";
        this.render();
        $.ajax({
            type: "POST",
            url: url,
            data: data,
            success: function (response) {
                if (savedMembers) {
                    that.baseline = savedMembers.slice();
                    that.order = savedMembers.slice();
                    that.removed.clear();
                    that.version += 1;
                }
                if (response.status) {
                    that.taskComplete(response.status === "SUCCESS");
                } else if (response.celery_task) {
                    that.pollTask(response.celery_task);
                } else {
                    that.taskComplete(true);
                }
            },
            error: function () {
                that.status = "failed";
                that.render();
            },
        });
    }

    pollTask(celeryTask) {
        const that = this;
        this.status = "building";
        this.render();
        $.getJSON(Urls.job_state(celeryTask), function (data) {
            if (data.status === "SUCCESS") {
                that.taskComplete(true);
            } else if (data.status === "FAILURE") {
                that.taskComplete(false);
            } else {
                window.setTimeout(() => that.pollTask(celeryTask), COHORT_TASK_POLL_FREQUENCY);
            }
        });
    }

    taskComplete(success) {
        if (success) {
            this.status = "saved";
            this.render();
            // Details tab, related data and analyses all read the membership we just changed
            window.location.reload();
        } else {
            this.status = "failed";
            this.render();
        }
    }
}
