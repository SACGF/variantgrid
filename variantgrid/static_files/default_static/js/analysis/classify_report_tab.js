// @ts-check
// analysis/templates/analysis/classify_report_tab.html
$(document).ready(function() {
    const POPULATING_POLL_MS = 3000;
    const POPULATING_MAX_POLLS = 20;
    const container = $("#classify-report");
    const tabPane = container.parent();  // Reloaded in place, the way the other sample tabs do it
    // The wizard walks the rows still wanting a classification, one dialog at a time
    let wizard = null;
    let reloadWhenDone = false;

    // The modals live inside the tab, so replacing it under an open one strands Bootstrap's backdrop -
    // a dialog opened while the reload was in flight holds it off until that dialog closes
    function reloadTab() {
        $.get(container.data("tab-url"), function(html) {
            if ($(".modal.show", tabPane).length) {
                reloadWhenDone = true;
            } else {
                tabPane.html(html);
            }
        });
        $(document).trigger("classifyReportChanged");  // The tab label's counts follow the tab
    }

    // The walk is about making classifications - a row that already has one is waiting on "Clear tag"
    function unclassifiedTagIds() {
        return $(".classify-queue-row:not(.has-classification)", container).map(function() {
            return $(this).data("variantTagId");
        }).get();
    }

    // Classified rows drop out of the queue as they're done. Kept on the tab pane so it survives reloadTab()
    // but not a page load
    function updateHideClassified() {
        const hide = tabPane.data("showClassified") !== true;
        const rows = $(".classify-queue-row", container);
        $("#hide-classified", container).prop("checked", hide);
        container.toggleClass("hide-classified", hide);
        container.toggleClass("all-classified", rows.length > 0 && rows.not(".classified").length === 0);
    }

    container.on("change", "#hide-classified", function() {
        tabPane.data("showClassified", !$(this).is(":checked"));
        updateHideClassified();
    });
    updateHideClassified();

    function openClassifyDialog(variantTagId) {
        const row = $(`.classify-queue-row[data-variant-tag-id='${variantTagId}']`, container);
        const url = row.data("dialogUrl");
        const modal = $("#classify-modal");
        $(".modal-body", modal).load(url, function() {
            // Skip is only meaningful while stepping through the queue
            $(".wizard-skip", modal).toggle(Boolean(wizard));
            if (wizard) {
                const position = wizard.indexOf(variantTagId) + 1;
                $(".wizard-position", modal).text(`tag ${position} of ${wizard.length}`);
            }
        });
        modal.data("variantTagId", variantTagId).modal("show");
    }

    // Always forwards through the queue, so a skipped tag isn't offered again on the next step
    function nextInWizard(afterTagId) {
        if (!wizard) return null;
        const unclassified = new Set(unclassifiedTagIds());
        for (let i = wizard.indexOf(afterTagId) + 1; i < wizard.length; i++) {
            if (unclassified.has(wizard[i])) {
                return wizard[i];
            }
        }
        return null;
    }

    function markRowClassified(variantTagId, data) {
        const row = $(`.classify-queue-row[data-variant-tag-id='${variantTagId}']`, container);
        const action = $(".classify-action", row);
        action.html($("<span>", {class: "text-success mr-2", text: "✓ Classified"})).append(
            $("<a>", {href: data.url, target: "_blank", class: "classification-link", text: data.label}));
        row.addClass("has-classification");
        // Made for a sample other than the tagging's own, it stays outstanding until someone says this is the person
        if (data.resolved) {
            row.addClass("classified table-success");
            updateHideClassified();
        } else {
            action.append($("<button>", {
                type: "button", class: "btn btn-outline-success btn-sm ml-2 resolve-tag",
                title: "Say this classification is what the tag was asking for", text: "Clear tag"
            }));
        }
        reloadWhenDone = true;
    }

    // Whether stepping on or finishing, the modal closes on the last tag of the walk
    function stepWizard(afterTagId) {
        const next = nextInWizard(afterTagId);
        if (next) {
            openClassifyDialog(next);
        } else {
            $("#classify-modal").modal("hide");
        }
    }

    function afterClassified(variantTagId, data) {
        markRowClassified(variantTagId, data);
        stepWizard(variantTagId);
    }

    // Copy a previously curated record into a new one for this sample - the whole record when the
    // allele has been curated before, otherwise just the gene content
    function applyPrevious(button, fieldName) {
        const form = $("#classify-form");
        const variantTagId = $("#classify-modal").data("variantTagId");
        const sampleSelect = form.find("[name=sample_id]");
        if (!sampleSelect.val()) {
            sampleSelect.addClass("is-invalid").focus();
            return;
        }
        form.find(`[name=${fieldName}]`).val(button.data("vcmId"));
        $("button", form).prop("disabled", true);
        button.text("Creating\u2026");
        $.post(form.attr("action"), form.serialize(), function(data) {
            afterClassified(variantTagId, data);
        }).fail(function() {
            $("button", form).prop("disabled", false);
            window.alert("Could not create the classification - please use the full form.");
        });
    }

    container.on("click", ".apply-previous", function() {
        applyPrevious($(this), "copy_from_vcm_id");
    });

    container.on("click", ".apply-gene", function() {
        applyPrevious($(this), "copy_gene_from_vcm_id");
    });

    container.on("click", ".resolve-tag", function() {
        const button = $(this);
        const row = button.closest(".classify-queue-row");
        button.prop("disabled", true);
        $.post($(".classify-action", row).data("resolveUrl"), function() {
            button.remove();
            row.addClass("classified table-success");
            updateHideClassified();
            reloadWhenDone = true;
        }).fail(function() {
            button.prop("disabled", false);
            window.alert("Could not clear the tag.");
        });
    });

    // The full form opens in its own tab, so the row is picked up when they come back to this one.
    // The flag lives on the tab pane, which outlives reloadTab() replacing the tab's contents
    container.on("click", ".classify-full-form", function() {
        tabPane.data("reloadOnReturn", true);
        if (wizard) {
            // Sending a tag off to the full form doesn't classify it here - the walk moves on, as Skip does
            stepWizard($("#classify-modal").data("variantTagId"));
        } else {
            $("#classify-modal").modal("hide");
        }
    });

    $(window).off("focus.classifyReport").on("focus.classifyReport", function() {
        if (tabPane.data("reloadOnReturn") && !$("#classify-modal").hasClass("show")) {
            tabPane.data("reloadOnReturn", false);
            reloadTab();
        }
    });

    // Editing, submitting or resolving a record's errors happens on its own page, so the rows are picked up on the way back
    container.on("click", ".report-fix-link, .classification-link", function() {
        tabPane.data("reloadOnReturn", true);
    });

    container.on("click", ".classify-tag", function() {
        wizard = null;
        openClassifyDialog($(this).closest(".classify-queue-row").data("variantTagId"));
    });

    container.on("click", "#classify-all", function() {
        const unclassified = unclassifiedTagIds();
        if (unclassified.length) {
            wizard = unclassified;
            openClassifyDialog(unclassified[0]);
        }
    });

    container.on("click", ".wizard-skip", function() {
        stepWizard($("#classify-modal").data("variantTagId"));
    });

    // A record created from the queue only reaches Classifications once celery has published it. A timer
    // outliving its tab (replaced by some other reload) lets the newer tab do the polling
    if (container.data("populating")) {
        const polls = (tabPane.data("populatingPolls") || 0) + 1;
        if (polls <= POPULATING_MAX_POLLS) {
            tabPane.data("populatingPolls", polls);
            setTimeout(function() {
                if ($.contains(document, container[0])) {
                    reloadTab();
                }
            }, POPULATING_POLL_MS);
        }
    } else {
        tabPane.data("populatingPolls", 0);
    }

    $("#classify-modal").on("hidden.bs.modal", function() {
        wizard = null;
        if (reloadWhenDone) {
            reloadTab();
        }
    });

    function updateReportSelection() {
        const selected = $(".report-select:checked", container).length;
        $("#report-selected-count", container).text(`${selected} selected`);
        $("#build-case-report", container).prop("disabled", selected === 0)
            .text(`Build multiple variant report (${selected})`);
    }

    container.on("change", ".report-select", updateReportSelection);
    container.on("change", "#report-select-all", function() {
        $(".report-select:not(:disabled)", container).prop("checked", $(this).is(":checked"));
        updateReportSelection();
    });
    updateReportSelection();

    // The case report modal is server rendered throughout: the ticked classifications are
    // posted to the build dialog, changing the template re-posts it (a different template
    // asks for different case_fields), and building replaces the body with the preview
    const caseReportModal = $("#case-report-modal");

    // load() only POSTs an object, and the dialog is POST only - so serializeArray, not serialize
    function loadCaseReportDialog(url, form) {
        $(".modal-body", caseReportModal).load(url, form.serializeArray());
    }

    container.on("click", "#build-case-report", function() {
        loadCaseReportDialog($(this).data("dialogUrl"), $("#report-form", container));
        caseReportModal.modal("show");
    });

    caseReportModal.on("change", ".case-report-template", function() {
        const form = $("#case-report-form", caseReportModal);
        loadCaseReportDialog(form.data("dialogUrl"), form);
    });

    caseReportModal.on("click", ".build-case-report", function() {
        const button = $(this);
        const form = $("#case-report-form", caseReportModal);
        button.prop("disabled", true).text("Building\u2026");
        $.post(form.data("buildUrl"), form.serialize(), function(html) {
            $(".modal-body", caseReportModal).html(html);
            reloadWhenDone = true;
        }).fail(function() {
            button.prop("disabled", false).text("Build multiple variant report");
            window.alert("Could not build the report.");
        });
    });

    // Finalise / Rebuild documents / New version all act on one report and then reload the
    // tab, so the Reports card is redrawn from the database rather than patched here
    container.on("click", ".case-report-action", function() {
        const button = $(this);
        const confirmMessage = button.data("confirm");
        if (confirmMessage && !window.confirm(confirmMessage)) {
            return;
        }
        button.prop("disabled", true);
        $.post(button.data("actionUrl"), {csrfmiddlewaretoken: $("[name=csrfmiddlewaretoken]", container).first().val()},
               function(data) {
            if (data.unsubmitted && data.unsubmitted.length) {
                window.alert("Left alone - these records have unsubmitted changes, so publishing them " +
                             "would push out edits nobody has submitted:\n" + data.unsubmitted.join("\n"));
            }
            reloadTab();
        }).fail(function() {
            button.prop("disabled", false);
            window.alert("Could not action this report.");
        });
    });

    container.on("submit", ".case-report-lis-form", function(event) {
        event.preventDefault();
        const form = $(this);
        $.post(form.attr("action"), form.serialize(), function() {
            reloadTab();
        }).fail(function() {
            window.alert("Could not save the report ID.");
        });
    });

    $("#case-report-modal").on("hidden.bs.modal", function() {
        if (reloadWhenDone) {
            reloadTab();
        }
    });
});
