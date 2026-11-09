// @ts-check
/* global $, Cookies */

/**
 * The classification form's Change sample dialog - @see views_classification_sample.py
 * Opened by the #sample-change-open button vc_form.js puts in the summary card's Sample row, only for a user
 * who may reassign (the only page this script is on).
 * @param {{sampleUrl: string, sampleSelectId: string}} options
 */
function initClassificationSample(options) {
    const modal = $('#sample-change-modal');
    const sampleSelect = $(document.getElementById(options.sampleSelectId));
    const preview = $('#sample-change-preview', modal);
    const repopulate = $('#sample-change-repopulate', modal);
    const overrideRow = $('#sample-change-override-row', modal);
    const override = $('#sample-change-override', modal);
    const confirmButton = $('#sample-change-confirm', modal);

    function updateConfirm() {
        const ready = $('#sample-change-preview-ready', preview);
        const needsOverride = ready.length > 0 && ready.data('hasGenotype') === false;
        overrideRow.toggleClass('d-none', !needsOverride);
        confirmButton.prop('disabled', ready.length === 0 || (needsOverride && !override.prop('checked')));
    }

    function loadPreview() {
        const sampleId = sampleSelect.val();
        preview.empty();
        override.prop('checked', false);
        updateConfirm();
        if (!sampleId) {
            return;
        }
        preview.html('<span class="text-muted">Loading…</span>');
        preview.load(options.sampleUrl + '?' + $.param({sample_id: sampleId}), function(response, status) {
            if (status === 'error') {
                preview.html('<div class="alert alert-danger">Could not load the changes for this sample.</div>');
            }
            updateConfirm();
        });
    }

    $('#sample-change-form', modal).on('submit', function(event) {
        event.preventDefault();
    });

    // Searches start within the current sample's patient (the forwarded hidden input) - clearing it widens
    // the next search to every sample the user can see
    $('#sample-change-show-all', modal).on('click', function(event) {
        event.preventDefault();
        $('#sample-change-patient', modal).val('');
        $('#sample-change-patient-scope', modal).remove();
        sampleSelect.val(null).trigger('change');
        sampleSelect.select2('close');
        sampleSelect.select2('open');
    });

    // vc_form.js re-renders the summary card, so the button is bound by delegation
    $('#vcform-summary').on('click', '#sample-change-open', function() {
        modal.modal('show');
    });
    sampleSelect.on('change', loadPreview);
    override.on('change', updateConfirm);
    repopulate.on('change', function() {
        preview.toggleClass('d-none', !repopulate.prop('checked'));
    });

    confirmButton.on('click', function() {
        confirmButton.prop('disabled', true).text('Changing…');
        /** @type {Object<string, any>} */
        const data = {
            sample_id: sampleSelect.val(),
            csrfmiddlewaretoken: Cookies.get('csrftoken'),
        };
        if (repopulate.prop('checked')) {
            data.repopulate = 'on';
        }
        if (override.prop('checked')) {
            data.confirm_no_genotype = 'on';
        }
        $.post(options.sampleUrl, data, function() {
            window.location.reload();
        }).fail(function() {
            confirmButton.text('Confirm');
            updateConfirm();
            window.alert('Could not change the sample.');
        });
    });
}
