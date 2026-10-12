// @ts-check
// classification/templates/classification/classification_import_tool.html
let waitingForKeys = true;

function generateRequest() {
    if (waitingForKeys) {
        return;
    }
    const test = $('#test').prop('checked');
    let lab = $('#lab').val();
    const operation = $('#operation').val();
    let record_id = $('#record_id').val();
    const delete_mode = $('#delete').prop('checked');
    const delete_reason = $('#delete_reason').val();
    const publish = $('#publish').val();
    let data = $('#data').val().trim();
    if (data.length === 0) {
        data = '{}';
    }
    const genome_build = $('#genome_build option:selected').val();
    const c_hgvs = $('#c_hgvs').val();
    let valid = lab && record_id;
    const data_fields = $('.data-field');
    const data_field_values = {};
    for (let data_field of data_fields.toArray()) {
        data_field = $(data_field);
        const value = data_field.val();
        if (value.length) {
            data_field_values[data_field.attr('id')] = value;
        }
    }
    const custom_field_values = {};
    const custom_fields = $('.custom-value');
    for (let custom_field of custom_fields.toArray()) {
        custom_field = $(custom_field);
        const dataIndex = custom_field.attr('data-index');
        const custom_key = $(`#custom-key-${dataIndex}`);
        const key = custom_key.val();
        let value = custom_field.val();
        if (value.length == 0) {
            value = null;
        }
        if (key.length) {
            custom_field_values[key] = value;
        }
    }

    record_id = record_id || 'no_record_id';
    lab = lab || 'no_org/no_lab';

    let envelope = {
        "id": `${lab}/${record_id}`
    };
    if (test) {
        envelope["test"] = true;   
    }
    if (delete_mode) {
        envelope["delete"] = true;
        envelope["delete_reason"] = delete_reason;
    } else {
        try {
            const special_values = {
                "c_hgvs": c_hgvs,
                "genome_build": genome_build
            };
            data = Object.assign({}, special_values, data_field_values, custom_field_values, JSON.parse(data));
        } catch (e) {
            data = {
                'invalid_json': e.message
            };
            valid = false;
        }
        envelope[operation] = data;

        if (publish) {
            envelope["publish"] = publish;
        }
    }
    const prettyHtml = formatJson(envelope);
    const previewDom =  $('#request_preview');
    previewDom.html(prettyHtml);

    if (!valid) {
        previewDom.addClass('invalid');
    } else {
        previewDom.removeClass('invalid');
    }

    envelope = {
        "records": [envelope]
    };
    if ($('#import_recording').prop('checked')) {
        envelope["import_id"] = "classification_import_tool";
        envelope["status"] = "complete";
    }

    return [JSON.stringify(envelope), valid];
}

function run() {
    const [text, valid] = generateRequest();
    if (!valid) {
        alert('Please fix any errors.');
        return;
    }
    $('#preview').LoadingOverlay('show');
    $.ajax({
        headers: {
            'Accept' : 'application/json',
            'Content-Type' : 'application/json'
        },
        url: 'api/classifications/v2/record/',
        type: 'POST',
        data: text,
        error: (call, status, text) => {
            $('#response_preview').text(text || `An error occurred ${status}`);
        },
        success: (record) => {
            $('#response_preview').html(formatJson(record));
            let link = null;

            if (record.results) {
                record = record.results[0];
                if (record.meta && record.meta.id) {
                    const id = record.meta.id;
                    link = `/classification/classification/${id}`;

                    $('#record-link-pair').show();
                    $('#record-link').attr('href', link);
                    $('#record-link').text(`Open classification ${id}`);
                } else {
                    $('#record-link-pair').hide();
                }
            }
        },
        complete: () => {
            $('#preview').LoadingOverlay('hide');
        }
    });
}

$(document).ready(() => {
    $('#record_id').val('test_' + Math.floor(Date.now() / 1000));
    $('#import-fields input, #import-fields select, #import-fields textarea').change(() => {
       generateRequest(); 
    }).click(() => {
        generateRequest(); 
    }).keyup(() => {
        generateRequest();
    });

    EKeys.load().then(ekeys => {
        $('.custom-key').each((index, ck) => {
            const select = $(ck);
            select.html(ekeys.keySelectOptions());
            select.chosen({
                allow_single_deselect: true,
                width: '250px',
                placeholder_text_single: `Key`
            }).change(() => {
                generateRequest();
            });
        });
        $('.custom-value').keyup(() => {generateRequest();});

        waitingForKeys = false;
        generateRequest();
    });

    /*
    EKeys.load().then(ekeys => {
        $('.custom-key').each((index, ck) => {
            let humanIndex = index+1;
            let containerDiv  = $('<label>').appendTo(ck);
            let keySelect = $('<select>', {id: `custom-key-${humanIndex}`, html: ekeys.keySelectOptions()}).appendTo(containerDiv);
            keySelect.chosen({
                allow_single_deselect: true,
                width: '300px',
                placeholder_text_single: `Key ${humanIndex}`,
            }).change(() => {
                generateRequest();
            });
            $('<input>', {id: `custom-value-${humanIndex}`, type:'text', placeholder:`value ${humanIndex}`, class:'custom-value', 'data-index':humanIndex}).keyup(() => {generateRequest()}).appendTo(ck);

        });
    });
     */
});
