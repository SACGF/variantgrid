// @ts-check
// upload/templates/upload/upload.html
const RUNNING_STATUSES = ['C', 'P'];
const INITIAL_RELOAD_SECS = 1;
const MAX_RELOAD_SECS = 16;
let reloadSecs = INITIAL_RELOAD_SECS;
let reloadTimer = null;

function fileUploadDatatableFilter(data) {
    data.file_type = $("#file-type-select").val();
}

function reloadUploadsGrid() {
    clearTimeout(reloadTimer);
    reloadTimer = null;
    $("#file-uploads-datatable").DataTable().ajax.reload(null, false);
}

// Keep reloading, backing off, while any row on the page is still importing
function scheduleReloadWhileRunning(json) {
    const running = (json.data || []).some(row => row.status && RUNNING_STATUSES.includes(row.status.status));
    if (!running) {
        reloadSecs = INITIAL_RELOAD_SECS;
        return;
    }
    if (!reloadTimer) {
        reloadTimer = setTimeout(reloadUploadsGrid, reloadSecs * 1000);
        reloadSecs = Math.min(reloadSecs * 2, MAX_RELOAD_SECS);
    }
}

function renderUploadStatus(data, type, row) {
    if (!data || !data.icon) {
        return '';
    }
    let icon = $('<i>', {class: `fas ${data.icon} ${data.css || ''}`, title: data.title});
    if (data.url) {
        icon = $('<a>', {href: data.url}).append(icon);
    }
    if (data.summary) {
        return $('<span>', {title: data.title}).append(icon, $('<span>', {class: 'text-danger small ml-1', text: data.summary}));
    }
    return icon;
}

function renderUploadFileType(data, type, row) {
    const dom = $('<span>', {class: 'upload-file-type'}).append(data.icon || '');
    if (data.label) {
        if (data.url) {
            dom.append($('<a>', {href: data.url, class: 'hover-link', text: data.label}));
        } else {
            dom.append($('<span>', {text: data.label}));
        }
    }
    return dom;
}

$(document).ready(function() {
    $("#file-uploads-datatable").on('xhr.dt', function(e, settings, json) {
        if (json) {
            scheduleReloadWhileRunning(json);
        }
    });
    $("#file-type-select").change(reloadUploadsGrid);

    const userAgent = navigator.userAgent.toLowerCase();
    const isFirefox = userAgent.indexOf('firefox') > -1;
    const isLinux = userAgent.indexOf('linux') > -1;
    if (isFirefox && isLinux) {
        $("#linux-firefox-warning").show();
    }
});
