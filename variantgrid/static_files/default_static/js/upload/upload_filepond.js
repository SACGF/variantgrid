// @ts-check
// upload/templates/upload/upload.html
/* global getCookie */ // global.js
/* global FilePond */ // lib/filepond
/* global reloadUploadsGrid */ // upload/upload.js
(function() {
    const inputElement = document.querySelector('#upload-file');
    if (!inputElement) {
        return;
    }
    const csrftoken = getCookie('csrftoken');

    const pond = FilePond.create(inputElement, {
        allowMultiple: true,
        credits: false,
        dropOnPage: true,
        dropOnElement: false,
        server: {
            process: {
                url: Urls.upload_file(),
                method: 'POST',
                headers: {'X-CSRFToken': csrftoken},
            },
            revert: (uniqueFileId, load, error) => {
                fetch("/upload/upload_file_delete/" + uniqueFileId, {
                    method: 'DELETE',
                    headers: {'X-CSRFToken': csrftoken},
                }).then((resp) => {
                    if (resp.ok) {
                        load();
                        reloadUploadsGrid();
                    } else {
                        error('Could not revert');
                    }
                }).catch(() => error('Could not revert'));
            },
        },
    });

    function progressRowId(fileId) {
        return 'upload-progress-' + fileId.replace(/[^a-zA-Z0-9_-]/g, '');
    }

    pond.on('addfile', (error, file) => {
        if (error) {
            return;
        }
        const tr = $('<tr/>', {id: progressRowId(file.id)});
        tr.append('<td><i class="fas fa-spinner fa-spin" title="Uploading"></i></td>');
        const nameTd = $('<td/>').append($('<div/>').text(file.filename));
        nameTd.append($('<div/>').addClass('progress')
            .append($('<div/>').addClass('progress-bar bg-success').css('width', '0%')));
        tr.append(nameTd);
        $("#upload-progress-table tbody").prepend(tr);
    });

    pond.on('processfileprogress', (file, progress) => {
        $('#' + progressRowId(file.id)).find('.progress-bar').css('width', (progress * 100) + '%');
    });

    pond.on('processfile', (error, file) => {
        $('#' + progressRowId(file.id)).remove();
        if (!error) {
            pond.removeFile(file.id);
            reloadUploadsGrid();
        }
    });

    pond.on('processfileabort', (file) => {
        $('#' + progressRowId(file.id)).remove();
    });

    const browseButton = document.querySelector('#upload-file-browse');
    if (browseButton) {
        browseButton.addEventListener('click', () => pond.browse());
    }
})();
