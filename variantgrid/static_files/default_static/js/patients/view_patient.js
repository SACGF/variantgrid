// @ts-check
// patients/templates/patients/view_patient.html
/* global getCookie */ // global.js
/* global FilePond, FilePondPluginFilePoster */ // lib/filepond
(function() {
    const inputElement = document.querySelector('#patient-file-upload');
    if (!inputElement) {
        return;
    }
    if (typeof FilePondPluginFilePoster !== 'undefined') {
        FilePond.registerPlugin(FilePondPluginFilePoster);
    }
    const data = readJsonData("view-patient-data");
    const csrftoken = getCookie('csrftoken');
    const pond = FilePond.create(inputElement, {
        allowMultiple: true,
        credits: false,
        server: {
            process: {
                url: Urls.patient_file_upload(data.patient_id),
                method: 'POST',
                headers: {'X-CSRFToken': csrftoken},
            },
            revert: (uniqueFileId, load, error) => {
                fetch("/patients/patient_file_delete/" + uniqueFileId, {
                    method: 'DELETE',
                    headers: {'X-CSRFToken': csrftoken},
                }).then((resp) => {
                    if (resp.ok) {
                        load();
                    } else {
                        error('Could not revert');
                    }
                }).catch(() => error('Could not revert'));
            },
        },
    });

    const browseButton = document.querySelector('#patient-file-browse');
    if (browseButton) {
        browseButton.addEventListener('click', () => pond.browse());
    }
})();
