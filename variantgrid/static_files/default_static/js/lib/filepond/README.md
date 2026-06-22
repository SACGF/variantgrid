# FilePond vendored assets

Vendored copies of FilePond and a single plugin, served directly from
`{% static %}`. Each library lives in its own version subdirectory so the
version is visible in the path (matches the `bootstrap/4.5.2`,
`codemirror/5.48.4`, `datatables/1.11.4` layout used elsewhere in this
tree).

| Package | Version | Path | Source |
| --- | --- | --- | --- |
| `filepond` (core JS + CSS) | 4.32.12 | `4.32.12/filepond.min.{js,css}` | https://unpkg.com/filepond@4.32.12/dist/ |
| `filepond-plugin-file-poster` (thumbnail poster for image attachments) | 2.5.2 | `file-poster/2.5.2/filepond-plugin-file-poster.min.{js,css}` | https://unpkg.com/filepond-plugin-file-poster@2.5.2/dist/ |

When upgrading, drop a new versioned folder beside the existing one, retarget
the `{% static %}` paths in the three consumer templates
(`upload/templates/upload/upload.html`,
`patients/templates/patients/view_patient.html`,
`classification/templates/classification/classification.html`), update this
table, then delete the old folder.

`filepond-overrides.css` (next to this README) is our own — local styling
on top of FilePond — and is not versioned.

Replaces the previous `django-jfu` + blueimp/jQuery-File-Upload stack
(issue #1566).
