// @ts-check
// snpdb/templates/snpdb/settings/api_tokens.html, api_tokens_admin.html, api_token_detail.html
$(document).ready(() => {
    $("form.api-token-revoke").submit(function() {
        return window.confirm(`Revoke API token '${$(this).data("name")}'? Anything using it will stop working.`);
    });
    $("#copy-api-token-key").click(function() {
        const button = $(this);
        const keyInput = /** @type {HTMLInputElement} */ (document.getElementById("new-api-token-key"));
        const copied = () => button.text("Copied");
        const copyBySelection = () => {
            keyInput.select();
            if (document.execCommand("copy")) {
                copied();
            }
        };
        // navigator.clipboard only exists in a secure context (https or localhost); some deployments are plain http
        if (navigator.clipboard) {
            navigator.clipboard.writeText(keyInput.value).then(copied, copyBySelection);
        } else {
            copyBySelection();
        }
    });
});
