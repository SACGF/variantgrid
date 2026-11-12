// @ts-check
// snpdb/templates/snpdb/settings/api_tokens.html, api_tokens_admin.html, api_token_detail.html
$(document).ready(() => {
    $("form.api-token-revoke").submit(function() {
        return window.confirm(`Revoke API token '${$(this).data("name")}'? Anything using it will stop working.`);
    });
    $("#copy-api-token-key").click(function() {
        navigator.clipboard.writeText(String($("#new-api-token-key").val()));
        $(this).text("Copied");
    });
});
