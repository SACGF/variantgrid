// @ts-check
// snpdb/templates/snpdb/data/group_permissions.html
function openDelete() {
   $("#delete-confirm-box").slideDown();
}

function initGroupPermissions(className, instanceId, deleteRedirectUrl) {
    $("button#delete-object").click(function() {
        const delete_obj_url = Urls.group_permissions_object_delete(className, instanceId);
        $.ajax({
            type: "POST",
            url: delete_obj_url,
            success: function(data) {
                window.location = deleteRedirectUrl;
            },
            error: function(data) {
                const errorMessageUl = createMessage("error", data.responseText);
                $("#delete-container").empty().append(errorMessageUl);
            }
        });    

    });


    $("button#no-delete").click(function() {
       $("#delete-confirm-box").slideUp();
    });

    const options = {
        target: '#permissions-embedded-page'
    };
    $('form#group-permission-form').ajaxForm(options);
}
