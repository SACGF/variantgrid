// @ts-check
// snpdb/templates/snpdb/data/group_permissions.html
function openDelete() {
   $("#delete-confirm-box").slideDown();
}

$(document).ready(function() {
    const pageData = readJsonData("group-permissions-data");
    $("button#delete-object").click(function() {
        const delete_obj_url = Urls.group_permissions_object_delete(pageData.class_name, pageData.instance_id); 
        $.ajax({
            type: "POST",
            url: delete_obj_url,
            success: function(data) {
                window.location = pageData.delete_redirect_url;
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
});
