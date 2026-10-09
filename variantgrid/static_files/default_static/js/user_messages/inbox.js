// @ts-check
// user_messages/templates/user_messages/inbox.html
function deleteSelectedMessages() {
    const selectedMessages = [];
    $("input.message-checkbox:checked").each(function() {
        selectedMessages.push($(this).attr("message_id"));
    });
    if (selectedMessages) {
        const data = 'message_ids=' + JSON.stringify(selectedMessages);;

        $.ajax({
            type: "POST",
            data: data,
            url: Urls.messages_bulk_delete(),
            success: function(data) {
                window.location.reload();
            }
        });
    }
}

$(document).ready(function() {
    $("#select-all-checkbox").click(function() {
        const checked = $(this).is(":checked");
        $("input.message-checkbox").each(function() {
            $(this).prop('checked', checked);
        });
    });
});
