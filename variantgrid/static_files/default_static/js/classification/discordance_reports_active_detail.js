// @ts-check
// classification/templates/classification/discordance_reports_active_detail.html
$(document).ready(() => {
    $('input[name=discordances-filter]').click(() => {
        const selected = $("input[name=discordances-filter]:checked").val();

        $('.contact-details').hide();
        $(`.contact-details[data-lab=${selected}]`).show();
        $('.discordance-row').each((index, row) => {
            row = $(row);
            const data_labs = row.attr('data-labs');
            let show = false;
            if (selected == "all" || !data_labs) {
                show = true;
            } else {
                if (selected == "internal") {
                    show = data_labs == 'internal';
                } else {
                    const labs = data_labs.split(";");
                    show = labs.includes(selected);
                }
            }
            if (show) {
                row.show();
            } else {
                row.hide();
            }
        });
    });
});
