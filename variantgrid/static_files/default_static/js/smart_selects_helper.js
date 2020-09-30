// Work around for issue https://github.com/digi604/django-smart-selects/issues/219
function initSmartSelectsInForm(formSelect) {
    function initItem(item) {
        var chainfield = "#id_" + $(item).attr("data-chainfield");
        var url = $(item).attr("data-url");
        var id = "#" + $(item).attr("id");
        var value = JSON.parse($(item).attr("data-value"));
        var auto_choose = $(item).attr("data-auto_choose");
        if($(item).hasClass("chained-fk")) {
            var empty_label = $(item).attr("data-empty_label");
            chainedfk.init(chainfield, url, id, value, empty_label, auto_choose);
        } else if ($(item).hasClass("chained")) {
            chainedm2m.init(chainfield, url, id, value, auto_choose);
        }
    }

    $.each($(".chained-fk", formSelect), function(index, item) {
        initItem(item);
    });
}