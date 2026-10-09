// @ts-check
// classification/templates/classification/overlap_action.html
function updateCalculator() {
    const selectedValues = $('.value-change option:selected').map(function() {
        return this.value;
    }).get();
    const nonInteractiveValues = $('.non-interactive').map(function() {
        return this.value;
    }).get();

    const uniqueValues = [...new Set(selectedValues.concat(nonInteractiveValues))];
    const valueStr = uniqueValues.join(',');
    const valueType = readJsonData("overlap-action-data").value_type;
    const baseUrl = `${Urls.overlap_calc()}?value_type=${valueType}`;
    const ajaxDiv = $('<div>', {'data-toggle': 'ajax', 'href': `${baseUrl}&values=${valueStr}`});
    const parentDiv = $('#resulting-status');
    parentDiv.html(ajaxDiv);
}

$(document).ready(() => {
    $('.value-change').change(() => {updateCalculator();});
    updateCalculator();
});

function showFormCheck() {
    const value = $("input[name='outcome']:checked").val();
    const allValues = {
        "#pending-changes": value == "agree",
        "#postponed-changes": value == "postpone"
    };
    Object.entries(allValues).forEach(entry => {
        const [key, value] = entry;
        if (value) {
            $(key).slideDown();
        } else {
            $(key).slideUp();
        }
    });
}

$(document).ready(() => {
   $("input[name='outcome']").change(showFormCheck);
   showFormCheck();
});
