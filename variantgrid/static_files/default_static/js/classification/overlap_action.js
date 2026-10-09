// @ts-check
// classification/templates/classification/overlap_action.html
function updateCalculator(valueType) {
    const selectedValues = $('.value-change option:selected').map(function() {
        return this.value;
    }).get();
    const nonInteractiveValues = $('.non-interactive').map(function() {
        return this.value;
    }).get();

    const uniqueValues = [...new Set(selectedValues.concat(nonInteractiveValues))];
    const valueStr = uniqueValues.join(',');
    const baseUrl = `${Urls.overlap_calc()}?value_type=${valueType}`;
    const ajaxDiv = $('<div>', {'data-toggle': 'ajax', 'href': `${baseUrl}&values=${valueStr}`});
    const parentDiv = $('#resulting-status');
    parentDiv.html(ajaxDiv);
}

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

function initOverlapAction(valueType) {
    $('.value-change').change(() => updateCalculator(valueType));
    updateCalculator(valueType);

    $("input[name='outcome']").change(showFormCheck);
    showFormCheck();
}
