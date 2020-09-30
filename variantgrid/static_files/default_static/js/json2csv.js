// From https://stackoverflow.com/a/24643992/295724
function JSON2CSV(objArray, header, quote) {
    var getValue = function(val, defaultValue) {
        return (typeof val !== 'undefined') ?  val : defaultValue;
    }

    var header = getValue(header, false);
    var quote = getValue(quote, false);

    var array = typeof objArray != 'object' ? JSON.parse(objArray) : objArray;
    var str = '';
    var line = '';

    if (header) {
        var head = array[0];
        if (quote) {
            for (var index in array[0]) {
                var value = index + "";
                line += '"' + value.replace(/"/g, '""') + '",';
            }
        } else {
            for (var index in array[0]) {
                line += index + ',';
            }
        }

        line = line.slice(0, -1);
        str += line + '\r\n';
    }

    for (var i = 0; i < array.length; i++) {
        var line = '';

        if (quote) {
            for (var index in array[i]) {
                var value = array[i][index] + "";
                line += '"' + value.replace(/"/g, '""') + '",';
            }
        } else {
            for (var index in array[i]) {
                line += array[i][index] + ',';
            }
        }

        line = line.slice(0, -1);
        str += line + '\r\n';
    }
    return str;
}