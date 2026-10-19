function defaultFor(arg, val) {
	return typeof arg !== 'undefined' ? arg : val;
}

/** Pass width=null for a plot that fills (and resizes with) its container */
function defaultLayout(title, width, height) {
    const layout = {
      title: {text: title},
      'height': defaultFor(height, 400),
      xaxis: { tickmode: 'linear' },
    };
    if (width === null) {
        layout.autosize = true;
    } else {
        layout.width = defaultFor(width, 600);
    }
    return layout;
}


function dictToData(dict, type) {
    const labels = [];
    const values = [];
    for (const k in dict) {
        labels.push(k);
        values.push(dict[k]);
    }

    const data = {
      values: values,
      labels: labels,
      type: type,
    };
    return data;    
}


function plotPieDict(selector, title, dict, width, height) {
    const data = dictToData(dict, 'pie');
    const layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, [data], layout);
}


function plotBoxDict(selector, title, dict, width, height) {
    const keys = Object.keys(dict);
    return plotBoxKeysDict(selector, title, keys, dict, width, height);
}
    

function getBoxDataFromDict(keys, dict) {
    const data = [];

    for (let i=0 ; i<keys.length ; i++) {
        const k = keys[i];
        const d = {
            name: k,
            type: 'box',
            y: dict[k],
            boxpoints: 'all',
            pointpos: 0,
            jitter: 0.3,
        };
        data.push(d);
    }
    return data;
}

function plotBoxKeysDict(selector, title, keys, dict, width, height) {
    const data = getBoxDataFromDict(keys, dict);
    const layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, data, layout);
}

function getBoxDataFromLabelsMatrix(labels, matrix) {
    const data = [];

    for (let i=0 ; i<matrix.length ; i++) {
        const label = labels[i];
        const d = {
            name: label,
            type: 'box',
            y: matrix[i],
            boxpoints: 'all',
            pointpos: 0,
            jitter: 0.3,
        };
        data.push(d);
    }
    return data;
}

function plotBarArrays(selector, title, x, y, width, height, color) {
    const data = {
        x: x,
        y: y,
        type: 'bar',
    };
    
    if (color) {
        data["marker"] = {color: color};
    }
    const layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, [data], layout);
}



function plotHBarArrays(selector, title, x, y, width, height, color, margin) {
    const data = {
        x: x,
        y: y,
        type: 'bar',
        orientation: 'h',
    };
    
    if (color) {
        data["marker"] = {color: color};
    }
    const layout = defaultLayout(title, width, height);
    layout["xaxis"] = {tickmode: 'auto'};

    if (margin) {
        layout["margin"] = margin;
    }
    Plotly.newPlot(selector, [data], layout, {responsive: layout.autosize === true});
}


// series_colors (optional) maps a series name to its colour; unnamed series keep Plotly's defaults
function showStackedBar(elementId, title, width, height, named_data, x_labels, series_colors) {
    const x = [];
    for (let i=0 ; i<x_labels.length ; ++i) {
        x.push(i);
    } 
    
    const data = [];
    for(let i=0 ; i<named_data.length ; ++i) {
        const nd = named_data[i];
        const trace = { 'x' : x,
                        'y' : nd[1],
                        'name' : nd[0],
                        'type' : 'bar'};
        const color = series_colors && series_colors[nd[0]];
        if (color) {
            trace.marker = {color: color};
        }
        data.push(trace);
    }

    const layout = defaultLayout(title, width, height);
    layout.xaxis = {
        tickvals: x,
        ticktext: x_labels,
        tickmode: 'array',
    };
    layout.barmode = 'relative';

    $("#" + elementId).empty();
    Plotly.newPlot(elementId, data, layout);
}
