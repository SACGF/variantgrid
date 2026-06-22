function defaultFor(arg, val) {
	return typeof arg !== 'undefined' ? arg : val;
}

function defaultLayout(title, width, height) {
    const layout = {
      title: title,
      'width': defaultFor(width, 600),
      'height': defaultFor(height, 400),
      xaxis: { autotick: false },
    };
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

function plotBoxLabelsMatrix(selector, title, labels, matrix, width, height) {
    const data = getBoxDataFromLabelsMatrix(labels, matrix);
    const layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, data, layout);
}

function plotBarDict(selector, title, dict, width, height) {
    const data = dictToData(dict, 'bar');

    const layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, [data], layout);
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
    layout["xaxis"] = {autotick: true};

    if (margin) {
        layout["autosize"] = false;
        layout["margin"] = margin;
    }
    Plotly.newPlot(selector, [data], layout);
}


function plotLineArrays(selector, x, y, layout) {
    const data = [
      {
        x: x,
        y: y,
        type: 'line',
      }
    ];

    Plotly.newPlot(selector, data, layout);
}


function showStackedBar(elementId, title, width, height, named_data, x_labels) {
    const x = [];
    for (let i=0 ; i<x_labels.length ; ++i) {
        x.push(i);
    } 
    
    const data = [];
    for(let i=0 ; i<named_data.length ; ++i) {
        const nd = named_data[i];
        data.push({ 'x' : x,
                    'y' : nd[1],
                    'name' : nd[0],
                    'type' : 'bar'});
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


function showHeatMap(elementId, title, x, y, z, labels) {
    const data = [{
      x: x,
      y: y,
      z: z,
      type: 'heatmap',
    }];
    
    const layout = {
      title: title,
      titlefont: {
        size: 32,
      },
      annotations: [],
      xaxis: {
        title: 'Old Category',
        ticks: '',
        side: 'top'
      },
      yaxis: {
        title: 'New Category',
        ticks: '',
        ticksuffix: ' ',
        width: 700,
        height: 700,
        autosize: false
      }
    };

    let textColor;
    for ( let i = 0; i < y.length; i++ ) {
      for ( let j = 0; j < x.length; j++ ) {
        const currentValue = labels[i][j];
        if (currentValue != 0.0) {
          textColor = 'white';
        } else {
          textColor = 'black';
        }
        const result = {
          xref: 'x1',
          yref: 'y1',
          x: x[j],
          y: y[i],
          text: currentValue + " (" + (z[i][j]).toFixed(2) + "%)",
          font: {
            family: 'Arial',
            size: 12,
            // color: 'rgb(50, 171, 96)'
          },
          showarrow: false,
          font: {
            color: textColor
          }
        };
        layout.annotations.push(result);
      }
    }
    
    Plotly.newPlot(elementId, data, layout);
    
    // Shift title up a bit
    const TITLE_Y = -20;
    const textTitle = $("text.gtitle", "#" + elementId);
    y = parseInt(textTitle.attr("y"));
    textTitle.attr("y", y + TITLE_Y);
}



