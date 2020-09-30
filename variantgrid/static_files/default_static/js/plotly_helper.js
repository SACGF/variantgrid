function defaultFor(arg, val) {
	return typeof arg !== 'undefined' ? arg : val;
}

function defaultLayout(title, width, height) {
    var layout = {
      title: title,
      'width': defaultFor(width, 600),
      'height': defaultFor(height, 400),
      xaxis: { autotick: false },
    };
    return layout;
}


function dictToData(dict, type) {
    var labels = [];
    var values = [];
    for (var k in dict) {
        labels.push(k);
        values.push(dict[k]);
    }

    var data = {
      values: values,
      labels: labels,
      type: type,
    };
    return data;    
}


function plotPieDict(selector, title, dict, width, height) {
    var data = dictToData(dict, 'pie');
    var layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, [data], layout);
}


function plotBoxDict(selector, title, dict, width, height) {
    var keys = Object.keys(dict);
    return plotBoxKeysDict(selector, title, keys, dict, width, height);
}
    

function getBoxDataFromDict(keys, dict) {
    var data = [];

    for (var i=0 ; i<keys.length ; i++) {
        var k = keys[i];
        var d = {
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
    var data = getBoxDataFromDict(keys, dict);
    var layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, data, layout);
}

function getBoxDataFromLabelsMatrix(labels, matrix) {
    var data = [];

    for (var i=0 ; i<matrix.length ; i++) {
        var label = labels[i];
        var d = {
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
    var data = getBoxDataFromLabelsMatrix(labels, matrix);
    var layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, data, layout);
}

function plotBarDict(selector, title, dict, width, height) {
    var data = dictToData(dict, 'bar');

    var layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, [data], layout);
}


function plotBarArrays(selector, title, x, y, width, height, color) {
    var data = {
        x: x,
        y: y,
        type: 'bar',
    };
    
    if (color) {
        data["marker"] = {color: color};
    }
    var layout = defaultLayout(title, width, height);
    Plotly.newPlot(selector, [data], layout);
}



function plotHBarArrays(selector, title, x, y, width, height, color, margin) {
    var data = {
        x: x,
        y: y,
        type: 'bar',
        orientation: 'h',
    };
    
    if (color) {
        data["marker"] = {color: color};
    }
    var layout = defaultLayout(title, width, height);
    layout["xaxis"] = {autotick: true};

    if (margin) {
        layout["autosize"] = false;
        layout["margin"] = margin;
    }
    Plotly.newPlot(selector, [data], layout);
}


function plotLineArrays(selector, x, y, layout) {
    var data = [
      {
        x: x,
        y: y,
        type: 'line',
      }
    ];

    Plotly.newPlot(selector, data, layout);
}


function showStackedBar(elementId, title, named_data, x_labels) {
    var x = [];
    for (var i=0 ; i<x_labels.length ; ++i) {
        x.push(i);
    } 
    
    var data = [];
    for(var i=0 ; i<named_data.length ; ++i) {
        var nd = named_data[i];
        data.push({ 'x' : x,
                    'y' : nd[1],
                    'name' : nd[0],
                    'type' : 'bar'});
    }

    var layout = {
      xaxis: {  tickvals: x,
                ticktext: x_labels,
                tickmode: 'array',},
      barmode: 'relative',
      title: title
    };

    $("#" + elementId).empty();
    Plotly.newPlot(elementId, data, layout);
}


function showHeatMap(elementId, title, x, y, z, labels) {
    var data = [{
      x: x,
      y: y,
      z: z,
      type: 'heatmap',
    }];
    
    var layout = {
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
    
    for ( var i = 0; i < y.length; i++ ) {
      for ( var j = 0; j < x.length; j++ ) {
        var currentValue = labels[i][j];
        if (currentValue != 0.0) {
          var textColor = 'white';
        }else{
          var textColor = 'black';
        }
        var result = {
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
    var TITLE_Y = -20;
    var textTitle = $("text.gtitle", "#" + elementId);
    var y = parseInt(textTitle.attr("y"));
    textTitle.attr("y", y + TITLE_Y);
}



