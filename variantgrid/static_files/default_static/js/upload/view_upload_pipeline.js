// @ts-check
// upload/templates/upload/view_upload_pipeline.html
const COLORS = [  '#e6194b',
                '#3cb44b',
                '#ffe119',
                '#0082c8',
                '#f58231',
                '#911eb4',
                '#46f0f0',
                '#f032e6',
                '#d2f53c',
                '#fabebe',
                '#008080',
                '#e6beff',
                '#aa6e28',
                '#800000',
                '#aaffc3',
                '#808000',
                '#ffd8b1',
                '#000080',
                '#808080',
                '#000000',
];

function plotIntervalsGraph(selector, title, step_order, step_start_end_lines) {
    // There's no way to put shapes in the legend, so we'll make fake
    // traces, then only show the legend for them
    const data = [];
    const shapes = [];

    let color_index = 0;
    let y = 0;
    const y_step = -1;
    for (let i=0 ; i<step_order.length ; ++i) {
        const color = COLORS[color_index];
        ++color_index;
        const step = step_order[i];
        const step_lines = step_start_end_lines[step];

        const fake_trace = {
            x: [null],
            y: [null],
            mode: 'markers',
            name: step,
            marker: {'color': color},
            visible: 'legendonly',
        };
        data.push(fake_trace);

        const running_markers = {
            x: [],
            y: [],
            type: 'scatter',
            mode: 'markers',
            marker: {'color': color, symbol: 'triangle-right'},
            showlegend: false,
        };
        data.push(running_markers);

        const crash_markers = {
            x: [],
            y: [],
            type: 'scatter',
            mode: 'markers',
            marker: {'color': color, symbol: 'x'},
            showlegend: false,
        };
        data.push(crash_markers);

        for (let j=0 ; j<step_lines.length ; ++j) {
            const step_line = step_lines[j];

            for (let s=0 ; s<step_line.length ; s++) {
                const start_stop = step_line[s];
                const x_start = start_stop[0];
                const x_end = start_stop[1];
                const status = start_stop[2];

                const shape = {
                    x0: x_start,
                    x1: x_end,
                    y0: y,
                    y1: y,
                    type: 'line',
                    line: {
                        'color': color,
                        width: 3,
                    }
                };
                shapes.push(shape);

                if (status === 'P') { // processing - still running
                    running_markers.x.push(x_end);
                    running_markers.y.push(y);
                } else if (status === 'E') { // error
                    crash_markers.x.push(x_end);
                    crash_markers.y.push(y);
                }
            }
            y += y_step;
        }
    }


    const layout = {
        'title': {text: title},
        'xaxis': {
            showgrid: false,
            zeroline: false,
            title: {text: 'seconds'}
        },
        'yaxis': {
            showgrid: false,
            showticklabels: false,
            linewidth: 0,
            zeroline: false,
            showline: false,
        },
        'showlegend': true,
        'shapes': shapes,
    };

    Plotly.newPlot(selector, data, layout);
}

$(document).ready(function() {
    const pageData = readJsonData("view-upload-pipeline-data");
    $("button#delete-button").click(function() {
        const delete_obj_url = Urls.group_permissions_object_delete('upload.models.models.FileUpload', pageData.file_upload_id);
        $.ajax({
            type: "POST",
            url: delete_obj_url,
            success: function(data) {
                window.location = Urls.upload();
            },
            error: function(data) {
                console.log("Error: ");
                console.log(data);
                const errorMessageUl = createMessage("error", data.responseText);
                $("#delete-container").empty().append(errorMessageUl);
            }
        });

    });

    if (!$.isEmptyObject(pageData.step_start_end_lines)) {
        plotIntervalsGraph('upload-steps-graph', "Job Times", pageData.step_order, pageData.step_start_end_lines);
    }
});
