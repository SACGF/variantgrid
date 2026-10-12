// @ts-check
// analysis/templates/analysis/analysis_editor_and_grid.html
/* global GRID:writable, COMPONENTS_REQUIRED:writable, inFormOrLink:writable, grid_and_editor_registry:writable */
/* global layout_analysis_editor_and_grid */ // analysis.js
// Re-executed each time the fragment loads, so the state is plain assignments (a top-level let would throw)
GRID = "grid";
EDITOR = "editor";
COMPONENTS_REQUIRED = [EDITOR, GRID];
inFormOrLink = undefined;
grid_and_editor_registry = {}; // Stores by unique_code (node_id_variant_id)

function hasRequiredComponents(data) {
    let has_required_components = true;
    for (let i=0 ; i<COMPONENTS_REQUIRED.length ; ++i) {
        const component = COMPONENTS_REQUIRED[i];
        has_required_components &= component in data;
    }
    return has_required_components;
}

// Grids and editors load separately, we can make registers here so that
// we can call them when both are triggered.
function registerComponent(unique_code, name, exec_function) {
    if (!exec_function) {
        exec_function = function() { }; // Do nothing
    }
    let data = grid_and_editor_registry[unique_code];
    if (!data) {
        data = {};
        grid_and_editor_registry[unique_code] = data;
    }

    let funcs = data[name];
    if (!funcs) {
        funcs = [];
        data[name] = funcs;
    }
    //console.log("Registering code: " + unique_code + " name: " + name);
    funcs.push(exec_function);

    if (hasRequiredComponents(data)) {
        for (const k in data) {
            funcs = data[k];
            for (let i=0 ; i<funcs.length ; ++i) {
                const func = funcs[i];
                func();
            }
        }
        delete grid_and_editor_registry[unique_code];
    }
}


$(document).ready(function() {
    layout_analysis_editor_and_grid();

    $('a').on('click', function() { inFormOrLink = true; });
    $('form').on('submit', function() { inFormOrLink = true; });

    $(window).on("beforeunload", function() {
        if (!inFormOrLink) {
            getAnalysisWindow().secondWindowClosing();
        }
    });

});
