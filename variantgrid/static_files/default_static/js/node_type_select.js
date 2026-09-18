// Node class icons and the themed node type dropdown - the analysis add-node toolbar and the analysis issues filter.
// Styles are in analysis_nodes.css, icon symbols in uicore/tags/svg_icon_sprite.html

// icon is a NodeIcon dict: FontAwesome classes, or a symbol id in svg_icon_sprite.html
function renderNodeIcon(icon) {
    if (icon && icon.symbol) {
        const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("class", "node-icon");
        const use = document.createElementNS("http://www.w3.org/2000/svg", "use");
        use.setAttribute("href", "#" + icon.symbol);
        svg.appendChild(use);
        return $(svg);
    }
    return $("<i/>", {class: "node-icon " + ((icon && icon.fa) || "")});
}

// Bootstrap dropdown over a <select> of node class names - callers still read/listen to the select.
// nodeTypes is node_types.get_node_display_data_by_class_name() - icon and source/filter colour
function setupNodeTypeSelect(select, nodeTypes) {
    function renderNodeTypeItem(className, label) {
        const nodeType = nodeTypes[className];
        // Class name on the row picks up the node's accent colour - see analysis_nodes.css
        const wrapper = $("<div>", {"class": "node-type-item " + ((nodeType && nodeType.class_name) || "")});
        if (nodeType) {
            wrapper.attr("node_classification", nodeType.classification);
        }
        renderNodeIcon(nodeType && nodeType.icon).appendTo(wrapper);
        $("<span>", {text: label}).appendTo(wrapper);
        return wrapper;
    }

    select = $(select);
    if (!select.length) {
        return;
    }
    select.hide();
    const buttonId = select.attr("id") + "-button";
    const button = $("<button>", {id: buttonId, type: "button", "class": "dropdown-toggle node-type-button",
                                  "data-toggle": "dropdown", "aria-haspopup": "true", "aria-expanded": "false"});
    const menu = $("<div>", {"class": "dropdown-menu node-type-menu", "aria-labelledby": buttonId});

    function addMenuItem(option) {
        $("<a>", {"class": "dropdown-item", href: "javascript:void(0)", "data-value": option.val()})
            .append(renderNodeTypeItem(option.val(), option.text()))
            .appendTo(menu);
    }

    select.children().each(function() {
        const child = $(this);
        if (child.is("optgroup")) {
            $("<h6>", {"class": "dropdown-header", text: child.attr("label")}).appendTo(menu);
            child.children("option").each(function() { addMenuItem($(this)); });
        } else {
            addMenuItem(child);
        }
    });

    function showSelected() {
        const value = select.val();
        const label = select.find("option:selected").text();
        button.empty().append(renderNodeTypeItem(value, label).addClass("node-type-button-text"));
        $(".dropdown-item", menu).removeClass("active")
            .filter("[data-value='" + value + "']").addClass("active");
    }

    menu.on("click", ".dropdown-item", function() {
        select.val($(this).data("value")).trigger("change");
    });
    select.on("change.nodeTypeSelect", showSelected);  // namespaced so a caller can refresh the button alone

    $("<div>", {"class": "dropdown"}).append(button, menu).insertAfter(select);
    showSelected();
}
