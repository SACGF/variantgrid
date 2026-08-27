const endpointColor = "#121212";
const ACTIVE_CLASS = "ui-selected";
const DELETING_CLASS = "node-deleting";  // fade out transition, see analysis_nodes.css
const ACTIVE_NODE_COUNT_CLASS = "node-counts-selected"; // Needs to be different so not grabbed with multi-draggable
const NODE_COUNT_TOTAL = "T";  // snpdb.models.models_enums.BuiltInFilters.TOTAL
const SHOW_NODE_IDS_IN_TOOLTIPS = true;

// Connector/endpoint look - shared between the instance defaults and the per-node endpoints below
const CONNECTOR_STYLE = {stroke: endpointColor, strokeWidth: 3};
const CONNECTOR_SPEC = {type: "Bezier", options: {curviness: 63}};
const ENDPOINT_STYLE = {fill: endpointColor};
const DOT_ENDPOINT = {type: "Dot", options: {radius: 11}};

// Which way the graph flows. Everything else about the endpoints (uuids, styles, MergeNode multi-input,
// invalid-target highlighting) is orientation-independent - see setupConnections
const NODE_ORIENTATIONS = {
	vertical: {
		instanceAnchors: ["Bottom", "Top"],
		input: "Top",
		output: "Bottom",
		// Venn takes 2 inputs, spread along the input edge
		vennInputs: [[0.25, 0, 0, -1], [0.75, 0, 0, -1]],
	},
	horizontal: {
		instanceAnchors: ["Right", "Left"],
		input: "Left",
		output: "Right",
		// The card's rings are rotated anti-clockwise (see analysis_nodes.css), which puts the left
		// ring at the bottom - so the left input takes the lower endpoint
		vennInputs: [[0, 0.75, -1, 0], [0, 0.25, -1, 0]],
	},
};

function isHorizontalNodeFlow() {
	return typeof ANALYSIS_HORIZONTAL_MODE !== "undefined" && ANALYSIS_HORIZONTAL_MODE;
}

function getNodeOrientation() {
	return NODE_ORIENTATIONS[isHorizontalNodeFlow() ? "horizontal" : "vertical"];
}

// @jsplumb/browser-ui instance, created in variantgridPipeline.init once #analysis is in the DOM
let jsPlumbInstance = null;

function getNode(nodeId) {
	return $("#analysis-node-" + nodeId);
}

// Nodes should have a div of class "node-overlay" that has its CSS set by nodeData.overlay_css_classes
// Nodes should also create a "updateState(args)" method, this will be called after creation
//
// Everything else on the card - badge icon, class strip, chips - comes from the rendering dict, see
// analysis/models/nodes/node_display.py

// icon is a NodeIcon dict: FontAwesome classes, or a symbol id in node_icon_sprite.html
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

function renderChip(chip) {
	const span = $("<span/>", {class: "node-chip"});
	if (chip.css_class) {
		span.addClass(chip.css_class);
	}
	if (chip.title) {
		span.attr("title", chip.title);
	}
	if (chip.icon) {
		$("<i/>", {class: chip.icon}).appendTo(span);
	}
	$("<span/>", {class: "node-chip-text", text: chip.text}).appendTo(span);
	return span;
}

function createDefaultNode() {
	const div = $('<div/>').addClass("window design-a-node");
	const nodeOverlay = $('<div/>').addClass("node-overlay");
	$("<span/>", {class: "node-badge node-icon-badge"}).appendTo(nodeOverlay);
	$("<span/>", {class: "node-klass"}).appendTo(nodeOverlay);
	$("<div/>", {class: "node-name-holder"}).append($("<span/>", {class: "node-name"})).appendTo(nodeOverlay);
	$("<div/>", {class: "node-chips"}).appendTo(nodeOverlay);
	$("<div />", {class: "node-color-overlay"}).appendTo(nodeOverlay);
	div.append(nodeOverlay);
	$("<div/>", {class: "node-counts-strip"}).appendTo(div);
	div[0].updateState = function(args) { };
	return div;
}

const VENN_WIDTH = 56;
const VENN_HEIGHT = 40;

function createVennNode() {
	// The rings say "venn" better than a badge and class strip can, so the card carries neither
	const div = createDefaultNode();
	$(".node-badge, .node-klass", div).remove();
	const vennHolder = $('<div/>', {class: "node-venn"});
	venn2(vennHolder[0], VENN_WIDTH, VENN_HEIGHT);
	// venn2 leaves a wide margin around the rings - crop to them so the card can shrink to the widget
	$("svg", vennHolder).attr("viewBox", "4 4 44 32");
	vennHolder.insertBefore($(".node-counts-strip", div));

	div[0].updateState = function(args) {
		venn_select(this, args['venn_flag']);
	};
	return div;
}

function createMergeNode() {
	// Merging is all the node does, so the glyph in the middle says it better than a name and class strip
	const div = createDefaultNode();
	$(".node-badge, .node-klass", div).remove();
	const graphic = $('<div/>', {class: "node-graphic"});
	graphic.html('<svg viewBox="0 0 32 30"><use href="#node-icon-merge"></use></svg>');
	graphic.insertBefore($(".node-counts-strip", div));
	return div;
}

function createNodeFromData(nodeData) {
	const NODE_FACTORIES = {
		"VennNode" : createVennNode,
		"MergeNode" : createMergeNode,
	};

	const nodeClass = nodeData['node_class'];
	let factory = NODE_FACTORIES[nodeClass];
	if (!factory) {
		factory = createDefaultNode;
	}
	const node = factory();
	updateNodeFromData(node, nodeData);

	return node;
}

// The name is clamped to 3 lines on the card, so hovering has to give the whole thing
function updateNodeTitle(node) {
	const parts = [$(".node-name", node).text()];
	let nodeHelp = NODE_HELP[node.attr("node_class")];
	if (nodeHelp) {
		if (SHOW_NODE_IDS_IN_TOOLTIPS) {
			nodeHelp += " (node #" + node.attr("node_id") + ")";
		}
		parts.push(nodeHelp);
	}
	node.attr("title", parts.filter(Boolean).join("\n\n"));
}

function isManagedByJsPlumb(el) {
	return jsPlumbInstance !== null && el.hasAttribute(jsPlumb.ATTRIBUTE_MANAGED);
}

// Chips and count rows change the card height, so jsPlumb needs fresh offsets for the Top/Bottom anchors
function repaintNode(node) {
	node.each(function() {
		if (isManagedByJsPlumb(this)) {
			jsPlumbInstance.revalidate(this);
		}
	});
}

function updateNodeFromData(node, nodeData) {
	node.addClass(nodeData.node_class);
	node.attr(nodeData.attributes);
	const nodeOverlay = $(".node-overlay", node);
	nodeOverlay.attr("class", nodeData.overlay_css_classes);
	$(".node-name", node).text(nodeData.name);
	$(".node-klass", node).text(nodeData.class_label_short);
	$(".node-badge", node).empty().append(renderNodeIcon(nodeData.icon));

	const chipsHolder = $(".node-chips", node).empty();
	$.each(nodeData.chips || [], function() {
		chipsHolder.append(renderChip(this));
	});

	updateNodeTitle(node);
	node.each(function() { this.updateState(nodeData['args']); });
	repaintNode(node);
}

// Wait for previous update to come back (so we don't end up with race conditions on the server)
let waiting_for_message_callback = false;
const update_message_buffer = [];

function updateNode(nodeId, op, params, on_success_function) {
	const message = [nodeId, op, params, on_success_function];
	update_message_buffer.push(message);
	checkSendMessage();
}


function checkSendMessage() {
	if (!waiting_for_message_callback) {
		const message = update_message_buffer.shift();
		if (message) {
			const nodeId = message[0];
			const op = message[1];
			const params = message[2];
			const old_on_success_function = message[3];
			const on_success_function = function (data) {
				if (old_on_success_function) {
					old_on_success_function(data);
				}
				waiting_for_message_callback = false;
				checkSendMessage();
			};
			sendUpdateNodeMessage(nodeId, op, params, on_success_function);
		}
	}
}

function sendUpdateNodeMessage(nodeId, op, params, on_success_function) {
	waiting_for_message_callback = true;
	const data = '&op=' + op + '&params=' + JSON.stringify(params);
	$.ajax({
	    type: "POST",
	    data: data,
	    url: Urls.node_update(ANALYSIS_ID, nodeId),
	    success: function(data) {
	       on_success_function(data.dirty_nodes);
        },
	});
}

function addNodesToDOM(selector, nodeDataArray, readOnly) {
	selector = $(selector);
	for (let i=0 ; i<nodeDataArray.length ; ++i) {
		const node = createNodeFromData(nodeDataArray[i]);
		selector.append(node);
	}
}


function getEndpoint(id, endpoint_type, side) {
	const el = document.getElementById(id);
	if (!el)
		return null;

	const endpoints = jsPlumbInstance.getEndpoints(el);
	if (!endpoints)
        return null;    

	for (let i=0 ; i<endpoints.length ; ++i) {
		const ep = endpoints[i];
		if (endpoint_type === 'source') {
			if (ep.isSource) {
				return ep;
			}
		} else if (endpoint_type === 'target') {
			if (ep.isTarget) {
				if (side) {
					if (side !== getEndpointSide(ep)) {
						continue;
					}
				}
				return ep;
			}
		}
	}		 
}

function setupHideInvalidConnectionsOnDrag() {
		// Endpoints we turned red for the connection drag currently in progress
		const invalidEndpoints = [];

		function setFillStyle(endpoint, color) {
			endpoint.setPaintStyle(Object.assign({}, endpoint.getPaintStyle(), {fill: color}));
			jsPlumbInstance.revalidate(endpoint.element);
		}

		function setEndpointInvalid(endpoint) {
		    endpoint.enabled = false;
		    setFillStyle(endpoint, "#FF0000");
			invalidEndpoints.push(endpoint);
		}

		function resetInvalidEndpoints() {
			while (invalidEndpoints.length) {
				const endpoint = invalidEndpoints.pop();
				endpoint.enabled = true;
				setFillStyle(endpoint, endpointColor);
			}
		}

		function getAncestors(id) {
			let ancestors = [id];
			const endpoint = getEndpoint(id, 'target');
			if (endpoint) {
				for (let i=0 ; i<endpoint.connections.length ; ++i) {
					const source = endpoint.connections[i].sourceId;
					ancestors = ancestors.concat(getAncestors(source));
				}
			}
			return ancestors;
		}

		jsPlumbInstance.bind("connection:drag", function(connection) {
			const ancestors = getAncestors(connection.endpoints[0].elementId);

			for (let i=0 ; i<ancestors.length ; ++i) {
				const endpoint = getEndpoint(ancestors[i], 'target');
				if (endpoint) {
					setEndpointInvalid(endpoint);
				}
			}		
		});

		// There's no single "connection drag finished" event - dropping a connection back where it
		// started fires nothing - so make everything visible again once the drop has been handled.
		$(document).on("mouseup", function() {
			if (invalidEndpoints.length) {
				setTimeout(resetInvalidEndpoints, 0);
			}
		});
}


function add_click_overlay(endpoint) {
	const create_overlay = function (component) {
		return $('<span><span class="attach-label">+</span></span>')[0];
	};

	// location 0 puts the label's top left on the endpoint's, so the "+" sits over the dot
	jsPlumbInstance.addOverlay(endpoint, {
		type: "Custom",
		options: {
			create: create_overlay,
			location: 0,
		}
	});
}

function add_delete_overlay(connection) {
	const create_overlay = function (component) {
		const overlay = $('<span><svg width=18 height=18><circle cx="9" cy="9" r="9" version="1.1" xmlns="http://www.w3.org/1999/xhtml" style="" stroke="none"></circle></svg><span class="detach-label cancel">-</span></span>');
		$('circle', overlay).css('fill', endpointColor);
		return overlay[0];
	};

	jsPlumbInstance.addOverlay(connection, {
		type: "Custom",
		options: {
			create: create_overlay,
			location: 0.5,
			events: {
				"click": function (params) {
					jsPlumbInstance.deleteConnection(params.overlay.component);
				}
			}
		}
	});
}


function addConnection(sourceId, targetId, side, readOnly) {
	const source = getEndpoint(sourceId, 'source');
	const target = getEndpoint(targetId, 'target', side);

	const connection = jsPlumbInstance.connect({source: source, target: target, fireEvent: false, detachable: !readOnly});
	if (!readOnly) {
		add_delete_overlay(connection);
	}
	connection.setReattach(false);
	return connection;
}


function attatchAnalysisNodeConnections(connections, readOnly) {
	for (let i=0 ; i<connections.length ; ++i) {
		const conn = connections[i];
		addConnection(conn["source_id"], conn["target_id"], conn["side"], readOnly);
	}
}


function addNewNodeToPage(data) {
	const newNode = createNodeFromData(data);
	newNode.appendTo($("#analysis"));
	setupNodes(newNode);
	return newNode;
}

function addNewNodeAndFlash(data) {
	const newNode = addNewNodeToPage(data);
	const endPoints = jsPlumbInstance.getEndpoints(newNode[0]);
	const endpointVis = function (visible) {
		$.each(endPoints, function () {
			this.setVisible(visible);
		});
	};
	endpointVis(false);

	bringNodeToFront(newNode);
	$(newNode).fadeOut(500).fadeIn(500, function() {
		endpointVis(true);
	});
}

function addNode() {
	// Get node type from select
	const nodeType = $("select#id_node_types").val();
	$.ajax({
	    type: "POST",
	    url: Urls.node_create(ANALYSIS_ID, nodeType),
	    success: addNewNodeAndFlash,
	});
}

function getActiveNodesIds() {
	const nodes = [];
	$("." + ACTIVE_CLASS).each(function() {
		const nodeId = $(this).attr('node_id');
		if (nodeId) {
		  nodes.push(nodeId);
		}
	});
	return nodes;
}

function copyNode() {
	const nodes = getActiveNodesIds();
	unselectActive(); // Unselect old, so we can select new

	const addNewNodesToPage = function (data) {
		const nodes_data = data['nodes'];
		for (let i = 0; i < nodes_data.length; ++i) {
			const node_data = nodes_data[i];
			const node = addNewNodeToPage(node_data);
			node.addClass(ACTIVE_CLASS);
		}
		const edges = data['edges'];
		attatchAnalysisNodeConnections(edges);
	};

	const data = 'nodes=' + encodeURIComponent(JSON.stringify(nodes));
	$.ajax({
	    type: "POST",
	    data: data,
	    url: Urls.nodes_copy(ANALYSIS_ID),
	    success: addNewNodesToPage,
	});
}


// Overlays (the '+' and '-' handles) are separate absolutely positioned elements to the connector canvas
function getJsPlumbElements(component) {
	const elements = [];
	// Endpoints and connections keep their element on whatever paints them
	const painted = component.endpoint || component.connector;
	if (painted && painted.canvas) {
		elements.push(painted.canvas);
	}
	const overlays = component.getOverlays ? component.getOverlays() : {};
	for (let key in overlays) {
		const overlay = overlays[key];
		if (overlay.canvas) {  // Arrows are drawn into the connector canvas and don't have one
			elements.push(overlay.canvas);
		}
	}
	return elements;
}


// Endpoints and connections are their own absolutely positioned elements, so they have to fade along with the card
function fadeOutAndRemoveNode(node) {
	const DELETE_FADE_MS = 300;
	if (node.hasClass(DELETING_CLASS)) {
		return;  // already on its way out
	}

	let fadingElements = node.get();
	$.each(jsPlumbInstance.getEndpoints(node[0]) || [], function(_, endpoint) {
		fadingElements = fadingElements.concat(getJsPlumbElements(endpoint));
		$.each(endpoint.connections || [], function(_, connection) {
			fadingElements = fadingElements.concat(getJsPlumbElements(connection));
		});
	});

	// CSS transition rather than jQuery animation - jsPlumb's teardown blocks the main thread
	$(fadingElements).addClass(DELETING_CLASS);
	setTimeout(function() {
		// Detach without events so we don't tell the server about connections going away with the node
		jsPlumbInstance.deleteConnectionsForElement(node[0], {fireEvent: false});
		jsPlumbInstance.unmanage(node[0], true);
	}, DELETE_FADE_MS);
}


function deleteNodesFromDOM(nodes, data) {
    let clearGridAndEditor = false;

    // Start every fade before anything else, so the transition isn't held up by the tear down below
    for (let i=0 ; i<nodes.length ; ++i) {
		const nodeId = nodes[i];
		const node = getNode(nodeId);

		const node_version_select = "#" + nodeId + "_" + node.attr("version_id");
		const loadedNode = $(node_version_select);
		if (node.hasClass(ACTIVE_CLASS) || loadedNode.length) {
            clearGridAndEditor = true;
        }

        fadeOutAndRemoveNode(node);
        messagePoller.delete_node(nodeId);

        // remove AnalysisVariable if exists
        if (analysisNodeVariables[nodeId]) {
        	delete analysisNodeVariables[nodeId];
        	$(".analysis-variable-node[node-id=" + nodeId + "]").remove();
		}
    }

    // Emptying a loaded grid is slow enough to be seen - let the browser paint the fade first
    setTimeout(function() {
        if (clearGridAndEditor) {
            loadNodeData(); // empty
        }
        update_dirty_nodes(data);
    }, 0);
}


function deleteNodes(nodes) {
	// Remove from the page straight away - the server call only tells us which other nodes became dirty
	deleteNodesFromDOM(nodes, []);

	const data = 'nodes=' + encodeURIComponent(JSON.stringify(nodes));
	$.ajax({
	    type: "POST",
	    data: data,
	    url: Urls.nodes_delete(ANALYSIS_ID),
	    success: update_dirty_nodes,
	});
}


function deleteActiveNodes() {
	const nodes = getActiveNodesIds();
	if (nodes.length > 0) {
		deleteNodes(nodes);
	}
}


// From: http://stackoverflow.com/a/2901298
function intWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function setVariantCount(variant_count_selector, count) {
	const countValue = $('.count-value', variant_count_selector);
	countValue.html(count);
	variant_count_selector.show();
}

// A node reading mutable tables (eg internal frequency counts) has an advisory count - a handful of
// variants cross the filter threshold with every import, so show the magnitude rather than a number
// that implies more precision than we have (#235)
function formatNodeCount(count, deterministic) {
	if (deterministic || count < 10000) {
		return intWithCommas(count);
	}
	const abbreviated = count >= 1e6 ? (count / 1e6).toFixed(1) + 'M' : Math.round(count / 1e3) + 'K';
	return '~' + abbreviated;
}

function markLiveDataCount(node_counts, deterministic, liveDataSources) {
	const totalCount = $(".node-count-" + NODE_COUNT_TOTAL, node_counts);
	totalCount.toggleClass("live-data-count", !deterministic);
	if (deterministic) {
		totalCount.removeAttr("title");
		return;
	}
	const sources = Object.keys(liveDataSources || {}).map(function(k) {
		return k + " #" + liveDataSources[k];
	}).join(", ");
	totalCount.attr("title", "Count reflects live data - " + sources + ". May differ slightly from the grid.");
	$(".count-value", totalCount).append(" <span class='live-data-marker'>&#9889;</span>");
}

function updateDirtyNode(node, refresh) {
	const node_id = node.attr("node_id");

	// Set to unknown while waiting for update
	const node_counts = $(".node-counts", node);
	$(".count-value", node_counts).empty();
	const variant_count = $(".node-count-__total", node_counts);
	setVariantCount(variant_count, '?');
    node.attr("loading", "true"); // #616 - Don't flash red when loading - this will stop next cycle of shadow setting

	const asyncUpdateNode = function (data) {
		// Flash the card border between its normal colour and shadowColor. The border is currentColor,
		// so animating the node's colour drives it - see .window.design-a-node
		const DEFAULT_COLOR = "#aaa";

		// Stopping animation ended up breaking "new node flash" so just let it time out
		//node.stop(); // any previous colours

		const nodeVersion = data["version"];
		node.attr("version_id", nodeVersion);
		node.removeAttr("loading");
		const shadowColor = data["shadow_color"];

		if (shadowColor) {
			function flashShadowColor() {
				const myNode = getNode(node_id); // get latest version
				const version_id = myNode.attr("version_id");
				const loading = myNode.attr("loading");

				if (version_id == nodeVersion && !loading) {
					myNode.animate({color: shadowColor}, 1000)
						.animate({color: DEFAULT_COLOR}, 1000, flashShadowColor);
				}
			}

			flashShadowColor();
		}

		if (data.valid) {
			const counts = data.counts;
			const deterministic = data.deterministic !== false;
			for (const c in counts) {
				const vc = $(".node-count-" + c, node_counts);
				const count = counts[c];
				if (count > 0 || vc.hasClass("show-zero")) {
					setVariantCount(vc, formatNodeCount(count, deterministic));
				} else {
					vc.hide();
				}
			}
			markLiveDataCount(node_counts, deterministic, data.live_data_sources);
		} else {
			setVariantCount(variant_count, "");
		}
		repaintNode(node);
	};
	messagePoller.observe_node(node_id, "count", asyncUpdateNode);
}

function clickCounter(evt) {
    evt.stopPropagation(); // Don't pass through click and reload node
	const node = $(this).parents(".window");
	const nodeId = node.attr("node_id");
	const countType = $(this).attr("count_type");

	unselectActive();
	setActiveNode(node);
    $(this).addClass(ACTIVE_NODE_COUNT_CLASS);
	loadNodeData(nodeId, countType, true);
} 


function attachVariantCounters(nodes_selector, nodeCountTypes) {
    drawCountLegend(nodeCountTypes);

	nodes_selector.filter("[output_endpoint=true]").each(function() {
		const strip = $(".node-counts-strip", this).empty();
		const node_counts = $("<span class='node-counts'></span>").appendTo(strip);

		for (let i=0 ; i<nodeCountTypes.length ; ++i) {
			const node_count_type = nodeCountTypes[i];
			const name = node_count_type[0];
			const data = node_count_type[1];
			const label = data["label"];
			const nodeCountContent = "<div class='user-tag-colored count-value'>?</div>";
			const nodeCountDiv = $("<div count_type=" + name + " title='" + label + "' class='node-count node-count-" + name + "'>" + nodeCountContent + "</div>");
			if (data["link"]) {
				$(nodeCountDiv).addClass("clickable-count");
				$(nodeCountDiv).click(clickCounter);
			}
			if (data["show_zero"]) {
				$(nodeCountDiv).addClass("show-zero");
			}
			window.count_data = data;
			nodeCountDiv.appendTo(node_counts);
		}
		updateDirtyNode($(this));
	});
}


function setupConnections(nodes_selector, readOnly) {
	const inputEndpoint = {
		endpoint: DOT_ENDPOINT,
		paintStyle: ENDPOINT_STYLE,
		reattachConnections: true,
		connectorStyle: CONNECTOR_STYLE,
		target: true,
	};

	const outputEndpoint = {
		endpoint: DOT_ENDPOINT,
		paintStyle: ENDPOINT_STYLE,
		source: true,
		enabled: !readOnly,
		connectorStyle: CONNECTOR_STYLE,
		connector: CONNECTOR_SPEC,
		maxConnections: -1,
	};

	const orientation = getNodeOrientation();

	nodes_selector.each(function() {
		// Nodes without endpoints still need to be managed, so they can be dragged
		jsPlumbInstance.manage(this);

		const uuid = "input-endpoint-" + $(this).attr("node_id");
		const existingEndpoint = jsPlumbInstance.getEndpoint(uuid);

		if ($(this).attr("input_endpoint") == "true") {
			if ($(this).hasClass("VennNode")) {
				const inputLeft = uuid + "-left";
				const inputRight = uuid + "-right";
				jsPlumbInstance.addEndpoint(this, { anchor: orientation.vennInputs[0], uuid: inputLeft}, inputEndpoint);
				jsPlumbInstance.addEndpoint(this, { anchor: orientation.vennInputs[1], uuid: inputRight}, inputEndpoint);
            } else if ($(this).hasClass("MergeNode")) {
				const commonInputEndpoint = $.extend({}, inputEndpoint);
				commonInputEndpoint["endpoint"] = {type: "Dot", options: {radius: 18}};
                commonInputEndpoint["maxConnections"] = -1;
                jsPlumbInstance.addEndpoint(this, { anchor: orientation.input, uuid: uuid}, commonInputEndpoint);
			} else {
			    // Keep existing, so we don't
			    if (!existingEndpoint) {
				    jsPlumbInstance.addEndpoint(this, { anchor: orientation.input, uuid: uuid}, inputEndpoint);
				}
			}
		} else {
            if (existingEndpoint) {
                jsPlumbInstance.deleteEndpoint(uuid);
            }
		}

		if ($(this).attr("output_endpoint") == "true") {
			const e = jsPlumbInstance.addEndpoint(this, {anchor: orientation.output}, outputEndpoint);
			if (!readOnly) {
				add_click_overlay(e);
			}
		}
	});
}

function loadNodeWhenReady(node_id) {
	const updateNode = function (node) {
		// console.log("loadNodeWhenReady.updateNode:" + node);
		// Reload the data container if it's showing for this node
		const gew = getGridAndEditorWindow();
		if (node_id == gew.getLoadedNodeId()) {
			loadNodeData(node_id);
		}
	};
	messagePoller.observe_node(node_id, "ready", updateNode);
}


function setActiveNode(node) {
	node.addClass(ACTIVE_CLASS);
	bringNodeToFront(node);
}


function unselectActive() {
	const container = $("#analysis-container");
	$("." + ACTIVE_CLASS, container).removeClass(ACTIVE_CLASS);
    $("." + ACTIVE_NODE_COUNT_CLASS, container).removeClass(ACTIVE_NODE_COUNT_CLASS);
}


function setupNodes(nodes_selector, readOnly) {
	nodes_selector.click(function() {
		unselectActive();
		setActiveNode($(this));
		const nodeId = $(this).attr('node_id');
		loadNodeData(nodeId, null, true);
	});
	
	setupConnections(nodes_selector, readOnly);
	
	const nodeCountTypes = ANALYSIS_SETTINGS['node_count_types'];
	if (nodeCountTypes) {
		attachVariantCounters(nodes_selector, nodeCountTypes);
	}
	
	if (!readOnly) {
		setupNodeModifications();
	}
}

function getMaxZIndex(nodes_selector) {
	let maxZ = null;
	nodes_selector.each(function() {
		const zIndex = parseFloat($(this).css('z-index'));
		maxZ = (zIndex > maxZ) ? zIndex : maxZ;
	});
	return maxZ;
}

function bringNodeToFront(node) {
		const zIndex = getMaxZIndex($(".window")) + 1;
		node.css("z-index", zIndex);
}


function setupNodeModifications() {
	// The marquee only starts on empty canvas - grabbing a node hands the mouse to jsPlumb's dragging.
	// unselectActive on start to also unselect node counts
	$('#analysis-container').selectable({
		filter: "div.window",
		start: function () {
			// clear editor/grid so people don't get confused about active node
			replaceEditorWindow();
			loadGridAndEditorForNode();
		},
		cancel: '.cancel, div.window'
	});
}


// jsPlumb moves everything in its drag selection, so mirror the currently selected nodes into it as
// the drag begins - grabbing a node that isn't selected moves only that node.
function syncDragSelection(params) {
	jsPlumbInstance.clearDragSelection();
	const dragElement = params.drag.getDragElement();
	if (dragElement.classList.contains(ACTIVE_CLASS)) {
		jsPlumbInstance.addToDragSelection.apply(jsPlumbInstance, $(".window." + ACTIVE_CLASS, "#analysis").toArray());
	}
}


function getEndpointSide(ep) {
	const uuid = ep.getUuid();
	return uuid.split("-").reverse()[0];
}


const updateConnections = function (info, remove) {
	const sourceId = $(info.source).attr("node_id");
	const target = $(info.target);
	const targetId = target.attr("node_id");
	const params = {parent_id: sourceId, remove: remove};

	if (!remove && target.hasClass("VennNode")) {
		const ep = info.connection.endpoints[1];
		params["side"] = getEndpointSide(ep);
	}

	const on_success_function = function () {
		const gew = getGridAndEditorWindow();
		if (targetId == gew.getLoadedNodeId()) {
			loadNodeData(targetId);
		}
		checkAndMarkDirtyNodes();
	};

	updateNode(targetId, 'update_connection', params, on_success_function);
};


function saveDraggedNodePositions(params) {
	for (let i=0 ; i<params.elements.length ; ++i) {
		const dragged = params.elements[i];
		const nodeId = dragged.el.getAttribute("node_id");
		updateNode(nodeId, 'move', {'x': dragged.pos.x, 'y': dragged.pos.y});
	}

	// Avoid click event at end of drag
	// http://stackoverflow.com/questions/3486760/how-to-avoid-jquery-ui-draggable-from-also-triggering-click-event/13973319#13973319
	if (params.e) {
		$(params.e.target).one('click', function(e) { e.stopImmediatePropagation(); });
	}
}


window.variantgridPipeline = {init : function(readOnly) {
	jsPlumbInstance = jsPlumb.newInstance({
		container: document.getElementById("analysis"),
		elementsDraggable: !readOnly,
		endpoint: DOT_ENDPOINT,
		endpointStyle: ENDPOINT_STYLE,
		connector: CONNECTOR_SPEC,
		paintStyle: CONNECTOR_STYLE,
		anchors: getNodeOrientation().instanceAnchors,
		dragOptions: {cursor: 'pointer', zIndex: 2000, beforeStart: syncDragSelection},
	});

	// bind to connection/detach events, and update the list of connections on screen.
	jsPlumbInstance.bind("connection", function(info, originalEvent) {
		updateConnections(info);
		add_delete_overlay(info.connection);
	});
	jsPlumbInstance.bind("connection:detach", function(info, originalEvent) {
		updateConnections(info, true);
	});

	// The dragged node - and anything dragged along with it - sits above the rest of the graph
	jsPlumbInstance.bind("drag:start", function(params) {
		bringNodeToFront($(params.el));
		$(".window." + ACTIVE_CLASS, "#analysis").each(function() {
			bringNodeToFront($(this));
		});
	});
	jsPlumbInstance.bind("drag:stop", saveDraggedNodePositions);

	setupHideInvalidConnectionsOnDrag();

	const nodes_selector = $(".window");
	setupNodes(nodes_selector, readOnly);
}};


function drawCountLegend(nodeCountTypes) {
	const legend = $("#node-count-legend");
	legend.empty();

    if (nodeCountTypes && nodeCountTypes.length) {
        $("<div><b>Node Counts:</b></div>").appendTo(legend);
        for (let i=0; i<nodeCountTypes.length ; ++i) {
			const node_count_type = nodeCountTypes[i];
			const name = node_count_type[0];
			const data = node_count_type[1];
			const row = $("<div class='legend-row node-count-legend-" + name + "'><div class='user-tag-colored legend-color-box'></div><div class='legend-label'>" + data.label + "</div><div class='clear'></div></div>");
			legend.append(row);
        }
    }

    if (!isHorizontalNodeFlow()) {
        // Absolute positioning can sometimes get thrown off after resizing contents
        const rowHeight = 31;
        const numRows = 1 + nodeCountTypes.length;
        legend.height(rowHeight * numRows);
    }
}


function get_array_second_dimension_elements(array, element_id) {
	const element_values = [];
	for (let i=0 ; i<array.length ; ++i) {
        element_values.push(array[i][element_id]);
    }
    return element_values; 
}

function get_array_element_keys(array, key_id) {
	const element_values = [];
	for (let i=0 ; i<array.length ; ++i) {
		const value = array[i][key_id];
		element_values.push(value);
    }
    return element_values; 
}


function changeAnalysisSettings(oldAnalysisSettings) {
	let requireReload = false;
	variantTagStaleDays = ANALYSIS_SETTINGS.variant_tag_stale_days;
	const oldAnnotationVersion = oldAnalysisSettings.annotation_version;
	const newAnnotationVersion = ANALYSIS_SETTINGS.annotation_version;
	requireReload = oldAnnotationVersion != newAnnotationVersion;

	const oldNodeCountTypes = oldAnalysisSettings.node_count_types;
	const newNodeCountTypes = ANALYSIS_SETTINGS.node_count_types;
	const arraysAreDifferent = true; // TODO: Add a way to test this? Could save a flash of attaching/reattaching counters
	if (arraysAreDifferent) {
        // If we just removed some it's ok to just change the counters
        // If there are new ones, need to force reload 
		const oldCountData = get_array_second_dimension_elements(oldNodeCountTypes, 1);
		const newCountData = get_array_second_dimension_elements(newNodeCountTypes, 1);
		const oldTypes = get_array_element_keys(oldCountData, "label");
		const newTypes = get_array_element_keys(newCountData, "label");

		for (let i=0 ; i<newTypes.length ; ++i) {
            const count_type = newTypes[i];
            if ($.inArray(count_type, oldTypes) == -1) {
                requireReload = true;
                break;
            }
        }

		const nodes_selector = $(".window");
		attachVariantCounters(nodes_selector, newNodeCountTypes);
    }
    
    if (requireReload) {
        reloadNodes();
    }        
}


function update_dirty_nodes(dirty_nodes) {
	const nodeCountTypes = ANALYSIS_SETTINGS.node_count_types;
	if (nodeCountTypes) {
		for (let i=0 ; i<dirty_nodes.length ; ++i) {
			const node_id = dirty_nodes[i];
			const node = getNode(node_id);
			if (node.length) {
				updateDirtyNode(node, true);
			}
		}
	}
}

function updateNodeAppearance(data) {
    const nodeId = data.attributes.node_id;
    const node = getNode(nodeId);
    updateNodeFromData(node, data);
    setupConnections(node);
}


function retrieveAndUpdateNodeAppearances(nodeList) {
    for (let i=0 ; i<nodeList.length ; ++i) {
        const nodeId = nodeList[i];
        $.ajax({
            url: Urls.node_data(ANALYSIS_ID, nodeId),
            success: updateNodeAppearance,
        });
    }
}


function addConnectedNode(data) {
	const node = addNewNodeToPage(data);
	const sourceId = data.node_id;
	const targetId = node.attr('id');
	addConnection(sourceId, targetId);
}


function loggedOutHandler() {
    showReloadPageErrorDialog($("#error-dialog"), "You have been logged out.", true);
}


function checkAndMarkDirtyNodes(aWin) {
    if (typeof aWin == 'undefined') {
        aWin = getAnalysisWindow();        
    }

    $.ajax({
        url: Urls.analysis_node_versions(ANALYSIS_ID),
        success: function(data) {
			const dirty_nodes = [];
			const appearance_update_nodes = [];

			for (let i=0 ; i<data.node_versions.length ; ++i) {
                const nodeData = data.node_versions[i];
                const nodeId = nodeData[0];
                const nodeVersion = nodeData[1];
                const nodeAppearanceVersion = nodeData[2];
                const node = getNode(nodeId);
                const localVersion = node.attr("version_id");
                if (localVersion != nodeVersion) {
                    dirty_nodes.push(nodeId);
                }
                const localAppearanceVersion = node.attr("appearance_version_id");
                if (localAppearanceVersion != nodeAppearanceVersion) {
                    appearance_update_nodes.push(nodeId);
                }
            }
            if (dirty_nodes.length) {
                aWin.update_dirty_nodes(dirty_nodes);
            }
            if (appearance_update_nodes.length) {
                aWin.retrieveAndUpdateNodeAppearances(appearance_update_nodes);
            }
        },
    });
}


function reloadNodeAndData(node_id) {
	const aWin = getAnalysisWindow();
	checkAndMarkDirtyNodes(aWin);
    aWin.loadNodeData(node_id);
    // An editor save is the one time we want to leave the editor - show the results being re-computed
    if (typeof showBottomPaneTab === "function") {
        showBottomPaneTab(BOTTOM_PANE_GRID);
    }
}
