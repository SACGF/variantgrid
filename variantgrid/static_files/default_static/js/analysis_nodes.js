var endpointColor = "#121212";
var ACTIVE_CLASS = "ui-selected";
var ACTIVE_NODE_COUNT_CLASS = "node-counts-selected"; // Needs to be different so not grabbed with multi-draggable
var SHOW_NODE_IDS_IN_TOOLTIPS = true;

function getNode(nodeId) {
	return $("#analysis-node-" + nodeId);
}

// Nodes should also create a "updateState(args)" method, this will be called after creation  

function createDefaultNode(nodeData) {
	var div = jQuery('<div></div>', nodeData["attributes"]);
	div.addClass("default-node-container");
	div.append("<div class='user-tag-colored node-overlay'><span class='node-name'>" + nodeData["name"] + "</span></div>");
	div[0].updateState = function(args) { };
	return div;
}

function createVennNode(nodeData) {
	var div = jQuery('<div></div>', nodeData["attributes"]);
	div.addClass("default-node-container");
	venn2(div[0], 64, 45);

	overlay_style = {	width: '100%',
						height: '100%',
		    			position: 'absolute',
	    				top: 0,
	 	 	 			left: 0,
	};
	$('.' + VENN_TOGGLE_WIDGET_CLASS, div[0]).css(overlay_style);
	span = $("<span class='node-name'>" + nodeData["name"] + "</span>");
	span.css(overlay_style);
	span.css("z-index", 30);
	div.append(span);
	div[0].updateState = function(args) {
		venn_select(this, args['venn_flag']);
	};
	return div;
}

function createNodeFromData(nodeData) {
	const NODE_FACTORIES = {
		"SampleNode" : createSampleNode,
		"VennNode" : createVennNode
	};

	let nodeClass = nodeData['node_class'];
	let factory = NODE_FACTORIES[nodeClass];
	if (!factory) {
		factory = createDefaultNode;
	}
	let node = factory(nodeData);
	let nodeColorOverlay = $("<div />").attr({class: "node-color-overlay"});
	nodeColorOverlay.appendTo(node);

	node.attr("node_name", nodeData.name);
	node.each(function() { this.updateState(nodeData['args']); });
	return node;
}

// Wait for previous update to come back (so we don't end up with race conditions on the server)
var waiting_for_message_callback = false;
var update_message_buffer = [];

function updateNode(nodeId, op, params, on_success_function) {
	var message = [nodeId, op, params, on_success_function];
	update_message_buffer.push(message);
	checkSendMessage();
}


function checkSendMessage() {
	if (!waiting_for_message_callback) {
		var message = update_message_buffer.shift();
		if (message) {
			var nodeId = message[0];
			var op = message[1];
			var params = message[2];
			var old_on_success_function = message[3]; 
			var on_success_function = function(data) {
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
	var data = '&op=' + op + '&params=' + JSON.stringify(params);
	$.ajax({
	    type: "POST",
	    data: data,
	    url: Urls.node_update(nodeId),
	    success: function(data) {
	       on_success_function(data.dirty_nodes);
        },
	});
}

function addNodesToDOM(selector, nodeDataArray, readOnly) {
	selector = $(selector);
	for (var i=0 ; i<nodeDataArray.length ; ++i) {
		node = createNodeFromData(nodeDataArray[i]);
		selector.append(node);
	}
}


function getEndpoint(id, endpoint_type, side) {
	var selector = '#' + id;
	var endpoints = jsPlumb.getEndpoints($(selector));
    if (!endpoints)
        return null;    

	for (var i=0 ; i<endpoints.length ; ++i) {
		ep = endpoints[i];
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
		// TODO: could possibly use endpointDropForbiddenClass etc to do this rather than switching fill styles
		function setFillStyle(endpoint, color) {
			var paintStyle = endpoint.getPaintStyle();
			paintStyle.fillStyle = color;
			endpoint.setPaintStyle(paintStyle);
		}
		function setEndpointInvalid(endpoint) {
		    endpoint.setEnabled(false);
		    setFillStyle(endpoint, "#FF0000");
		}

		// Returns true if endpoint was "dirty"
		function resetEndpoint(endpoint) {
			var dirty = endpoint.isEnabled();
			endpoint.setEnabled(true);
			setFillStyle(endpoint, endpointColor);
			return dirty;
		}

		function getAncestors(id) {
			var ancestors = [id];
			var endpoint = getEndpoint(id, 'target');
			if (endpoint) {
				for (var i=0 ; i<endpoint.connections.length ; ++i) {
					source = endpoint.connections[i].sourceId;
					ancestors = ancestors.concat(getAncestors(source));
				}
			}
			return ancestors;
		}

		jsPlumb.bind("connectionDrag", function(c) {
		    var id = c.endpoints[0].elementId;
			var ancestors = getAncestors(id);

			for (var i=0 ; i<ancestors.length ; ++i) {
				id = ancestors[i];
				var endpoint = getEndpoint(id, 'target');
				if (endpoint) {
					setEndpointInvalid(endpoint);
				}
			}		
		});

		// Make everything visible again
		jsPlumb.bind("connectionDragStop", function(c) {
			$(".window").each(function() {
				var endpoints = jsPlumb.getEndpoints($(this));
				var dirty = false;
				for (var i=0 ; i<endpoints.length ; ++i) {
					var endpoint = endpoints[i];
					dirty |= resetEndpoint(endpoint);
				}
				if (dirty) {
					jsPlumb.repaint(this);
				}
			});
		});

}


function add_click_overlay(endpoint) {
	var create_overlay = function(component) {
		return $('<span><span class="attach-label">+</span></span>');
	};

	var overlay = ["Custom", {
		create: create_overlay,
		location: [ 0, 0 ],
	}];
	endpoint.addOverlay(overlay);
}

function add_delete_overlay(connection) {
	var create_overlay = function(component) {
		var overlay = $('<span><svg width=18 height=18><circle cx="9" cy="9" r="9" version="1.1" xmlns="http://www.w3.org/1999/xhtml" style="" stroke="none"></circle></svg><span class="detach-label cancel">-</span></span>');
		$('circle', overlay).css('fill', endpointColor);
		return overlay;
	};

	var overlay = ["Custom", {
		create: create_overlay,	    			        				 
		label: "-",
		location: 0.5,
		events: {
			"click": function(overlay, evt) {
				jsPlumb.detach(overlay.component);
			}
	}}];
	connection.addOverlay(overlay);
}


function addConnection(sourceId, targetId, side, readOnly) {
	var source = getEndpoint(sourceId, 'source');
	var target = getEndpoint(targetId, 'target', side);

	var connection = jsPlumb.connect({source: source, target: target, fireEvent: false, detachable: !readOnly});
	if (!readOnly) {
		add_delete_overlay(connection);
	}
	connection.setReattach(false);
	return connection;
}


function attatchAnalysisNodeConnections(connections, readOnly) {
	for (var i=0 ; i<connections.length ; ++i) {
		var conn = connections[i];
		addConnection(conn["source_id"], conn["target_id"], conn["side"], readOnly);
	}
}


function addNewNodeToPage(data) {
	newNode = createNodeFromData(data);
	newNode.appendTo($("#analysis"));
	setupNodes(newNode);
	return newNode;
}

function addNewNodeAndFlash(data) {
	var newNode = addNewNodeToPage(data);
	var endPoints = jsPlumb.getEndpoints(newNode);
	var endpointVis = function(visible) {
		$.each(endPoints, function() { this.setVisible(visible); });
	};
	endpointVis(false);
	
	$(newNode).fadeOut(500).fadeIn(500, function() {
		endpointVis(true);
	});
}

function addNode() {
	// Get node type from select
	var nodeType = $("select#id_node_types").val();

	$.ajax({
	    type: "POST",
	    url: Urls.node_create(ANALYSIS_ID, nodeType),
	    success: addNewNodeAndFlash,
	});
}

function getActiveNodesIds() {
	var nodes = [];
	$("." + ACTIVE_CLASS).each(function() {
		var nodeId = $(this).attr('node_id');
		if (nodeId) {
		  nodes.push(nodeId);
		}
	});
	return nodes;
}

function copyNode() {
	let nodes = getActiveNodesIds();
	unselectActive(); // Unselect old, so we can select new

	var addNewNodesToPage = function(data) {
		let nodes_data = data['nodes'];
		for (let i=0 ; i<nodes_data.length ; ++i) {
			let node_data = nodes_data[i];
			let node = addNewNodeToPage(node_data);
			node.addClass(ACTIVE_CLASS);
		}
		let edges = data['edges'];
		attatchAnalysisNodeConnections(edges);
	};

	let nodeId = $(this).attr('node_id');
	let data = 'nodes=' + encodeURIComponent(JSON.stringify(nodes));
	$.ajax({
	    type: "POST",
	    data: data,
	    url: Urls.nodes_copy(ANALYSIS_ID),
	    success: addNewNodesToPage,
	});
}


function deleteNodesFromDOM(nodes, data) {
    for (var i=0 ; i<nodes.length ; ++i) {
        var nodeId = nodes[i];
        var node = getNode(nodeId);
        
        // Delete grid and editor if it's open
        var node_version_select = "#" + nodeId + "_" + node.attr("version_id");
        var loadedNode = $(node_version_select);
        if (node.hasClass(ACTIVE_CLASS) || loadedNode.length) {
            loadNodeData(); // empty
        }

        // Detatch connections first (without triggering events) so we don't send anything to server upon deletion
        jsPlumb.detachAllConnections(node, {fireEvent: false});
        jsPlumb.remove(node);

        messagePoller.delete_node(nodeId);

        // remove AnalysisVariable if exists
        if (analysisNodeVariables[nodeId]) {
        	delete analysisNodeVariables[nodeId];
        	$(".analysis-variable-node[node-id=" + nodeId + "]").remove();
		}

    }
    update_dirty_nodes(data);
}


function deleteNodes(nodes) {
	var data = 'nodes=' + encodeURIComponent(JSON.stringify(nodes));
	$.ajax({
	    type: "POST",
	    data: data,
	    url: Urls.nodes_delete(ANALYSIS_ID),
	    success: function(data) {
	        deleteNodesFromDOM(nodes, data);
	    },
	});
}


function deleteActiveNodes() {
	var nodes = getActiveNodesIds();
	if (nodes.length > 0) {
		deleteNodes(nodes);
	}
}


// From: http://stackoverflow.com/a/2901298
function intWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function setVariantCount(variant_count_selector, count) {
	var countValue = $('.count-value', variant_count_selector);
	countValue.html(count);
	variant_count_selector.show();
}

function updateDirtyNode(node, refresh) {
	var node_id = node.attr("node_id");

	// Set to unknown while waiting for update
	var node_counts = $(".node-counts", node);
	$(".count-value", node_counts).empty();
	var variant_count = $(".node-count-__total", node_counts);
	setVariantCount(variant_count, '?');
    node.attr("loading", "true"); // #616 - Don't flash red when loading - this will stop next cycle of shadow setting

	var asyncUpdateNode = function(data) {
		var DEBUG = 0;
		// Flash between normal shadow and shadowColor
		var DEFAULT_COLOR = "#aaa";

		// Stopping animation ended up breaking "new node flash" so just let it time out
		//node.stop(); // any previous colours
		var shadow = $(".sample-shadow", node);
		shadow.css({"fill" : DEFAULT_COLOR});

		var nodeVersion = data["version"];
		node.attr("version_id", nodeVersion);
		node.removeAttr("loading");
		var shadowColor = data["shadow_color"];

		if (shadowColor) {
			if (node.hasClass("default-node-container")) {
				function defaultRunIt() {
					if (DEBUG) { console.log("defaultRunIt()"); }
					var myNode = getNode(node_id); // get latest version
					var version_id = myNode.attr("version_id");
					var loading = myNode.attr("loading");
					if (version_id == nodeVersion && !loading) {
						if (DEBUG) { console.log(".animate()"); }
				        myNode.animate({ color: shadowColor }, 1000)
				        	  .animate({ color: DEFAULT_COLOR }, 1000, defaultRunIt);
					}
					if (DEBUG) { console.log("end defaultRunIt()"); }

			    }
			    defaultRunIt();
		   	} else if (node.hasClass("SampleNode")) {
				function sampleRunIt(i) {
					if (DEBUG) { console.log("sampleRunIt()"); }
					var SHADOW_COLORS = [shadowColor, DEFAULT_COLOR];
					var myNode = getNode(node_id); // get latest version
					var version_id = myNode.attr("version_id");
                    var loading = myNode.attr("loading");
					if (version_id == nodeVersion && !loading) {
						setTimeout(function() {
							shadow.css({"fill" : SHADOW_COLORS[i%2],
										"transition" : '1.0s'});
							sampleRunIt(i+1);
					    }, 1000);
					}
			    }
			    sampleRunIt(0);
		   	}
		}
	
		if (DEBUG) { console.log("done setting shadowColor"); }

		if (data.valid) {
			var counts = data.counts;
			for (var c in counts) {
				var vc = $(".node-count-" + c, node_counts);
				var count = counts[c];
				if (count > 0 || vc.hasClass("show-zero")) {
					setVariantCount(vc, intWithCommas(counts[c]));
				} else {
					vc.hide();				
				}
			}
		} else {
			console.log("Error: Node " + node_id);
			console.log(data);
			setVariantCount(variant_count, "");
		}
	};
	messagePoller.observe_node(node_id, "count", asyncUpdateNode);
}

function clickCounter(evt) {
    evt.stopPropagation(); // Don't pass through click and reload node
	var node = $(this).parents(".window");
	var nodeId = node.attr("node_id");
	var countType = $(this).attr("count_type");

	unselectActive();
    node.addClass(ACTIVE_CLASS);
    $(this).addClass(ACTIVE_NODE_COUNT_CLASS);
	loadNodeData(nodeId, countType, true);
} 


function attachVariantCounters(nodes_selector, nodeCountTypes) {
    drawCountLegend(nodeCountTypes);
	const COUNTER_SIZE = 20;
	const LOCK_SIZE = 16;
	
	nodes_selector.filter(".outputEndpoint").each(function() {
		var counts_size = COUNTER_SIZE * nodeCountTypes.length;

		var node_width = $(this).width();
		var count_overlay = $("<div class='count-overlay'><span class='node-counts'></span></div>");
		count_overlay.css({width: node_width + COUNTER_SIZE, height: $(this).height() + counts_size});
		count_overlay.appendTo(this);
		var node_counts = $(".node-counts", count_overlay);
		node_counts.css({height: counts_size});

		for (var i=0 ; i<nodeCountTypes.length ;++i) {
			var node_count_type = nodeCountTypes[i];
			var name = node_count_type[0];
			var data = node_count_type[1];
			var label = data["label"];
			var nodeCountContent = "<div class='user-tag-colored count-value'>?</div>";
			var nodeCountDiv = $("<div count_type=" + name + " title='" + label + "' class='node-count node-count-" + name + "'>" + nodeCountContent + "</div>");
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
	// configure some drop options for use by all endpoints.
	let exampleDropOptions = {
		tolerance:"touch",
		hoverClass:"dropHover",
		activeClass:"dragActive"
	};

	let connectorStyle = {
		lineWidth : 3,
		strokeStyle: endpointColor,
	};
	
	let inputEndpoint = {
		endpoint: ["Dot", { radius:11 }],
		paintStyle:{ width:25, height:21, fillStyle: endpointColor },
		reattach:true,
		connectorStyle : connectorStyle,
		//maxConnections: -1,
		isTarget:true,
		dropOptions : exampleDropOptions
	};			
	
	var outputEndpoint = {
		endpoint:["Dot", { radius:11 }],
		paintStyle:{ fillStyle: endpointColor },
		isSource: true,
		enabled: !readOnly,
		connectorStyle : connectorStyle,
		connector: ["Bezier", { curviness:63 } ],
		maxConnections: -1,
		dropOptions : exampleDropOptions
	};

	nodes_selector.each(function() {
        var uuid = "input-endpoint-" + $(this).attr("node_id");
        var existingEndpoint = jsPlumb.getEndpoint(uuid);

		if ($(this).hasClass("inputEndpoint")) {
			if ($(this).hasClass("VennNode")) {
				var topLeft = uuid + "-left";
				var topRight = uuid + "-right";
				jsPlumb.addEndpoint(this, { anchor: [0.25, 0, 0, -1], uuid: topLeft}, inputEndpoint);
				jsPlumb.addEndpoint(this, { anchor: [0.75, 0, 0, -1], uuid: topRight}, inputEndpoint);
            } else if ($(this).hasClass("MergeNode")) {
                var commonInputEndpoint = $.extend({}, inputEndpoint);
                commonInputEndpoint["endpoint"] = ["Dot", { radius: 18 }];
                commonInputEndpoint["maxConnections"] = -1;
                jsPlumb.addEndpoint(this, { anchor:"TopCenter", uuid: uuid}, commonInputEndpoint);
			} else {
			    // Keep existing, so we don't 
			    if (!existingEndpoint) {
				    jsPlumb.addEndpoint(this, { anchor:"TopCenter", uuid: uuid}, inputEndpoint);
				}
			}
		} else {
            if (existingEndpoint) {
                jsPlumb.deleteEndpoint(uuid);
            }
		}
		
		if ($(this).hasClass("outputEndpoint")) {
			var e = jsPlumb.addEndpoint(this, { anchor:"BottomCenter" }, outputEndpoint);
			if (!readOnly) {
				add_click_overlay(e);
			}
		}
	});
}

function loadNodeWhenReady(node_id) {
	var updateNode = function(node) {
		// console.log("loadNodeWhenReady.updateNode:" + node);
		// Reload the data container if it's showing for this node
		var gew = getGridAndEditorWindow();
		if (node_id == gew.getLoadedNodeId()) {
			loadNodeData(node_id);
		}
	};
	messagePoller.observe_node(node_id, "ready", updateNode);
}


function unselectActive() {
    var container = $("#analysis-container");
	$("." + ACTIVE_CLASS, container).removeClass(ACTIVE_CLASS);
    $("." + ACTIVE_NODE_COUNT_CLASS, container).removeClass(ACTIVE_NODE_COUNT_CLASS);
}


function setupNodes(nodes_selector, readOnly) {
	nodes_selector.click(function() {
		unselectActive();
		$(this).addClass(ACTIVE_CLASS);
		let nodeId = $(this).attr('node_id');
		loadNodeData(nodeId, null, true);
	});
	
	setupConnections(nodes_selector, readOnly);
	
	let nodeCountTypes = ANALYSIS_SETTINGS['node_count_types'];
	if (nodeCountTypes) {
		attachVariantCounters(nodes_selector, nodeCountTypes);
	}
	
	if (Object.keys(NODE_HELP).length) { // Empty if no tooltips
		nodes_selector.each(function () {
			// test global if we should assign
			let nodeClass = $(this).attr("node_class");
			let node_help = NODE_HELP[nodeClass];

			if (SHOW_NODE_IDS_IN_TOOLTIPS) {
				let nodeId = $(this).attr('node_id');
				node_help += " (node #" + nodeId + ")";
			}
			$(this).attr('title', node_help);
		});
	}

	if (!readOnly) {
		setupNodeModifications(nodes_selector);
	}
}


function setupNodeModifications(nodes_selector) {
	function dragStop(event) {
		let nodeId = $(this).attr("node_id");
		let position = $(this).position();
		let params = {'x' : position.left, 'y' : position.top};
		updateNode(nodeId, 'move', params);

		// Avoid click event at end of drag
		// http://stackoverflow.com/questions/3486760/how-to-avoid-jquery-ui-draggable-from-also-triggering-click-event/13973319#13973319
		$( event.toElement ).one('click', function(e){ e.stopImmediatePropagation(); } );
	}

	let params = {
		activeClass: ACTIVE_CLASS,
		stopNative: dragStop,
		stopAll: dragStop
	};
	nodes_selector.multiDraggable(params);

	// unselectActive on start to also unselect node counts
	$('#analysis-container').selectable({
		filter: "div.window",
		start: function () {
			// clear editor/grid so people don't get confused about active node
			replaceEditorWindow();
			loadGridAndEditorForNode();
		},
		cancel: '.cancel'
	});
}


function getEndpointSide(ep) {
	var uuid = ep.getUuid();
	var side = uuid.split("-").reverse()[0];
	return side;
}


jsPlumb.ready(function() {
	var _initialised = false;

	var updateConnections = function(info, remove) {
        var sourceId = $("#" + info.sourceId).attr("node_id");
        var target = $("#" + info.targetId);
		var targetId = target.attr("node_id");
		var params = {parent_id: sourceId, remove: remove};

		if (!remove && target.hasClass("VennNode")) {
			var ep = info.connection.endpoints[1];
			params["side"] = getEndpointSide(ep); 
		} 

		var on_success_function = function() {
            var gew = getGridAndEditorWindow();
            if (targetId == gew.getLoadedNodeId()) {
                loadNodeData(targetId);
            }
            checkAndMarkDirtyNodes();
		};

		updateNode(targetId, 'update_connection', params, on_success_function);
	};

    window.variantgridPipeline = {init : function(readOnly) {
		// setup jsPlumb defaults.
		jsPlumb.importDefaults({
			DragOptions : { cursor: 'pointer', zIndex:2000 },
			PaintStyle : { strokeStyle:'#666' },
			EndpointStyle : { width:20, height:16, strokeStyle:'#666' },
			Endpoint : "Rectangle",
			Anchors : ["TopCenter", "TopCenter"]
		});												

		// bind to connection/connectionDetached events, and update the list of connections on screen.
		jsPlumb.bind("connection", function(info, originalEvent) {
			updateConnections(info);
			add_delete_overlay(info.connection);
		});
		jsPlumb.bind("connectionDetached", function(info, originalEvent) {
			updateConnections(info, true);
		});

		setupHideInvalidConnectionsOnDrag();		

		var nodes_selector = $(".window");
		setupNodes(nodes_selector, readOnly);
		
		if (!_initialised) {
			$(".drag").click(function() {
				var s = jsPlumb.toggleDraggable($(this).attr("rel"));
				$(this).html(s ? 'disable dragging' : 'enable dragging');
				if (!s)
					$("#" + $(this).attr("rel")).addClass('drag-locked');
				else
					$("#" + $(this).attr("rel")).removeClass('drag-locked');
				$("#" + $(this).attr("rel")).css("cursor", s ? "pointer" : "default");
			});

			$(".detach").click(function() {
				jsPlumb.detachAllConnections($(this).attr("rel"));
			});

			_initialised = true;
		}
	}
};
});


function drawCountLegend(nodeCountTypes) {
    var legend = $("#node-count-legend");
    legend.empty();

    if (nodeCountTypes && nodeCountTypes.length) {
        $("<div><b>Node Counts:</b></div>").appendTo(legend);
        for (var i=0; i<nodeCountTypes.length ; ++i) {
            var node_count_type = nodeCountTypes[i];
            var name = node_count_type[0];
            var data = node_count_type[1];
            var row = $("<div class='legend-row node-count-legend-" + name + "'><div class='user-tag-colored legend-color-box'></div><div class='legend-label'>" + data.label + "</div><div class='clear'></div></div>");
            legend.append(row);
        }
    }
}


function get_array_second_dimension_elements(array, element_id) {
    var element_values = [];
    for (var i=0 ; i<array.length ; ++i) {
        element_values.push(array[i][element_id]);
    }
    return element_values; 
}

function get_array_element_keys(array, key_id) {
    var element_values = [];
    for (var i=0 ; i<array.length ; ++i) {
        var value = array[i][key_id];
        element_values.push(value);
    }
    return element_values; 
}


function changeAnalysisSettings(oldAnalysisSettings) {
    var requireReload = false;
    var oldAnnotationVersion = oldAnalysisSettings.annotation_version;
    var newAnnotationVersion = ANALYSIS_SETTINGS.annotation_version;
    requireReload = oldAnnotationVersion != newAnnotationVersion;

	var oldNodeCountTypes = oldAnalysisSettings.node_count_types;
	var newNodeCountTypes = ANALYSIS_SETTINGS.node_count_types;
	var arraysAreDifferent = true; // TODO: Add a way to test this? Could save a flash of attaching/reattaching counters
	if (arraysAreDifferent) {
        // If we just removed some it's ok to just change the counters
        // If there are new ones, need to force reload 
        var oldCountData = get_array_second_dimension_elements(oldNodeCountTypes, 1);
        var newCountData = get_array_second_dimension_elements(newNodeCountTypes, 1);
        var oldTypes = get_array_element_keys(oldCountData, "label");
        var newTypes = get_array_element_keys(newCountData, "label");
        
        for (var i=0 ; i<newTypes.length ; ++i) {
            count_type = newTypes[i];
            if ($.inArray(count_type, oldTypes) == -1) {
                requireReload = true;
                break;
            }
        }

        $(".count-overlay").remove();
        var nodes_selector = $(".window");
        attachVariantCounters(nodes_selector, newNodeCountTypes);
    }
    
    if (requireReload) {
        reloadNodes();
    }        
}


function update_dirty_nodes(dirty_nodes) {
	var nodeCountTypes = ANALYSIS_SETTINGS.node_count_types;
	if (nodeCountTypes) {
		for (i=0 ; i<dirty_nodes.length ; ++i) {
			var node_id = dirty_nodes[i];
			var node = getNode(node_id);
			if (node.length) {
				updateDirtyNode(node, true);
			}
		}
	}
}

function updateNodeAppearance(data) {
    let nodeId = data.attributes.node_id;
    let oldNode = getNode(nodeId);
    let wasActive = oldNode.hasClass(ACTIVE_CLASS);
    oldNode.remove();

    newNode = createNodeFromData(data);
    if (wasActive) {
    	newNode.addClass(ACTIVE_CLASS);
	}
    newNode.appendTo($("#analysis"));
    setupNodes(newNode);
}


function retrieveAndUpdateNodeAppearances(nodeList) {
    for (i=0 ; i<nodeList.length ; ++i) {
        let nodeId = nodeList[i];
        $.ajax({
            url: Urls.node_data(nodeId),
            success: updateNodeAppearance,
        });
    }
}


function addConnectedNode(data) {
	var node = addNewNodeToPage(data);

	var sourceId = data.node_id;
	var targetId = node.attr('id');
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
            var dirty_nodes = [];
            var appearance_update_nodes = [];

            for (var i=0 ; i<data.node_versions.length ; ++i) {
                let nodeData = data.node_versions[i];
                let nodeId = nodeData[0];
                let nodeVersion = nodeData[1];
                let nodeAppearanceVersion = nodeData[2];
                let node = getNode(nodeId);
                let localVersion = node.attr("version_id");
                if (localVersion != nodeVersion) {
                    dirty_nodes.push(nodeId);
                }
                let localAppearanceVersion = node.attr("appearance_version_id");
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
    var aWin = getAnalysisWindow();
    checkAndMarkDirtyNodes(aWin);
    aWin.loadNodeData(node_id);
}
