
function AnalysisMessagePoller(url) {
	this.url = url;
	this.update_frequency = 1000;
	this.observed_nodes = {};
	this.update = function() {
//		console.log("update!");
		var nodes = Object.keys(this.observed_nodes);
		if (nodes.length) {
			this.send_requests(nodes);
		}	
	};

	this.send_requests = function(nodes) {
		var that = this;
		var on_success_function = function(data, status, xhr) {
            var nodeStatusList = data.node_status;
            if (nodeStatusList) {
                that.process_nodes(nodeStatusList);
            } else {
                // This should have def been in the JSON response, so perhaps
                // we've been logged out.
                
                that.stop_polling();
                function loggedInHandler() {
                    var message = "AnalysisMessagePoller couldn't handle response.";
					showReloadPageErrorDialog($("#error-dialog"), message, true);
                }
                
                checkLoggedIn(loggedInHandler, loggedOutHandler);
            }
		};

		var data = 'nodes=' + JSON.stringify(nodes);
		$.ajax({
		    type: "GET",
		    data: data,
		    url: this.url,
		    success: on_success_function,
		});
	};

	this.process_nodes = function(nodeStatusList) {
        var nodeStatus;
        var node_id;
		while (nodeStatus = nodeStatusList.pop()) {
			node_id = nodeStatus.id;
			var node_actions = this.observed_nodes[node_id];
			var invalid = nodeStatus.valid == false;
			for (var action in node_actions) {
//				console.log("looking for action: " + action);
				if (invalid || this.can_dispatch[action](nodeStatus)) {
//					console.log("can dispatch");
					
					var callbacks = node_actions[action];
					for (var i=0 ; i<callbacks.length ; ++i) {
						callbacks[i](nodeStatus);
					}
					delete node_actions[action];
				}
			}
		}
		
		// Clean up any now empty nodes / actions
		for (node_id in this.observed_nodes) {
			var action_keys = Object.keys(this.observed_nodes[node_id]);
//			console.log("action keys for " + node_id + " are: " + action_keys + " of length " + action_keys.length);
			if (!action_keys.length) {
				this.delete_node(node_id);
			}
		}
	};

	this.observe_node = function(node_id, action, callback) {
		if (!(action in this.can_dispatch)) {
			console.log("Warning! Unknown action " + action);
			return;
		}
	
		if (!(node_id in this.observed_nodes)) {
			this.observed_nodes[node_id] = {};
		}

		var node_actions = this.observed_nodes[node_id];
		if (!(action in node_actions)) {
			node_actions[action] = [];
		}
		node_actions[action].push(callback);
	};

	this.delete_node = function(node_id) {
		delete this.observed_nodes[node_id];
	};

	this.can_dispatch_count = function(node) {
		return node.count !== null || node.status == 'C';
	};

	this.can_dispatch_ready = function(node) {
		return node.ready;
	};

	this.update_loop = function() {
//		console.log("reschedule_update");
		var that = this;
		
		var pollThis = function() {
			that.update();
			that.update_loop();
		};
		this.timeout = window.setTimeout(pollThis, this.update_frequency);
	};

	this.stop_polling = function() {
		window.clearTimeout(this.timeout);
	};

	this.can_dispatch = {
	   count: this.can_dispatch_count,
	   ready: this.can_dispatch_ready,
	};
}
