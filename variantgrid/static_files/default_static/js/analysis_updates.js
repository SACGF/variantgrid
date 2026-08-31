
function AnalysisMessagePoller(node_status_url, task_status_url) {
	const that = this;
	this.node_status_url = node_status_url;
	this.task_status_url = task_status_url;
	this.update_frequency = 1000;
	this.observed_nodes = {};
	this.counts_watches = {};
	this.update = function() {
//		console.log("update!");
		const nodes = Object.keys(this.observed_nodes);
		if (nodes.length) {
			this.send_requests(nodes);
		}
	};

	this.send_requests = function(nodes) {
		const that = this;
		const on_node_status_success_function = function (data, status, xhr) {
			const nodeStatusList = data.node_status;
			if (nodeStatusList) {
				that.process_nodes(nodeStatusList);
			} else {
				// This should have def been in the JSON response, so perhaps
				// we've been logged out.

				that.stop_polling();

				function loggedInHandler() {
					const message = "AnalysisMessagePoller couldn't handle response.";
					showReloadPageErrorDialog(message, true);
				}

				checkLoggedIn(loggedInHandler, loggedOutHandler);
			}
		};

		const error_function = function(data) {
			that.stop_polling();
		};

		const data = 'nodes=' + JSON.stringify(nodes);
		$.ajax({
		    type: "GET",
		    data: data,
		    url: this.node_status_url,
		    success: on_node_status_success_function,
			error: error_function,
		});

        if (this.task_status_url) {
            const on_task_status_success_function = function (data, status, xhr) {
                const analysisTasksSpan = $("#analysis-tasks");
                const errorSpan = $("#analysis-tasks-error", analysisTasksSpan);
                const analysisActiveSpan = $("#analysis-tasks-active", analysisTasksSpan);
                const analysisQueueSpan = $("#analysis-tasks-queued", analysisTasksSpan);

                if (data.error) {
                    errorSpan.text(data.error);
                    analysisActiveSpan.empty();
                    analysisQueueSpan.empty();
                    analysisTasksSpan.show();
                } else {
                    errorSpan.hide();
                    const numActive = data["ACTIVE"] || 0;
                    const numQueued = data["QUEUED"] || 0;
                    if (numActive || numQueued) {
                        analysisActiveSpan.text(numActive);
                        if (numQueued) {
                            analysisQueueSpan.text(` (Q: ${numQueued})`);
                        } else {
                            analysisQueueSpan.empty();
                        }
                        const title = `Active tasks: ${numActive}, Queued: ${numQueued}`;
                        analysisTasksSpan.attr("title", title);
                        analysisTasksSpan.show();
                    } else {
                        analysisTasksSpan.hide();
                    }
                }
            };

            $.ajax({
                type: "GET",
                url: this.task_status_url,
                success: on_task_status_success_function,
            });
        }
	};

	this.process_nodes = function(nodeStatusList) {
		let nodeStatus;
		let node_id;
		while (nodeStatus = nodeStatusList.pop()) {
			node_id = nodeStatus.id;
			const node_actions = this.observed_nodes[node_id];
			const invalid = nodeStatus.valid === false;
			for (const action in node_actions) {
//				console.log("looking for action: " + action);
				if (invalid || this.can_dispatch[action](nodeStatus)) {
//					console.log("can dispatch");

					const callbacks = node_actions[action];
					for (let i=0 ; i<callbacks.length ; ++i) {
						callbacks[i](nodeStatus);
					}
					delete node_actions[action];
				}
			}
		}
		
		// Clean up any now empty nodes / actions
		for (node_id in this.observed_nodes) {
			const action_keys = Object.keys(this.observed_nodes[node_id]);
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

		const node_actions = this.observed_nodes[node_id];
		if (!(action in node_actions)) {
			node_actions[action] = [];
		}
		node_actions[action].push(callback);
	};

	this.delete_node = function(node_id) {
		delete this.observed_nodes[node_id];
		delete this.counts_watches[node_id];
	};

	// Tag counts are recalculated in a task after the node has loaded, so they change without the node
	// version moving. Watch counts_modified to know when that landed - see nodes_status
	this.observe_counts_changed = function(node_id, baseline, timeout_ms, callback) {
		this.counts_watches[node_id] = {baseline: baseline, expires: Date.now() + timeout_ms};
		this.observe_node(node_id, "counts_changed", callback);
	};

	this.can_dispatch_count = function(node) {
		return node.count !== null || node.status === 'C';
	};

	this.can_dispatch_ready = function(node) {
		return node.ready;
	};

	this.can_dispatch_counts_changed = function(node) {
		const watch = that.counts_watches[node.id];
		if (!watch) {
			return true;
		}
		// Give up on the deadline and draw what we have - a node that errored or was deleted will
		// never write new counts, and we don't want to poll for it forever
		if (node.counts_modified != watch.baseline || Date.now() > watch.expires) {
			delete that.counts_watches[node.id];
			return true;
		}
		return false;
	};

	this.update_loop = function() {
//		console.log("reschedule_update");
		const that = this;

		const pollThis = function () {
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
	   counts_changed: this.can_dispatch_counts_changed,
	};
}
