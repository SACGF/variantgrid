freq = 1000;

function clearGraph(graph_selector) {
	graph_selector.empty();
	graph_selector.addClass('generated-graph graph-loading');
}


function poll_graph_status(graph_selector, poll_url, delete_url) {
	$.getJSON(poll_url, function(data) {
		if (data.status == "SUCCESS") {
			const img = new Image();
			$(img).hide();
	        $(img).attr('src', data.url);
            graph_selector.removeClass('graph-loading').append(img);
            $(img).fadeIn();
		} else if (data.status == 'FAILURE') {

			// Each graph retries itself - a page can show several (the node Summary tab's boxplots)
			const retryGenerateGraph = function (event) {
				event.preventDefault();
				event.stopPropagation();
				$.ajax({
					type: "POST",
					data: {cgf_id: data.cgf_id},
					url: delete_url,
					success: function () {
						graph_selector.empty();
						poll_graph_status(graph_selector, poll_url, delete_url);
					},
				});
			};

            graph_selector.removeClass('graph-loading').addClass("graph-failure").one("click", function(){
				const retryLink = $("<a>", {href: "#", text: "try again?"}).click(retryGenerateGraph);
				$(this).removeClass("graph-failure").empty().append(
					$("<p>", {text: "Graph generation failed with error message:"}),
					$("<b>", {text: data.exception}),
					$("<p>").append("Maybe you can ", retryLink));
            });
		} else {
			const retry_func = function () {
				poll_graph_status(graph_selector, poll_url, delete_url);
			};
			window.setTimeout(retry_func, freq);
		}
	});
}


function poll_cached_generated_file(poll_url, success_func, failure_func, update_func) {
	$.getJSON(poll_url, function (data) {
		if (data.status == "SUCCESS") {
			success_func(data);
		} else if (data.status == 'FAILURE') {
			failure_func(data);
		} else {
			if (update_func) {
				update_func(data.progress);
			}

			const retry_func = function () {
				poll_cached_generated_file(poll_url, success_func, failure_func, update_func);
			};
			window.setTimeout(retry_func, freq);
		}
	}).fail(failure_func);
}
