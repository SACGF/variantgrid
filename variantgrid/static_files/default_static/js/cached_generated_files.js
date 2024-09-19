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

			// Capture passed params and make global func so link generated below will work.
			const retry_generate_graph = function (cgf_id) {
				const clear_and_reload = function () {
					graph_selector.empty();
					poll_graph_status(graph_selector, poll_url, delete_url);
				};
				$.ajax({
					type: "POST",
					data: 'cgf_id=' + cgf_id,
					url: delete_url,
					success: clear_and_reload,
				});
			};
			window.retry_generate_graph = retry_generate_graph; 

            graph_selector.removeClass('graph-loading').addClass("graph-failure").click(function(){
				$(this).removeClass("graph-failure");
				const retryHTML = "<p>Maybe you can <a href='javascript:retry_generate_graph(" + data.cgf_id + ")'>try again?</a></p>";
				$(this).html("<p>Graph generation failed with error message:</p><b>" + data.exception + "</b>");
				$(this).append(retryHTML);
            });
		} else {
			const retry_func = function () {
				poll_graph_status(graph_selector, poll_url, delete_url);
			};
			window.setTimeout(retry_func, freq);
		}
	});
}


function poll_cached_generated_file(poll_url, success_func, failure_func) {
	$.getJSON(poll_url, function (data) {
		if (data.status == "SUCCESS") {
			success_func(data);
		} else if (data.status == 'FAILURE') {
			failure_func(data);
		} else {
			const retry_func = function () {
				poll_cached_generated_file(poll_url, success_func, failure_func);
			};
			window.setTimeout(retry_func, freq);
		}
	}).fail(failure_func);
}


class AnnotatedFileDownload {
	constructor(selector, pollUrl, fileType) {
		this.selector = $(selector);
		this.pollUrl = pollUrl;
		this.fileType = fileType;
		this.ucFileType = this.fileType.toUpperCase()
	}

	setGenerateDownloadLink() {
		// console.log("setGenerateDownloadLink");
		const that = this;
		const generateLink = $(`<a href="#" id="generate-${this.fileType}-download"><div class="icon24 left margin-r-5 ${this.fileType}-icon"></div> Annotated ${this.ucFileType}</a>`);
		generateLink.click(function(event) {
			event.preventDefault();
			that.setPolling(true);
		});
		this.selector.empty();
		this.selector.append(generateLink);
		// console.log("/setGenerateDownloadLink");
	}

	setDownloadLink(url) {
		// console.log("setDownloadLink");
		this.selector.html(`<a href="${url}"><div class="icon24 left margin-r-5 ${this.fileType}-icon"></div> Annotated ${this.ucFileType}</a>`);
	}

	setError() {
		// console.log("setError");
		this.selector.html('<i class="fas fa-xmark"></i> Error retrieving download...');
	}

	setPolling(download) {
		// console.log("setPolling");
		let that = this;
		function downloadFile(data) {
			that.setDownloadLink(data.url);
			if (download) {
				window.location.href = data.url;
			}
		}
		$(this.selector).html(`<span><i class="fas fa-spinner fa-spin"></i> Preparing download...</span>`);
		poll_cached_generated_file(this.pollUrl, downloadFile, this.setError)
	}
}
