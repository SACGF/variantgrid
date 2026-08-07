class AnnotatedFileDownload {
	constructor(selector, pollUrl, fileType) {
		this.selector = $(selector);
		this.pollUrl = pollUrl;
		this.fileType = fileType;
		this.ucFileType = this.fileType.toUpperCase();
	}

	setGenerateDownloadLink() {
		const that = this;
		const generateLink = $(`<a href="#" id="generate-${this.fileType}-download"><div class="icon24 left margin-r-5 ${this.fileType}-icon"></div> Annotated ${this.ucFileType}</a>`);
		generateLink.click(function(event) {
			event.preventDefault();
			that.setPolling(true);
		});
		this.selector.empty();
		this.selector.append(generateLink);
	}

	setDownloadLink(url) {
		this.selector.html(`<a href="${url}"><div class="icon24 left margin-r-5 ${this.fileType}-icon"></div> Annotated ${this.ucFileType}</a>`);
	}

	setError() {
		this.selector.html('<i class="fas fa-xmark"></i> Error retrieving download...');
	}

	setPolling(download) {
		const that = this;
		function downloadFile(data) {
			that.setDownloadLink(data.url);
			if (download) {
				window.location.href = data.url;
			}
		}
		const spinner = $(`<span><i class="fas fa-spinner fa-spin"></i> Preparing download... </span>`);
		const progressIndicator = $("<span></span>");
		spinner.append(progressIndicator);

		function updateProgress(progress) {
			const percent = Math.floor(100 * progress);
			progressIndicator.empty();
			progressIndicator.append(`${percent}% complete`);
		}

		function setError(data) {
			that.setError(data);
		}

		$(this.selector).html(spinner);
		poll_cached_generated_file(this.pollUrl, downloadFile, setError, updateProgress);
	}
}
