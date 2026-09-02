// Pill + count toggles that filter the grid below them - @see the tag_counts_summary template tag

function getTagCountsSummarySelected(container) {
	return $(".summary-count.selected", container).map(function() {
		return $(this).attr("data-tag");
	}).get();
}

function _tagCountsSummaryUpdated(container) {
	const selected = getTagCountsSummarySelected(container);
	container.toggleClass("filtered", selected.length > 0);
	$(".clear-tags", container).prop("hidden", selected.length === 0);
	return selected;
}

/* onChange(selectedTagIds) fires after every toggle - multiple tags mean "carries any of these" */
function setupTagCountsSummary(container, onChange) {
	container = $(container);
	const changed = function() {
		onChange(_tagCountsSummaryUpdated(container));
	};

	container.on("click", ".summary-count", function(event) {
		event.preventDefault();
		$(this).toggleClass("selected");
		changed();
	});
	$(".all-tags", container).click(function() {
		$(".summary-count", container).addClass("selected");
		changed();
	});
	$(".clear-tags", container).click(function() {
		$(".summary-count", container).removeClass("selected");
		changed();
	});
	_tagCountsSummaryUpdated(container);
}

/* The selection as a node count / extra_filters label - "default" when nothing is picked.
   @see TagFilter, which parses it back into tags server side */
function getTagCountsSummaryFilter(container) {
	const labels = $(".summary-count.selected", container).map(function() {
		return $(this).attr("data-tag-label");
	}).get();
	return labels.length ? labels.join(",") : "default";
}

// Deselect everything without firing onChange - for pages that clear their filters some other way
function clearTagCountsSummary(container) {
	container = $(container);
	$(".summary-count", container).removeClass("selected");
	_tagCountsSummaryUpdated(container);
}

/* Fill in counts ({tagId: count}) the server left blank, hiding the tags this view has none of -
   a selected tag stays put so you can always toggle it back off */
function setTagCountsSummaryCounts(container, counts) {
	$(".summary-count", container).each(function() {
		const summaryCount = $(this);
		const count = counts[summaryCount.attr("data-tag")];
		summaryCount.toggle(Boolean(count) || summaryCount.hasClass("selected"));
		$(".count", summaryCount).text(count === undefined ? "" : count.toLocaleString());
	});
}
