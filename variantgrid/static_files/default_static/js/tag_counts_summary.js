// Pill + count toggles that filter the grid below them - @see the tag_counts_summary template tag

function getTagCountsSummarySelected(container) {
	return $(".summary-count.selected", container).map(function() {
		return $(this).attr("data-tag");
	}).get();
}

function _tagCountsSummaryUpdated(container) {
	const selected = getTagCountsSummarySelected(container);
	container.toggleClass("filtered", selected.length > 0);
	// Disabled rather than hidden - the pills would shift along as the button came and went
	$(".clear-tags", container).prop("disabled", selected.length === 0);
	$(".all-tags", container).prop("disabled", selected.length === $(".summary-count", container).length);
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
