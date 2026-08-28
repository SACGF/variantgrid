VENN_TOGGLE_WIDGET_CLASS = "venn_toggle_widget";
// Rings draw in currentColor, so they match whatever the container is coloured - see analysis_nodes.css

function venn_select(selector, venn_flag) {
		$('.' + VENN_TOGGLE_WIDGET_CLASS, selector).each(function() {
			const toggled = venn_flag & $(this).attr("venn_bit");
			const widget = svgSelect(this);
			toggleSelect(widget, !!toggled);
		});
}

function vennAddToggleCallbacks(selector, callback) {
	get_venn_flag = function() {
		let venn_flag = 0;
		$('.' + VENN_TOGGLE_WIDGET_CLASS, selector).each(function() {
			const widget = svgSelect(this);
			toggled = widget.attr("toggled");
			if (toggled == "true") {
				venn_flag |= widget.attr("venn_bit");
			}
		});
		return venn_flag;
	};
	toggleColor = function() {
		toggleSelect(svgSelect(this));
		venn_flag = get_venn_flag();
		callback(venn_flag);
	};

	$('.' + VENN_TOGGLE_WIDGET_CLASS, selector).each(function() {
		const widget = svgSelect(this);
		widget.on("click", toggleColor);
	});
}

// Either set to value (if provided) or toggle (if not)
function toggleSelect(select, value) {
	toggled = value;
	if (toggled == null) {
		toggled = select.attr('toggled');
		toggled = toggled != "true";
	}
	select.attr('toggled', toggled)
	.style("fill", toggled? "#FF0000" : "#FFFFFF");
}

venn_id = 0;

function venn2(selector, w, h) {
	venn_id++;
	const circle1 = "circle1_" + venn_id;
	const circle2 = "circle2_" + venn_id;

	const radius = w * 0.25;
	const svg = svgSelect(selector).append("svg:svg")
	    .attr("width", w)
	    .attr("height", h);
	
	const defs = svg.append("svg:defs");

	const addCirc1 = function(selector) {
		return selector.append("svg:circle")
	    .attr("cx", w*0.36)
	    .attr("cy", h*0.5)
	    .attr("r", radius);
	};

	const addCirc2 = function(selector) {
		return selector.append("svg:circle")
	    .attr("cx", w*0.57)
	    .attr("cy", h*0.5)
	    .attr("r", radius);
	};

	
	let cp = defs.append("svg:clipPath").attr("id", circle1);
	addCirc1(cp);
	
	cp = defs.append("svg:clipPath").attr("id", circle2);
	addCirc2(cp);

	svg.append("svg:rect")
		.attr("class", VENN_TOGGLE_WIDGET_CLASS)
	    .attr("venn_bit", 1)
	    .attr("toggled", false)
	    .attr("clip-path", "url(#" + circle1 + ")")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#ffffff");
	
	svg.append("svg:rect")
		.attr("class", VENN_TOGGLE_WIDGET_CLASS)
	    .attr("venn_bit", 4)
	    .attr("toggled", false)
	    .attr("clip-path", "url(#" + circle2 + ")")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#ffffff");
	
	svg.append("svg:g")
	    .attr("clip-path", "url(#" + circle1 + ")")
	    .append("svg:rect")
	    .attr("clip-path", "url(#" + circle2 + ")")
		.attr("class", VENN_TOGGLE_WIDGET_CLASS)
	    .attr("venn_bit", 2)
	    .attr("toggled", false)
	    .attr("width", w)
	    .attr("height", h)
	    .style("stroke", "currentColor")
	    .style("stroke-width", 2)
	    .style("fill", "#ffffff");

	const setRing = function(selector) {
		return selector.style("stroke-width", 2)
			.style("fill-opacity", 0)
			.style("stroke", "currentColor")
			.style("pointer-events", 'none');
	};

	setRing(addCirc1(svg));
	setRing(addCirc2(svg));
}
