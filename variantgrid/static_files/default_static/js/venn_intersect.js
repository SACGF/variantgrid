VENN_TOGGLE_WIDGET_CLASS = "venn_toggle_widget";

function venn_select(selector, venn_flag) {
		$('.' + VENN_TOGGLE_WIDGET_CLASS, selector).each(function() {
			var toggled = venn_flag & $(this).attr("venn_bit");
			var widget = d3.select(this);
			toggleSelect(widget, !!toggled);
		});
}

function vennAddToggleCallbacks(selector, callback) {
	get_venn_flag = function() {
		var venn_flag = 0;
		$('.' + VENN_TOGGLE_WIDGET_CLASS, selector).each(function() {
			var widget = d3.select(this);
			toggled = widget.attr("toggled");
			if (toggled == "true") {
				venn_flag |= widget.attr("venn_bit");
			}
		});
		return venn_flag;
	};
	toggleColor = function() {
		toggleSelect(d3.select(this));
		venn_flag = get_venn_flag();
		callback(venn_flag);
	};

	$('.' + VENN_TOGGLE_WIDGET_CLASS, selector).each(function() {
		var widget = d3.select(this);
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
	var circle1 = "circle1_" + venn_id;
	var circle2 = "circle2_" + venn_id;

	var radius = w * 0.25;
	var svg = d3.select(selector).append("svg:svg")
	    .attr("width", w)
	    .attr("height", h);
	
	var defs = svg.append("svg:defs");

	var addCirc1 = function(selector) {
		return selector.append("svg:circle")
	    .attr("cx", w*0.36)
	    .attr("cy", h*0.5)
	    .attr("r", radius);
	};

	var addCirc2 = function(selector) {
		return selector.append("svg:circle")
	    .attr("cx", w*0.57)
	    .attr("cy", h*0.5)
	    .attr("r", radius);
	};

	
	var cp = defs.append("svg:clipPath").attr("id", circle1);
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
	    .style("stroke", "gray")
	    .style("stroke-width", 2)
	    .style("fill", "#ffffff");

	var setRing = function(selector) {
		return selector.style("stroke-width", 2)
			.style("fill-opacity", 0)
			.style("stroke", "gray")
			.style("pointer-events", 'none');
	}

	setRing(addCirc1(svg));
	setRing(addCirc2(svg));
}


function venn3(selector, w, h) {
	radius = w * 0.3

	var svg = d3.select(selector).append("svg:svg")
	    .attr("width", w)
	    .attr("height", h);
	
	var defs = svg.append("svg:defs");
	
	defs.append("svg:clipPath")
	    .attr("id", "circle1")
	  	.append("svg:circle")
	    .attr("cx", w*0.36)
	    .attr("cy", h*0.33)
	    .attr("r", radius);
	
	
	defs.append("svg:clipPath")
	    .attr("id", "circle2")
	  	.append("svg:circle")
	    .attr("cx", w*0.57)
	    .attr("cy", h*0.33)
	    .attr("r", radius);
	
	defs.append("svg:clipPath")
	    .attr("id", "circle3")
	  	.append("svg:circle")
	    .attr("cx", w*0.46)
	    .attr("cy", h*0.56)
	    .attr("r", radius);
	
	svg.append("svg:rect")
	    .attr("clip-path", "url(#circle1)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#ff0000")
	    .on("click", function() { alert(2) });
	
	svg.append("svg:rect")
	    .attr("clip-path", "url(#circle2)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#00ff00")
	    .on("click", function() { alert(2) });
	
	svg.append("svg:rect")
	    .attr("clip-path", "url(#circle3)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#0000ff")
	    .on("click", function() { alert(3) });
	
	svg.append("svg:g")
	    .attr("clip-path", "url(#circle1)")
	  	.append("svg:rect")
	    .attr("clip-path", "url(#circle2)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#ffff00")
	    .on("click", function() { alert("4") });
	
	svg.append("svg:g")
	    .attr("clip-path", "url(#circle2)")
	  	.append("svg:rect")
	    .attr("clip-path", "url(#circle3)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#00ffff")
	    .on("click", function() { alert("5") });
	
	svg.append("svg:g")
	    .attr("clip-path", "url(#circle3)")
	  	.append("svg:rect")
	    .attr("clip-path", "url(#circle1)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#ff00ff")
	    .on("click", function() { alert("6") });
	
	svg.append("svg:g")
	    .attr("clip-path", "url(#circle3)")
	  .append("svg:g")
	    .attr("clip-path", "url(#circle2)")
	  .append("svg:rect")
	    .attr("clip-path", "url(#circle1)")
	    .attr("width", w)
	    .attr("height", h)
	    .style("fill", "#ffffff")
	    .on("click", function() { alert("7") });
}

