SIDE_LENGTH = 60;
DEFAULT_SHADOW_COLOR = "#aaa";


function maleSVG(svg, sideLength) {
	const sideSize = sideLength * 0.88;
	const offset = (sideLength - sideSize) / 2;
	const translate = "translate(" + offset + "," + offset + ")";

	svg.append("svg:rect")
	   .attr("class", "sample-shadow")
	   .attr("width", sideSize)
	   .attr("height", sideSize)
	   .attr("transform", translate)
	   .style("fill", DEFAULT_SHADOW_COLOR)
	   .style("filter", "url(#dropshadow)");

	svg.append("svg:rect")
	   .attr("class", "sample-node")
	   .attr("width", sideSize)
	   .attr("height", sideSize)
	   .attr("transform", translate)
	   .style("fill", "#ffffff");
}	


function femaleSVG(svg, sideLength) {
	const radius = sideLength * 0.44;

	svg.append("svg:circle")
	   .attr("class", "sample-shadow")
	   .attr("cx", sideLength*0.5)
	   .attr("cy", sideLength*0.5)
	   .attr("r", radius)
	   .style("fill", DEFAULT_SHADOW_COLOR)
	   .style("filter", "url(#dropshadow)");

	svg.append("svg:circle")
	   .attr("class", "sample-node")
	   .attr("cx", sideLength*0.5)
	   .attr("cy", sideLength*0.5)
	   .attr("r", radius)
	   .style("fill", "#ffffff");
}	

function addDeceasedStroke(svg, sideLength) {
	svg.append("svg:line")
	   .attr("class", "sample-node")
	   .attr("x1", sideLength)
	   .attr("y1", 0)
	   .attr("x2", 0)
	   .attr("y2", sideLength)
	   .style("fill", "#ffffff")
	   .style("filter", "url(#dropshadow)");
}


function createSampleNode() {
	const sampleNode = $("<div/>").addClass("window");
	const nodeOverlay = $("<div class='node-overlay' />");
	const overlay_style = {
		width: '100%',
		height: '100%',
		position: 'absolute',
		top: 0,
		left: 0,
		'line-height': SIDE_LENGTH + 'px',
		'z-index' : 30,
	};
	nodeOverlay.css(overlay_style);

	const span = $("<span class='node-name'></span>");
	span.css({padding: '8px', display: 'inline-block', 'vertical-align' : 'middle', 'line-height' : '1em'});
	nodeOverlay.append(span);
	$("<div />", {class: "node-color-overlay"}).appendTo(nodeOverlay);
	sampleNode.append(nodeOverlay);
	sampleNode[0].updateState = sampleNodeUpdateState;
	return sampleNode;
}


function sampleNodeUpdateState(args) {
	// remove existing SVG
	$('svg', this).remove();
	const sideLength = SIDE_LENGTH;

	const patient = args['patient']
	if (patient) {
		const sex = patient['sex'];
		const deceased = patient['deceased'];

		const svg = d3.select(this).append("svg:svg")
			.attr("version", "1.1")
			.attr("width", sideLength)
			.attr("height", sideLength);

		if (sex == 'F') {
			femaleSVG(svg, sideLength);
		} else {
			maleSVG(svg, sideLength);
		}

		if (deceased) {
			addDeceasedStroke(svg, sideLength);			
		}
	}
}
