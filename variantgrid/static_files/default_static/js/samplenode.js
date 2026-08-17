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


function sourceBadge(node, source) {
	// Why a specimen node returns fewer rows than it used to - a VCF was archived - is a question
	// only the canvas can answer, so show how many VCFs the node is actually reading
	const badge = $("<div class='source-badge'/>");
	badge.attr("title", source['level'] + ": " + source['vcf_names'].join(", "));
	badge.css({
		position: 'absolute',
		bottom: 0,
		right: 0,
		'z-index': 40,
		'font-size': '10px',
		'line-height': '12px',
		background: '#ffffff',
		border: '1px solid #888',
		'border-radius': '3px',
		padding: '0px 2px',
	});
	badge.html("<i class='fas fa-file-alt'></i>&times;" + source['vcf_count']);
	badge.appendTo(node);
}


function sampleNodeUpdateState(args) {
	// remove existing SVG
	$('svg', this).remove();
	$('.source-badge', this).remove();
	const sideLength = SIDE_LENGTH;

	const source = args['source'];
	if (source && source['vcf_count']) {
		sourceBadge(this, source);
	}

	const patient = args['patient'];
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
