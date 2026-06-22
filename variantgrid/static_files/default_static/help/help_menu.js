var helpMenu = [
	["About", "about.html"],
    ["Terms & Conditions", "terms.html"],
];

function populateSideMenu(baseUrl, menu)	{
	for (var i=0;i<helpMenu.length ; ++i) {
		var row = helpMenu[i];
		var text = row[0];
		var url = baseUrl + row[1];
		menu.append('<li><a href="' + url + '">' + text + '</a></li>');					
	}
}

