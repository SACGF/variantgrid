// I have problems with formatters when re-importing the grid. Manage this via only importing if not already exists 
if (!$(document).jqGrid) {
    //console.log("Importing JQGrid stuff");
    var SCRIPTS = ['/static/js/lib/free-jqgrid/jquery.jqgrid.min.js'];
    var CSS = [ '/static/js/lib/free-jqgrid/css/ui.jqgrid.css' ];

    var headSelector = $("head");

    for (var i=0 ; i<CSS.length ; ++i) {
        var fileref = document.createElement("link");
        fileref.setAttribute("rel", "stylesheet");
        fileref.setAttribute("type", "text/css");
        fileref.setAttribute("href", CSS[i]);
        headSelector.append(fileref);
    }

    for (var i=0 ; i<SCRIPTS.length ; ++i) {
        var fileref = document.createElement('script');
        fileref.setAttribute("type","text/javascript");
        fileref.setAttribute("src", SCRIPTS[i]);
        headSelector.append(fileref);
    }
}
