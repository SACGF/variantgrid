const helpMenu = [
    ["About", "about.html"],
    ["Terms & Conditions", "terms.html"],
    ["FAQ", "faq.html"],
];

function populateSideMenu(baseUrl, menu)    {
    for (let i=0;i<helpMenu.length ; ++i) {
        const row = helpMenu[i];
        const text = row[0];
        const url = baseUrl + row[1];
        menu.append('<li><a href="' + url + '">' + text + '</a></li>');                 
    }
}