// @ts-check
// snpdb/templates/snpdb/tags/wiki_tag.html - included once per wiki, so several can share a page

/**
 * A wiki that links citations (the gene symbol wiki) returns their ids from the save, and any
 * [data-wiki-citations] container on the page is redrawn with them - citations.js fills each one in
 * @param {string[]} citationIds
 */
function showWikiCitations(citationIds) {
    const citationDivs = citationIds.map(citationId => $('<div>', {'data-citation-id': citationId}));
    $("[data-wiki-citations]").empty().append(citationDivs);
}

/**
 * @param {{wikiDomId: string, markdown: string, hasWritePermission: boolean, exists: boolean,
 *          className: string, uniqueKeyword: string, uniqueValue: string|number}} options
 */
function initWikiTag(options) {
    const wikiDiv = $("#" + options.wikiDomId);
    const wikiSaveButton = $(".wiki-save", wikiDiv);
    const saveUrl = Urls.wiki_save(options.className, options.uniqueKeyword, options.uniqueValue);
    let editor = null;

    function showEditor() {
        if (editor) {
            return;
        }
        /** @type {Object<string, Function>} */
        const events = {};
        if (options.hasWritePermission) {
            events['change'] = function() { wikiSaveButton.show(); };
        }

        editor = new toastui.Editor({
            el: $(".editSection", wikiDiv).get(0),
            previewStyle: 'vertical',
            initialValue: options.markdown,
            initialEditType: 'wysiwyg',
            height: '300px',
            events: events,
            usageStatistics: false,
        });

        $(".show-wiki-link", wikiDiv).hide();
        $(".editSection", wikiDiv).show();
    }

    wikiSaveButton.click(function() {
        const data = 'markdown=' + encodeURIComponent(editor.getMarkdown());
        const saveMessage = $(".saving-message", wikiDiv);
        wikiSaveButton.hide();
        saveMessage.show();

        /** @param {string} message */
        function showMessage(message) {
            saveMessage.hide();
            $(".saved-message", wikiDiv).html(message).fadeIn('fast').delay(2000).fadeOut('fast');
        }

        $.ajax({
            type: "POST",
            data: data,
            url: saveUrl,
            success: function(data) {
                // Should get back JSON. If we get text, it's a 302 so reload to go to login screen
                if (typeof data === "string") {
                    function loggedInHandler() {
                        console.log("Warning expected JSON response got back: " + data.length);
                    }
                    function loggedOutHandler() {
                        saveMessage.empty();
                        const message = createMessage("error", "You have been logged out and your wiki edit was NOT SAVED. Please copy/paste your change, reload this page, login, then re-enter your wiki entry.");
                        saveMessage.append(message).get(0).scrollIntoView();
                    }
                    checkLoggedIn(loggedInHandler, loggedOutHandler);
                } else {
                    showMessage("Save Successful");
                    if (data.citation_ids) {
                        showWikiCitations(data.citation_ids);
                    }
                }
            },
            error: function() {
                showMessage("Save Failed");
            }
        });
    });

    $(".show-wiki-link", wikiDiv).click(function() {
        showEditor();
        return false;
    });

    if (options.exists) {
        showEditor();
    }
}
