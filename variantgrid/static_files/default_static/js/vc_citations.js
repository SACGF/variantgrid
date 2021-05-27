let Citations = (function() {

    let Citations = function(dom) {
        if (dom) {
            this.dom = $(dom);
        }
        this.dbRefs = [];
        this.cache = {};
    };
    
    Citations.REQUESTED = 'requested';
    Citations.FAILED = 'failed';
    Citations.LOADING = 'loading';
    
    Citations.prototype = {
        
        request() {
            let request_ids = [];
            Object.keys(this.cache).forEach(k => {
                if (this.cache[k] == Citations.REQUESTED) {
                    this.cache[k] = Citations.LOADING;
                    request_ids.push(k);
                }
            });
            if (request_ids.length > 0) {
                $.ajax({
                    headers: {
                        'Accept' : 'application/json',
                        'Content-Type' : 'application/json'
                    },
                    //url: this.url + request_ids.join('/'),
                    url: Urls.citations_json(request_ids.join('/')),
                    type: 'GET',
                    error: (call, status, text) => {
                        console.log(text)
                    },
                    success: (record) => {
                        for (let citation of record['citations']) {
                            this.cache[citation.id] = citation;
                            this.cache[`${citation.source}:${citation.citation_id}`] = citation;
                        }
                        this.render();
                    }
                });
            }
        },

        refresh() {
            window.setTimeout(() => {
                let dbRefs = [];
                $('[data-citation]').each((index, dc) => {
                    let dbRef = Citations.parseId($(dc).attr('data-citation'));
                    if (dbRef) {
                        dbRefs.push(dbRef);
                    }
                });
                this.update(dbRefs);
            },1);
        },

        update(dbRefs) {
            this.dbRefs = [];
            let alreadyIncluded = {};
            let plainUrls = {'HTTP':true, 'HTTPS':true, 'FTP':true};
            for (let dbRef of dbRefs) {
                let citationId = Citations.effectiveId(dbRef);
                if (citationId) {
                    if (alreadyIncluded[citationId]) {
                        continue;
                    } else {
                        alreadyIncluded[citationId] = true;
                    }
                    if (!this.cache[citationId]) {
                        this.cache[citationId] = Citations.REQUESTED;
                    }
                    this.dbRefs.push(dbRef);
                } else if (plainUrls[dbRef.db]) {
                    this.dbRefs.push(dbRef);
                }
            }

            this.render();
            this.request();
        },

        renderDbRefState(dbRef, citDom) {
            citDom.empty();
            let effectiveId = Citations.effectiveId(dbRef);

            if (!effectiveId) {
                $('<p>', {html:$('<a>', {class: 'no-details external-link', text: dbRef.url, href: dbRef.url, target: '_blank'})}).appendTo(citDom);
            } else {
                let citation = this.cache[effectiveId];
                if (typeof(citation) == 'string') {
                    $('<p>', {html:$('<a>', {class: 'title external-link', text: dbRef.id || citation, href: dbRef.url, target: '_blank'})}).appendTo(citDom);
                    citDom.addClass(citation);
                } else {
                    citDom.addClass('loaded');

                    let sourceMap = {"PubMed": "PMID"};
                    let source = sourceMap[citation.source] || citation.source;
                    let sourceId = citation.citation_id;
                    let year = citation.year;

                    let authorShort = citation.authors_short;
                    let singleAuthor = null;
                    if (!authorShort && citation.authors) {
                        let authors = citation.authors.split(',');
                        if (authors.length) {
                            authorShort = authors[0];
                            singleAuthor = authors.length === 1;
                        }
                    } else if (authorShort && citation.authors) {
                        singleAuthor = authorShort == citation.authors;
                    }
                    let title = citation.title || 'Could not load title';

                    let linkRow = [];
                    if (citation.citation_link) {
                        linkRow.push($('<a>', {href: citation.citation_link, target: '_blank', class:'source external-link', text: `${source}: ${sourceId}`}));
                    }
                    if (authorShort) {
                        let text = authorShort;
                        if (singleAuthor == false) {
                            text += ' et al';
                        }
                        linkRow.push($('<span>', {class:'author', text: text}));
                    }
                    if (year) {
                        linkRow.push($('<span>', {class:'year', text: year}));
                    }
                    if (title) {
                        linkRow.push($('<span>', {class:'title', text: title}));
                    }
                    let linkDom = $('<span>', {class: 'title-row'});

                    linkDom = linkRow.reduce((prev, current) => {
                        prev.append(current);
                        prev.append((document.createTextNode(' ')));
                        return prev;
                    }, linkDom);
                    linkDom.appendTo(citDom);

                    if (citation.abstract || !singleAuthor || citation.journal) {
                        $('<a>', {class: 'toggle-link d-block', 'data-toggle':"collapse", href:`#detail-${citation.citation_id}`, text: 'Toggle detail'}).appendTo(citDom);
                        let detailContainer = $('<div>', {class: 'collapse', id:`detail-${citation.citation_id}`}).appendTo(citDom);
                        if (citation.journal) {
                            $('<p>', {class: 'journal', text: citation.journal}).appendTo(detailContainer);
                        }
                        if (!singleAuthor) {
                            $('<p>', {class: 'authors', text: citation.authors}).appendTo(detailContainer);
                        }
                        $('<p>', {class: 'abstract', text: citation.abstract && citation.abstract.length ? citation.abstract : 'Could not fetch abstract'}).appendTo(detailContainer);
                    }
                }
            }
        },
        
        render() {
            if (this.dom) {
                this.dom.empty();
                if (this.dbRefs.length) {
                    for (let dbRef of this.dbRefs) {
                        let citDom = $('<div>', {class: 'citation'});
                        this.renderDbRefState(dbRef, citDom);
                        this.dom.append(citDom);
                    }
                } else {
                    $('<span>', {text: 'No citations detected', class: 'no-value'}).appendTo(this.dom);
                }
            } else {
                for (let dbRef of this.dbRefs) {
                    let effectiveId = Citations.effectiveId(dbRef);
                    if (effectiveId) {
                        let citDom = $(`[data-citation='${effectiveId}']`);
                        if (citDom.length) {
                            this.renderDbRefState(dbRef, citDom);
                        } else {
                            console.log(`Could not find HTML element for {dbRef.internal_id}`);
                        }
                    }
                }
            }
        }
    };
    
    return Citations;
})();

// get an ID (that the server can convert to a citation) from a dbRef
Citations.effectiveId = function(dbRef) {
    if (dbRef.internal_id) {
        return dbRef.internal_id;
    } else if (dbRef.db == 'PubMed' && dbRef.idx) {
        return `${dbRef.db}:${dbRef.idx}`;
    } else {
        return null;
    }
};

// convert a string to a dbRef
Citations.parseId = function(id) {
    if (id == null || id == "null") {
        return null;
    }
    let parts = id.split(':');
    if (parts.length == 1) {
        return {"internal_id": parts[0]}
    } else if (parts.length == 2) {
        return {"db": parts[0], "idx": parts[1]}
    }
    return null;
}

Citations.renderDbRef = function(dbRef) {
    let effectiveId = Citations.effectiveId(dbRef);
    let citation = $('<div>', {class: 'citation', 'data-citation': effectiveId});
    if (dbRef.id && dbRef.url) {
        $('<a>', {class: 'no-details', text: dbRef.id, href: dbRef.url, target: '_blank'}).appendTo(citation);
    } else {
        $('<span>', {class: 'loading', text: `loading ${effectiveId}`}).appendTo(citation);
    }
    return citation;
};