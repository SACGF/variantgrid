const CitationsManager = (function() {

    function Deferred() {
        const self = this;
        this.promise = new Promise(function(resolve, reject) {
            self.reject = reject
            self.resolve = resolve
        });
    }

    const CitationsManager = function () {
        this.citationToDeferred = {};
        this.citationsToLoad = {};
    }

    CitationsManager.prototype = {

        citationPromise(citationId) {
            const existing = this.citationToDeferred[citationId];
            let deferred = null;
            if (existing) {
                return existing.promise;
            } else {
                deferred = new Deferred();
                this.citationToDeferred[citationId] = deferred;
            }
            this.citationsToLoad[citationId] = true;
            this.debounceRequestData();
            return deferred.promise;
        },

        requestData() {
            const loadIds = Object.keys(this.citationsToLoad);
            this.citationsToLoad = {};
            $.ajax({
                headers: {
                    'Accept' : 'application/json',
                    'Content-Type' : 'application/json'
                },
                //url: this.url + request_ids.join('/'),
                url: Urls.citations_json(loadIds.join('/')),
                type: 'GET',
                error: (call, status, text) => {
                    console.log(text)
                },
                success: (record) => {
                    const requestedIds = {};
                    for (const loadId of loadIds) {
                        requestedIds[loadId] = true;
                    }
                    for (const citationData of record['citations']) {
                        for (const requestingId of citationData.requested_using) {
                            delete requestedIds[requestingId];
                            const deferred = this.citationToDeferred[requestingId];
                            if (deferred) {
                                deferred.resolve(citationData);
                            }
                        }
                    }
                    for (const remainingId of Object.keys(requestedIds)) {
                        const deferred = this.citationToDeferred[remainingId];
                        if (deferred) {
                            deferred.resolve({
                                "title": `Could not load ${remainingId}`
                            })
                        }
                    }
                }
            });
        },

        // debounceRequestData: debounce(this.requestData),

        populate(dom) {
            const $dom = $(dom);
            const citationId = $dom.attr('data-citation-id');
            const prettyId = this.prettyId(citationId);
            $dom.text()
            $dom.addClass('citation');
            $dom.addClass('loading');
            $dom.html([
                $('<span>', {text:`${prettyId}`}),
                $('<span>', {text:' Loading...', class:'text-muted'})
            ]);

            this.citationPromise(citationId).then(data => {
                $dom.removeClass('loading').empty().append(this.renderData(data));
            });
        },

        prettyId(citation_id) {
            return `${citation_id}`.replace(':', ': ').replace('  ', ' ');
        },

        renderData(citation) {
            const citDom = $('<div>');

            if (citation.error) {
                if (citation.id) {
                    citDom.html([
                        $('<div>', {text: citation.id}),
                        $('<div>', {text: citation.error, class: 'text-danger'})
                    ])
                } else {
                    citDom.html($('<div>', {html: citation.error, class: 'text-danger'}));
                }
                return citDom;
            }

            const year = citation.year;

            const authorShort = citation.authors_short;
            const title = citation.title || 'Could not load title';
            const prettyId = this.prettyId(citation.id);

            const linkRow = [];
            if (citation.external_url) {
                linkRow.push($('<a>', {href: citation.external_url, target: '_blank', class:'source external-link', text: prettyId}));
            } else {
                linkRow.push($('<span>', {text: prettyId}));
            }
            if (authorShort) {
                let text = authorShort;
                if (citation.singleAuthor == false) {
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

            const safeId = citation.id.replace(':','_');
            if (!citation.error) {
                $('<a>', {'class':'d-block', 'data-toggle':'ajax-modal', 'data-href':Urls.view_citation_detail(citation.id), 'data-title': prettyId, 'text':'Show Detail'}).appendTo(citDom)
            } else {
                $('<p>', {class: 'abstract', text: citation.abstract && citation.abstract.length ? citation.abstract : 'Could not fetch abstract'}).appendTo(detailContainer);
            }
            return citDom;
        },

        citationDomFor(dbRef, renderNonCitations) {
            let citation_id = null;
            if (typeof(dbRef) == 'string' && dbRef.startsWith('PMID')) {
                citation_id = dbRef;
            } else {
                if (dbRef.internal_id) {
                    citation_id = dbRef.internal_id;
                } else if (['PMID', 'PMCID', 'Bookshelf ID'].indexOf(dbRef.db) !== -1) {
                    citation_id = dbRef.id.replace('\s*', '');
                }
            }
            if (citation_id) {
                return $('<div>', {'data-citation-id': citation_id, html: [
                    $('<span>', {text:`${this.prettyId(dbRef.id)}`}),
                    $('<span>', {text:' Loading...', class:'text-muted'})
                ]});
            } else if (renderNonCitations) {
                let text = dbRef.id;
                if (dbRef.db == "HTTP" || dbRef.db == "HTTPS") {
                    text = `${dbRef.db.toLowerCase()}:${dbRef.idx}`
                }
                const citationDom = $('<div>', {class: 'citation'});
                if (dbRef.id && dbRef.url) {
                    $('<a>', {class: 'no-details', text: text, href: dbRef.url, target: '_blank'}).appendTo(citationDom);
                }
                return citationDom;
            } else {
                return null;
            }
        }
    };

    return CitationsManager;
})();

CitationsManager.dbMigration = {
    "PUBMED": "PMID",
    "PMC": "PMCID",
    "NCBIBOOKSHELF": "Bookshelf ID"
}

CitationsManager.normalizeInPlace = function(dbRef) {
    let dbArray = null;
    if (Array.isArray(dbRef)) {
        dbArray = dbRef;
    } else {
        dbArray = [dbRef];
    }

    for (const dbRef of dbArray) {
        // only need to normalize old records, e.g. the ones with internal ID
        if (dbRef.db && dbRef.idx) {
            const migratedSource = CitationsManager.dbMigration[dbRef.db.toUpperCase()];
            if (migratedSource) {
                let idx = `${dbRef.idx}`;
                if (migratedSource == "PMID") {
                    // no special action
                } else if (migratedSource == "PMCID") {
                    idx = `PMC${idx}`;
                } else if (migratedSource == "Bookshelf ID") {
                    idx = `NBK${idx}`;
                }
                dbRef.db = migratedSource;
                dbRef.idx = idx;
            }
            dbRef.id = `${dbRef.db}:${dbRef.idx}`;
        }
    }
    return dbRef;
}

CitationsManager.prototype.debounceRequestData = debounce(CitationsManager.prototype.requestData);
CitationsManager.defaultManager = new CitationsManager();
