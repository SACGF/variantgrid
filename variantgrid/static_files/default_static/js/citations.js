let CitationsManager = (function() {

    function Deferred() {
        let self = this;
        this.promise = new Promise(function(resolve, reject) {
            self.reject = reject
            self.resolve = resolve
        });
    }

    let CitationsManager = function () {
        this.citationToDeferred = {};
        this.citationsToLoad = {};
    }

    CitationsManager.prototype = {

        citationPromise(citationId) {
            let existing = this.citationToDeferred[citationId];
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
            let loadIds = Object.keys(this.citationsToLoad);
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
                    let requestedIds = {};
                    for (let loadId of loadIds) {
                        requestedIds[loadId] = true;
                    }
                    for (let citationData of record['citations']) {
                        for (let requestingId of citationData.requested_using) {
                            delete requestedIds[requestingId];
                            let deferred = this.citationToDeferred[requestingId];
                            if (deferred) {
                                deferred.resolve(citationData);
                            }
                        }
                    }
                    for (let remainingId of Object.keys(requestedIds)) {
                        let deferred = this.citationToDeferred[remainingId];
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
            let $dom = $(dom);
            $dom.addClass('citation');
            $dom.addClass('loading');
            let citationId = $dom.attr('data-citation-id');
            this.citationPromise(citationId).then(data => {
                $dom.removeClass('loading').empty().append(this.renderData(data));
            });
        },

        renderData(citation) {
            let citDom = $('<div>');
            let year = citation.year;

            let authorShort = citation.authors_short;
            let title = citation.title || 'Could not load title';

            let linkRow = [];
            if (citation.external_url) {
                linkRow.push($('<a>', {href: citation.external_url, target: '_blank', class:'source external-link', text: `${citation.id}`}));
            } else {
                linkRow.push($('<span>', {text: `${citation.id}`}));
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

            if (citation.abstract || !citation.singleAuthor || citation.journal) {
                $('<a>', {class: 'toggle-link d-block', 'data-toggle':"collapse", href:`#detail-${citation.id}`, text: 'Toggle detail'}).appendTo(citDom);
                let detailContainer = $('<div>', {class: 'collapse', id:`detail-${citation.id}`}).appendTo(citDom);
                if (citation.journal) {
                    $('<p>', {class: 'journal', text: citation.journal}).appendTo(detailContainer);
                }
                if (!citation.singleAuthor) {
                    $('<p>', {class: 'authors', text: citation.authors}).appendTo(detailContainer);
                }
                $('<p>', {class: 'abstract', text: citation.abstract && citation.abstract.length ? citation.abstract : 'Could not fetch abstract'}).appendTo(detailContainer);
            }
            return citDom;

        }
    };

    return CitationsManager;
})();

CitationsManager.prototype.debounceRequestData = debounce(CitationsManager.prototype.requestData);

CitationsManager.defaultManager = new CitationsManager();
