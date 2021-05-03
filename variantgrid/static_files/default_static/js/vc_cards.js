const VCCard = (function() {
    'use strict';

    const VCCard = function (options) {
        this.classifications = options.classifications;
        this.eKeys = options.eKeys;
        this.isChild = options.isChild;
    };

    VCCard.RE_CHGVS = /^(.+?)[(](.+?)[)](.+?)$/;

    VCCard.prototype = {

        formatChgvs: function(cHgvs) {
            if (cHgvs) {
                let match = cHgvs.match(VCCard.RE_CHGVS);
                if (match) {
                    return $('<div>', {
                        class: 'variant-coordinate text-monospace', style: 'font-size:10px', html: [
                            $('<span>', {text: match[1]}),
                            $('<span>', {text: '('}),
                            $('<span>', {text: match[2], class: 'text-secondary'}),
                            $('<span>', {text: ')'}),
                            $('<span>', {text: match[3]})
                        ]
                    });
                } else {
                    return $('<div>', {class: 'variant-coordinate text-monospace font-small', text: cHgvs});
                }
            }
            return String(cHgvs);
        },
        formatClinSig: function(clinSig) {
            let formatted = this.eKeys.key(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(clinSig) || "Unclassified";
            let clinSigKey = (clinSig || '').toLowerCase();
            // needs to be a span due to css being on span
            return $('<div>', {class: `cs-card cs-${clinSigKey}`, text: formatted});
        },
        any: function() {
            return this.classifications[0];
        },
        valuesFor(key) {
            let allValues = new Set();
            for (let classification of this.classifications) {
                allValues.add(classification[key]);
            }
            let allValuesList = Array.from(allValues);
            return allValuesList.sort();
        },
        conditionTerms() {
            let idToTerm = {};
            let plainText = new Set();
            for (let classification of this.classifications) {
                let condition = classification.condition;
                if (condition.resolved_terms) {
                    for (let term of condition.resolved_terms) {
                        idToTerm[term.term_id] = term;
                    }
                } else {
                    if (condition.display_text) {
                        plainText.add(condition.display_text);
                    }
                }
            }
            let blocks = Object.values(idToTerm)
                .sort((a, b) => a.term_id.localeCompare(b.term_id))
                .map(term => {return {resolved_terms: [term]};})
                .map(VCForm.format_condition)
                .map(a => {a.addClass('d-block'); return a})
            blocks = blocks.concat(Array.from(plainText).sort());
            return blocks;
        },
        renderValues(values) {
            let container = $('<div>', {class:'mt-2'});
            for (let value of values) {
                if (value === null || value === "null") {
                    // $('<div>', {text: '-'}).appendTo(container);
                } else if (_.isString(value)) {
                    $('<div>', {text: value}).appendTo(container);
                } else {
                    value.appendTo(container);
                }
            }
            return container;
        },

        render: function() {
            let wrapper = $('<div>', {class: 'col-12 col-lg-4 col-md-6 mb-4'});
            let clinSig = (this.any().clinical_significance || '').toLowerCase();
            let card = $('<div>', {class: 'card'}).appendTo(wrapper);
            if (this.isChild) {
                card.addClass('fade-in-scale');
            }
            // clin sig is the card header
            let cardHeader = $('<div>', {class: 'card-header'}).appendTo(card);
            let cardBody = $('<div>', {class: `card-body cs-${clinSig}`}).appendTo(card);
            this.formatClinSig(this.any().clinical_significance).appendTo(card);

            let cardFooter = $('<a>', {class: 'card-footer hover-link'}).appendTo(card);

            let clinGrouping = this.any().clinical_grouping;
            let headerParts = [];
            if (clinGrouping !== 'default') {
                if (clinGrouping === 'unshared') {
                    card.css('opacity', 0.5);
                }
                headerParts.push($('<span>', {text: clinGrouping.toUpperCase(), class: 'font-weight-bold text-secondary'}));
                headerParts.push(" ");
            }
            headerParts.push(this.valuesFor('org').join(", "));
            if (this.classifications.length > 1) {
                headerParts.push($('<span>', {class: 'badge badge-light ml-2 text-monospace', text: 'x' + String(this.classifications.length)}));
            }

            $('<div>', {html:headerParts}).appendTo(cardHeader);
            this.renderValues(this.valuesFor('c_hgvs').map(this.formatChgvs.bind(this))).appendTo(cardBody);

            this.renderValues(this.conditionTerms()).appendTo(cardBody);
            if (this.classifications.length === 1) {
                $('<div>', {text:'View classification'}).appendTo(cardFooter);
                cardFooter.attr('href', Urls.view_classification(this.any().cid));
            } else {
                $('<div>', {class: 'toggle-link', html: `Expand ${this.classifications.length} records`}).appendTo(cardFooter);
                cardFooter.click(() => {this.expand();});
            }

            this.dom = wrapper;
            return wrapper;
        },

        expand: function() {
            this.dom.after(this.classifications.map(cm => new VCCard({classifications: [cm], eKeys: this.eKeys, isChild:true}).render()));
            this.dom.remove();
        }
    };

    return VCCard;

})();

VCCard.init = function(options) {
    'use strict';

    const classifications = options.classifications;
    const dom = options.dom;
    const clinGroupingOrder = Object.freeze({
       "unshared": 1,
       "default": -1
    });

    const byClinGrouping = function(data) {
        return Object.entries(_.groupBy(data, cm => cm.clinical_grouping)).sort((a, b) => {
            let aClin = a[0];
            let bClin = b[0];
            if (aClin === bClin) {
                return 0;
            }
            return (clinGroupingOrder[aClin] || 0) - (clinGroupingOrder[bClin] || 0);
        }).map(group => group[1])
    };

    const clinSigComparator = options.eKeys.key(SpecialEKeys.CLINICAL_SIGNIFICANCE).comparator();
    const byClinSig = function(data) {
        return Object.entries(_.groupBy(data, cm => cm.clinical_significance))
            .sort((a, b) => clinSigComparator(a[0], b[0]))
            .map(group => group[1]);
    };

    const byOrg = function(data) {
        return Object.entries(_.groupBy(data, cm => cm.org))
            .sort((a, b) => a[0].localeCompare(b[0]))
            .map(group => group[1]);
    }

    // FIXME sort groupings
    for (let group1 of byClinSig(classifications)) {
        for (let group2 of byClinGrouping(group1)) {
            for (let group3 of byOrg(group2)) {
                dom.append(new VCCard({classifications: group3, eKeys: options.eKeys}).render());
            }
        }
    }
};