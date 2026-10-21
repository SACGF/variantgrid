function groupTerms(terms, minPrefixLen = 15) {
  const lcp = (a, b) => {
    let i = 0, m = Math.min(a.length, b.length);
    while (i < m && a[i] === b[i]) i++;
    return a.slice(0, i);
  };

  const groups = [];
  for (const t of terms) {
    let bestGroup = null;
    let bestPrefix = '';
    for (const g of groups) {
      const p = lcp(g.prefix, t);
      if (p.length >= minPrefixLen && p.length > bestPrefix.length) {
        bestGroup = g;
        bestPrefix = p;
      }
    }
    if (bestGroup) {
      bestGroup.prefix = bestPrefix;
      bestGroup.items.push(t);
    } else {
      groups.push({ prefix: t, items: [t] });
    }
  }
  return groups;
}


function getOntologyTermObj(term_type, term, url) {
    const termSpan = $("<span>", {class: `grid-term-link ${term_type}`, title: term, term: term, term_type: term_type});
    if (url) {
        const termLink = $("<a>", {href: url, target: "_blank", text: term});
        termSpan.append(termLink);
    } else {
        termSpan.text(term);
    }
    return termSpan;
}

function getOntologyTermLinks(term_type, term_list, getUrl) {
    const MIN_PREFIX_LENGTH = 15;
    const NUM_TERMS_TO_ATTEMPT_COLLAPSE = 5;
    const MIN_COLLAPSES = 2;
    let links = '';
    if (term_list) {
        const terms = term_list.split('|');
        const termsContainer = $("<div>", {class: "ontology-terms-container"});
        if (terms.length >= NUM_TERMS_TO_ATTEMPT_COLLAPSE) {
            const groups = groupTerms(terms, MIN_PREFIX_LENGTH);
            for (const g of groups) {
                const fullTerms = [];
                for (const term of g.items) {
                    let url = null;
                    if (getUrl) {
                        url = getUrl(term_type, term);
                    }
                    fullTerms.push(getOntologyTermObj(term_type, term, url));
                }

                if (g.items.length >= MIN_COLLAPSES && g.prefix.length >= MIN_PREFIX_LENGTH) {
                    const prefix = g.prefix.replace(/[ ,]+$/, '');
                    const collapsedLabel = `${prefix} [${g.items.length} matches]`;
                    const collapsedTerm = getOntologyTermObj(term_type, collapsedLabel);
                    collapsedTerm.addClass('collapsed-term');
                    collapsedTerm.attr("title", "click to expand full terms");

                    const group = $("<span>", {class: "term-group"});
                    const fullTermsContainer = $("<span>", {
                        class: "full-terms",
                        style: "display:none"
                    });
                    fullTermsContainer.append(...fullTerms);

                    group.append(collapsedTerm);
                    group.append(fullTermsContainer);
                    termsContainer.append(group);
                } else {
                    termsContainer.append(...fullTerms);
                }
            }
        } else {
            // Just add the terms directly
            for (const term of terms) {
                let url = null;
                if (getUrl) {
                    url = getUrl(term_type, term);
                }
                const termObj = getOntologyTermObj(term_type, term, url);
                termObj.addClass("grid-term-link");
                termsContainer.append(termObj);
            }
        }
        links = termsContainer.prop("outerHTML");
    }
    return links;
}


// DataTables renderer - kwargs: {termType}, e.g. 'HP'
function ontologyTermsRenderer(data, type, row, ctx) {
    if (type !== 'display') {
        return data ? data.split('|').join(', ') : '';
    }
    return getOntologyTermLinks(ctx.kwargs.termType, data);
}


// This needs to be set in document.ready()
// '.ontology-terms-container .collapsed-term'
function expandCollapsedOntologyTerm() {
    const group = $(this).closest('.term-group');
    group.find('.collapsed-term').hide();
    group.find('.full-terms').show();
}
