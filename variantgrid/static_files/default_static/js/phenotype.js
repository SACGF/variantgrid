function groupTerms(terms, minPrefixLen = 15) {
  const lcp = (a, b) => {
    let i = 0, m = Math.min(a.length, b.length)
    while (i < m && a[i] === b[i]) i++
    return a.slice(0, i)
  }

  const groups = []
  for (const t of terms) {
    let bestGroup = null
    let bestPrefix = ''
    for (const g of groups) {
      const p = lcp(g.prefix, t)
      if (p.length >= minPrefixLen && p.length > bestPrefix.length) {
        bestGroup = g
        bestPrefix = p
      }
    }
    if (bestGroup) {
      bestGroup.prefix = bestPrefix
      bestGroup.items.push(t)
    } else {
      groups.push({ prefix: t, items: [t] })
    }
  }
  return groups;
}


function summariseTerms(terms, minPrefixLen = 15) {
  let groups = groupTerms(terms, minPrefixLen);
  console.log("groups:");
  console.log(groups);

  const out = []
  for (const g of groups) {
    if (g.items.length > 1 && g.prefix.length >= minPrefixLen) {
      const prefix = g.prefix.replace(/[ ,]+$/, '')
      out.push(`${prefix} [${g.items.length} matches]`)
    } else {
      out.push(...g.items)
    }
  }
  return out
}


function getOntologyTermObj(term_type, term, url) {
    let termSpan = $("<span>", {class: `grid-term-link ${term_type}`, title: term, term: term, term_type: term_type});
    if (url) {
        let termLink = $("<a>", {href: url, target: "_blank", html: term});
        termSpan.append(termLink);
    } else {
        termSpan.html(term);
    }
    return termSpan;
}

function getOntologyTermLinks(term_type, term_list, getUrl) {
    let MIN_PREFIX_LENGTH = 15;
    let NUM_TERMS_TO_ATTEMPT_COLLAPSE = 5;
    let MIN_COLLAPSES = 2;
    let links = '';
    if (term_list) {
        const terms = term_list.split('|');
        const termsContainer = $("<div>", {class: "ontology-terms-container"});
        if (terms.length >= NUM_TERMS_TO_ATTEMPT_COLLAPSE) {
            let groups = groupTerms(terms, MIN_PREFIX_LENGTH);
            for (const g of groups) {
                let fullTerms = [];
                for (const term of g.items) {
                    let url = null;
                    if (getUrl) {
                        url = getUrl(term_type, term);
                    }
                    fullTerms.push(getOntologyTermObj(term_type, term, url));
                }

                if (g.items.length >= MIN_COLLAPSES && g.prefix.length >= MIN_PREFIX_LENGTH) {
                    const prefix = g.prefix.replace(/[ ,]+$/, '');
                    let collapsedLabel = `${prefix} [${g.items.length} matches]`;
                    let collapsedTerm = getOntologyTermObj(term_type, collapsedLabel);
                    collapsedTerm.addClass('collapsed-term');
                    collapsedTerm.attr("title", "click to expand full terms");

                    let group = $("<span>", {class: "term-group"});
                    let fullTermsContainer = $("<span>", {
                        class: "full-terms",
                        style: "display:none"
                    })
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
                let termObj = getOntologyTermObj(term_type, term, url);
                termObj.addClass("grid-term-link");
                termsContainer.append(termObj);
            }
        }
        links = termsContainer.prop("outerHTML");
    }
    return links
}


// This needs to be set in document.ready()
// '.ontology-terms-container .collapsed-term'
function expandCollapsedOntologyTerm() {
    const group = $(this).closest('.term-group');
    group.find('.collapsed-term').hide();
    group.find('.full-terms').show();
}
