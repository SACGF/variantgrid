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

function getOntologyTermLinks(term_type, term_list) {
    let NUM_TERMS_TO_ATTEMPT_COLLAPSE = 5;
    let MIN_COLLAPSES = 200; // to be worth it
    let links = '';
    if (term_list) {
        const terms = term_list.split('|');
        const fullTermsContainer = $("<div>", {class: "container full-terms"});
        for (const term of terms) {
            let url = Urls.ontology_term_text(term_type, term);
            fullTermsContainer.append(getOntologyTermObj(term_type, term, url));
        }
        let container = fullTermsContainer;
        if (terms.length >= NUM_TERMS_TO_ATTEMPT_COLLAPSE) {
            let reducedTerms = summariseTerms(terms);
            if (reducedTerms.length < terms.length - MIN_COLLAPSES) {
                const collapsedTermsContainer = $("<div>", {class: "container collapsed-terms"});
                for (const term of reducedTerms) {
                    collapsedTermsContainer.append(getOntologyTermObj(term_type, term));
                }
                // TODO: We should make the original container hidden or something - then a link to expand
                container = collapsedTermsContainer;
            }
        }
        let numElements = container.length;
        console.log("Term elements: " + numElements);
        links = container.html();
    }
    return links;
}
