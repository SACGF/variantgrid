# Ontology

Whenever we refer to Ontology (or an Ontology Term) we are generally doing it in the context of a human condition.
The ontology sets we mirror within the product are:
* MONDO https://mondo.monarchinitiative.org/ (condition)
* OMIM https://www.omim.org/ (condition)
* HPO https://hpo.jax.org/ (phenotype)

As MONDO seeks to be a comprehensive ontology with a hierarchy, it most cases it's our preferred model.

We also have some basic support for a few other ontology sets that we don't store locally like MedGen, DOID, etc.
Gene Symbols are also stored as an Ontology, but only to support the relationships between standard ontology sets and Gene Symbols.

The main elements of the Ontology App are:

### OntologyTerm
A single term such as MONDO:0003432. It will typically have a name, an alias, a deprecated status.
These are not versioned.

### OntologyRelation
A relationship between two OntologyTerms. It might describe a parent/child relationship in the case
of MONDO, or it might describe that a gene symbol has a known relationship to a disease.
OntologyRelation's are versioned via OntologyImport.

### OntologyImport
Describes all the relationships imported when refreshing our internal cache of relationships.