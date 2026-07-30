"""
Canonical form for GeneAnnotationImport.url

The URL is how GeneAnnotationImportManager decides whether it has already seen a GTF/GFF,
so the same underlying annotation recorded under a different scheme or file format
produces a second GeneAnnotationImport row. Reducing URLs to a canonical form lets those
be recognised as the same import.
"""


def canonical_import_url(url: str) -> str:
    """ https, and Ensembl GTF rather than GFF3 - the form cdot records """
    scheme, separator, rest = url.partition("://")
    # UTA imports are a database connection, not a file - leave them exactly as they are
    if not separator or scheme == "postgresql":
        return url

    # Ensembl publish the same release as both GFF3 and GTF in sibling directories, and cdot reads
    # the GTF. Keyed on the directory rather than the extension because RefSeq archive files are
    # also named '.gff3.gz' but have no GTF equivalent.
    if "/gff3/" in rest:
        rest = rest.replace("/gff3/", "/gtf/").replace(".gff3.gz", ".gtf.gz")
    return "https://" + rest
