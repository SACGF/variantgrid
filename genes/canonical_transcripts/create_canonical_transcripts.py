import logging

from genes.gene_matching import GeneSymbolMatcher
from genes.models import CanonicalTranscriptCollection, CanonicalTranscript, Transcript
from library.database_utils import sql_delete_qs
from library.file_utils import file_md5sum
import pandas as pd

from snpdb.models import GenomeBuild

DELETE_IN_PARTS = False


def delete_old_canonical_transcript_collection_data(filename):
    """ Can only be 1 per kit so need to delete it """

    try:
        # This uses a lot of memory, so delete big tables in SQL
        ctc = CanonicalTranscriptCollection.objects.get(filename=filename)
        if DELETE_IN_PARTS:
            qs = ctc.genecoveragecanonicaltranscript_set.all()
            logging.info("Deleting gene coverage canonical transcripts")
            sql_delete_qs(qs)

            logging.info("Deleting canonical transcripts")
            qs = ctc.canonicaltranscript_set.all()
            sql_delete_qs(qs)

        ctc.delete()
    except CanonicalTranscriptCollection.DoesNotExist:
        pass
    except MemoryError:
        logging.error("Out of memory!")
        logging.error("This table is really large, perhaps it's easier to delete everything and start again.")
        logging.error("via running SQL: truncate genes_genecoveragecanonicaltranscript;")
        raise


def create_canonical_transcript_collection(genome_build: GenomeBuild, annotation_consortium, filename, gene_matcher=None):
    GENE_EXON_COLUMN = "Gene:Transcript"
    delete_old_canonical_transcript_collection_data(filename)

    file_md5hash = file_md5sum(filename)
    collection = CanonicalTranscriptCollection.objects.create(filename=filename,
                                                              description="",
                                                              genome_build=genome_build,
                                                              annotation_consortium=annotation_consortium,
                                                              file_md5sum=file_md5hash)

    if gene_matcher is None:
        gene_matcher = GeneSymbolMatcher()

    df = pd.read_csv(filename, index_col=None, sep='\t')
    gene_exon_series = df[GENE_EXON_COLUMN].str.split(":")
    del df

    known_transcript_ids = Transcript.known_transcript_ids(genome_build, annotation_consortium)
    canonical_transcript_list = []
    for (original_gene_symbol, original_transcript_id) in gene_exon_series:
        gene_symbol_id = gene_matcher.get_gene_symbol_id(original_gene_symbol)
        if original_transcript_id in known_transcript_ids:
            transcript = original_transcript_id
        else:
            transcript = None
        ct = CanonicalTranscript(collection=collection,
                                 gene_symbol_id=gene_symbol_id,
                                 transcript=transcript,
                                 original_gene_symbol=original_gene_symbol,
                                 original_transcript_id=original_transcript_id)
        canonical_transcript_list.append(ct)

    if canonical_transcript_list:
        CanonicalTranscript.objects.bulk_create(canonical_transcript_list)

    return collection
