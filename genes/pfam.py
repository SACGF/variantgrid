"""
PFam entries (used for annotating protein domains)
"""
import ftplib
import gzip
from io import BytesIO

import pandas as pd

from genes.models import Pfam


def store_pfam_from_web(cached_web_resource):
    """ Pfam-A.clans.tsv
      This file contains a list of all Pfam-A families that are in clans.
      The columns are: Pfam accession, clan accession, clan ID, Pfam
      ID, Pfam description. """

    PFAM_CLANS_COLUMNS = ["accession", "clan_accession", "clan_id", "pfam_id", "description"]
    ftp = ftplib.FTP("ftp.ebi.ac.uk")
    ftp.login("anonymous", "anonymous")
    buffer = BytesIO()
    ftp.retrbinary('RETR /pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz', buffer.write)
    buffer.seek(0)
    with gzip.GzipFile(fileobj=buffer) as f:
        df = pd.read_csv(f, sep='\t', names=PFAM_CLANS_COLUMNS, header=None)
        store_pfam_from_df(cached_web_resource, df)


def store_pfam_from_df(cached_web_resource, df):
    pfam_list = []
    for _, row in df.iterrows():
        pfam_list.append(Pfam(pk=Pfam.get_pk_from_accession(row["accession"]),
                              pfam_id=row["pfam_id"],
                              description=row["description"]))

    Pfam.objects.bulk_create(pfam_list)
    cached_web_resource.description = f"{len(pfam_list)} Pfam."
    cached_web_resource.save()
