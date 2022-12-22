import logging
from typing import Tuple, List

import pandas as pd
from dateutil import parser
from django.utils.timezone import make_aware

from genes.gene_matching import GeneSymbolMatcher
from genes.hgvs import HGVSMatcher
from genes.models import GeneSymbolWiki
from library.django_utils import UserMatcher
from library.pandas_utils import df_nan_to_none
from library.genomics.vcf_utils import write_vcf_from_tuples
from snpdb.models import ImportedWikiCollection, GenomeBuild, ImportedWiki, VariantWiki, Variant
from upload.models import UploadedWikiCollection, UploadStep
from upload.tasks.import_task import ImportTask
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from variantgrid.celery import app


def _process_imported_wiki(uploaded_file) -> Tuple[ImportedWikiCollection, List[ImportedWiki]]:
    df = pd.read_csv(uploaded_file.get_filename())
    df = df_nan_to_none(df)

    match_column_name = df.columns[0]
    if match_column_name == "Variant":
        genome_build_name = df["Genome Build"][0]  # Will all be the same
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    else:
        genome_build = None
    wiki_collection = ImportedWikiCollection.objects.create(match_column_name=match_column_name,
                                                            genome_build=genome_build)
    UploadedWikiCollection.objects.create(uploaded_file=uploaded_file,
                                          wiki_collection=wiki_collection)

    records = []
    for _, row in df.iterrows():
        created = parser.parse(row["Created"])
        created = make_aware(created)
        modified = parser.parse(row["Modified"])
        modified = make_aware(modified)
        iw = ImportedWiki(collection=wiki_collection,
                          match_column_value=row[wiki_collection.match_column_name],
                          markdown=row["Markdown"] or '',
                          username=row["User"],
                          created=created,
                          modified=modified)
        records.append(iw)
    if records:
        records = ImportedWiki.objects.bulk_create(records)

    return wiki_collection, records


def _create_or_modify_wiki(user_matcher, iw, wiki, created):
    user = user_matcher.get_user(iw.username)

    if created:
        wiki.created = iw.created
        wiki.markdown = iw.markdown
    else:
        if iw.markdown not in wiki.markdown:
            user_message = ""
            if user:
                user_message = f"by {user} on "
            modification_msg = f"-- Added from wiki import ({user_message}{iw.modified})"
            wiki.markdown += f"\n{modification_msg}\n{iw.markdown}\n"

    wiki.last_edited_by = user
    wiki.save()


class ImportGeneWikiCollection(ImportTask):
    def process_items(self, uploaded_file):
        wiki_collection, records = _process_imported_wiki(uploaded_file)

        matched_records = []
        gs_matcher = GeneSymbolMatcher()
        user_matcher = UserMatcher(uploaded_file.user)
        for iw in records:
            if gene_symbol_id := gs_matcher.get_gene_symbol_id(iw.match_column_value):
                gs_wiki, created = GeneSymbolWiki.objects.get_or_create(gene_symbol_id=gene_symbol_id)
                _create_or_modify_wiki(user_matcher, iw, gs_wiki, created)
                matched_records.append(iw)

        if matched_records:
            logging.info("Matched %d gene wiki records", len(matched_records))
            ImportedWiki.objects.bulk_update(matched_records, ["matched_wiki"], batch_size=2000)

        return len(records)


class VariantWikiCreateVCFTask(ImportVCFStepTask):
    """ Write a VCF with variants in variant wiki so they can go through normal insert pipeline """

    def process_items(self, upload_step):
        upload_pipeline = upload_step.upload_pipeline
        uploaded_file = upload_pipeline.uploaded_file

        wiki_collection, records = _process_imported_wiki(uploaded_file)
        matcher = HGVSMatcher(wiki_collection.genome_build)
        variant_tuples = []
        for iw in records:
            g_hgvs = iw.match_column_value
            if vt := matcher.get_variant_tuple(g_hgvs):
                variant_tuples.append(vt)

        if variant_tuples:
            write_vcf_from_tuples(upload_step.output_filename, variant_tuples)

        return len(records)


class VariantWikiInsertTask(ImportVCFStepTask):
    """ This is run after the VCF import data insertion stage.
        Variants will be in database at this stage """

    def process_items(self, upload_step: UploadStep):
        uploaded_file = upload_step.upload_pipeline.uploaded_file
        wiki_collection = uploaded_file.uploadedwikicollection.wiki_collection

        hgvs_matcher = HGVSMatcher(wiki_collection.genome_build)
        user_matcher = UserMatcher(uploaded_file.user)
        matched_records = []
        for iw in wiki_collection.importedwiki_set.all():
            if variant_tuple := hgvs_matcher.get_variant_tuple(iw.match_column_value):
                variant = Variant.get_from_tuple(variant_tuple, wiki_collection.genome_build)
                v_wiki, created = VariantWiki.objects.get_or_create(variant=variant)
                _create_or_modify_wiki(user_matcher, iw, v_wiki, created)
                matched_records.append(iw)

        if matched_records:
            logging.info("Matched %d variant wiki records", len(matched_records))
            ImportedWiki.objects.bulk_update(matched_records, ["matched_wiki"], batch_size=2000)

        return len(matched_records)


ImportGeneWikiCollection = app.register_task(ImportGeneWikiCollection())  # @UndefinedVariable
VariantWikiCreateVCFTask = app.register_task(VariantWikiCreateVCFTask())  # @UndefinedVariable
VariantWikiInsertTask = app.register_task(VariantWikiInsertTask())  # @UndefinedVariable
