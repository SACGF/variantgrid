import gzip
import json
import logging

import requests
from django.core.management import BaseCommand
from cdot.data_release import get_latest_combo_file_urls, get_latest_data_release_tag_name, _get_version_from_tag_name

from genes.management.commands.import_gene_annotation import Command as GeneAnnotationCommand
from genes.models import TranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.constants import MINUTE_SECS
from library.utils import get_single_element, invert_dict
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        This is a wrapper around 'import_gene_annotation'

        and makes use of files produced by a spin-off project: cdot @see https://github.com/SACGF/cdot
    """
    def add_arguments(self, parser):
        parser.add_argument('--genome-build', choices=self.genome_builds, required=False)
        parser.add_argument('--annotation-consortium', choices=self.annotation_consortia, required=False)
        parser.add_argument('--force', action='store_true', help="Force even if check says we've installed latest")

    @property
    def annotation_consortia(self):
        return [ac[1] for ac in AnnotationConsortium.choices]

    @property
    def genome_builds(self):
        return [gb.name for gb in GenomeBuild.builds_with_annotation()]

    def handle(self, *args, **options):
        force = options['force']
        if genome_build := options.get('genome_build'):
            genome_builds = [genome_build]
        else:
            genome_builds = self.genome_builds
        if annotation_consortium := options.get('annotation_consortium'):
            annotation_consortia = [annotation_consortium]
        else:
            annotation_consortia = self.annotation_consortia
        logging.info("Using genome builds: %s, annotation consortia: %s", genome_builds, annotation_consortia)

        tag_name = get_latest_data_release_tag_name()
        cdot_data_version = _get_version_from_tag_name(tag_name, data_version=True)

        ac_lookup = invert_dict(dict(AnnotationConsortium.choices))
        ac_codes = [ac_lookup[ac] for ac in annotation_consortia]
        tv_qs = TranscriptVersion.objects.filter(genome_build__in=genome_builds,
                                                 transcript__annotation_consortium__in=ac_codes)
        if last_tv := tv_qs.order_by("pk").last():
            our_latest_cdot = last_tv.data.get("cdot")
        else:
            our_latest_cdot = "No cdot data in system"

        logging.info("Most recent cdot data in our database: %s", our_latest_cdot)
        logging.info("Latest cdot release on GitHub: %s", cdot_data_version)
        needs_update = cdot_data_version != our_latest_cdot
        if not needs_update:
            if force:
               logging.info("Forcing update...")
            else:
                logging.info("No need to update. Existing. Use --force to bypass this warning")
                exit(0)

        for genome_build_name in genome_builds:
            for annotation_consortium_label in annotation_consortia:
                print(f"{genome_build_name} / {annotation_consortium_label}")
                # This is somewhat wasteful as we make an API call each time below, though that's not much compared to
                # the total download and insertion time
                combo_files = get_latest_combo_file_urls(annotation_consortia={annotation_consortium_label},
                                                         genome_builds=[genome_build_name])

                genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
                annotation_consortium = AnnotationConsortium(ac_lookup[annotation_consortium_label])

                combo_file_url = get_single_element(combo_files)
                logging.info(f"%s/%s - downloading: %s",
                             genome_build, annotation_consortium.label, combo_file_url)
                response = requests.get(combo_file_url, stream=True, timeout=2 * MINUTE_SECS)
                with gzip.GzipFile(fileobj=response.raw) as fz:
                    cdot_data = json.load(fz)
                    cdot_version = cdot_data.get("cdot_version")
                    GeneAnnotationCommand.import_cdot_data(genome_build, annotation_consortium, cdot_data, cdot_version)

