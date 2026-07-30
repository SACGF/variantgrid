import logging

from cdot.data_release import get_latest_combo_file_urls, get_latest_data_release_tag_name, _get_version_from_tag_name
from django.core.management import BaseCommand

from genes.cdot_data_release import download_cdot_json
from genes.management.commands.import_gene_annotation import Command as GeneAnnotationCommand
from genes.models import TranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.utils import get_single_element, invert_dict
from snpdb.models import GenomeBuild


def get_installed_cdot_data_version(genome_builds, annotation_consortia) -> str:
    """ The cdot version that produced the most recently inserted transcripts """
    ac_lookup = invert_dict(dict(AnnotationConsortium.choices))
    ac_codes = [ac_lookup[ac] for ac in annotation_consortia]
    tv_qs = TranscriptVersion.objects.filter(genome_build__in=genome_builds,
                                             transcript__annotation_consortium__in=ac_codes)
    if last_tv := tv_qs.order_by("pk").last():
        return last_tv.data.get("cdot")
    return "No cdot data in system"


def get_latest_cdot_data_version() -> str:
    tag_name = get_latest_data_release_tag_name()
    return _get_version_from_tag_name(tag_name, data_version=True)


def cdot_data_needs_update(genome_builds, annotation_consortia) -> tuple[bool, str, str]:
    """ Returns (needs update, ours, latest on GitHub) """
    our_latest_cdot = get_installed_cdot_data_version(genome_builds, annotation_consortia)
    cdot_data_version = get_latest_cdot_data_version()
    logging.info("Most recent cdot data in our database: %s", our_latest_cdot)
    logging.info("Latest cdot release on GitHub: %s", cdot_data_version)
    return cdot_data_version != our_latest_cdot, our_latest_cdot, cdot_data_version


def import_latest_combo_file(genome_build: GenomeBuild, annotation_consortium: AnnotationConsortium):
    """ Import the combo file - the latest copy of every transcript for a build, which is a superset
        of the per-GFF release files """
    combo_files = get_latest_combo_file_urls(annotation_consortia={annotation_consortium.label},
                                             genome_builds=[genome_build.name])
    combo_file_url = get_single_element(combo_files)
    logging.info("%s/%s - importing combo file", genome_build, annotation_consortium.label)
    with download_cdot_json(combo_file_url) as f:
        cdot_version = GeneAnnotationCommand.read_cdot_version(f)
        GeneAnnotationCommand.import_cdot_data_file(genome_build, annotation_consortium, f, cdot_version)


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

        needs_update, _ours, _latest = cdot_data_needs_update(genome_builds, annotation_consortia)
        if not needs_update:
            if force:
                logging.info("Forcing update...")
            else:
                logging.info("No need to update. Exiting. Use --force to bypass this warning")
                return

        ac_lookup = invert_dict(dict(AnnotationConsortium.choices))
        for genome_build_name in genome_builds:
            for annotation_consortium_label in annotation_consortia:
                print(f"{genome_build_name} / {annotation_consortium_label}")
                genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
                annotation_consortium = AnnotationConsortium(ac_lookup[annotation_consortium_label])
                # This is somewhat wasteful as we make an API call each time, though that's not much
                # compared to the total download and insertion time
                import_latest_combo_file(genome_build, annotation_consortium)
