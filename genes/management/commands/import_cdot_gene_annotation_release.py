import argparse
import json

from django.core.management import BaseCommand, call_command
from django.db.models import Max

from annotation.models import VariantAnnotationVersion, AnnotationVersion, GeneAnnotationVersion, \
    InvalidAnnotationVersionError
from genes.cdot_data_release import find_latest_release_asset, download_cdot_json
from genes.gene_matching import ReleaseGeneMatcher
from genes.management.commands.import_cdot_latest import cdot_data_needs_update, import_latest_combo_file
from genes.management.commands.import_gene_annotation import Command as GeneAnnotationCommand, \
    GeneAnnotationImportManager
from genes.models import GeneAnnotationRelease, GeneVersion, TranscriptVersion, ReleaseGeneVersion, \
    ReleaseTranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.utils.file_utils import open_handle_gzip
from snpdb.models import GenomeBuild


def create_gene_annotation_release(genome_build: GenomeBuild, annotation_consortium, release_version,
                                  cdot_data: dict, cdot_version) -> GeneAnnotationRelease:
    """ A GeneAnnotationRelease doesn't change/store transcript data, but does keep track of e.g. what
        symbols are used and how things are linked together """

    # A release should be from a single GTF - so all URLs should be the same, so take any one
    random_transcript = next(iter(cdot_data["transcripts"].values()))
    url = random_transcript["genome_builds"][genome_build.name]["url"]
    ga_import_manager = GeneAnnotationImportManager()
    # For a release, this must be there as it was imported before
    import_source = ga_import_manager.get_import_source_by_url(url)
    release, created = GeneAnnotationRelease.objects.update_or_create(version=release_version,
                                                                     genome_build=genome_build,
                                                                     annotation_consortium=annotation_consortium,
                                                                     defaults={
                                                                         "gene_annotation_import": import_source
                                                                     })
    if not created:
        print("Release exists - clearing existing data")
        release.releasegeneversion_set.all().delete()
        release.releasetranscriptversion_set.all().delete()

    gene_version_ids_by_accession = GeneVersion.id_by_accession(genome_build=genome_build,
                                                               annotation_consortium=annotation_consortium)
    transcript_version_ids_by_accession = TranscriptVersion.id_by_accession(genome_build=genome_build,
                                                                           annotation_consortium=annotation_consortium)

    release_transcript_version_list = []
    gene_versions_used_by_transcripts = set()
    for transcript_accession, tv_data in cdot_data["transcripts"].items():
        # cdot has some fake transcripts (e.g. 'fake-rna-ATP6') that import_cdot_data skips
        # because they have no version — they won't be in transcript_version_ids_by_accession.
        _, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        if version is None:
            continue
        transcript_version_id = transcript_version_ids_by_accession[transcript_accession]
        rtv = ReleaseTranscriptVersion(release=release, transcript_version_id=transcript_version_id)
        release_transcript_version_list.append(rtv)

        if gene_version := tv_data.get("gene_version"):
            gene_versions_used_by_transcripts.add(gene_version)

    # It's possible that we have gene versions in this release that we don't have normally in merged files
    # as a later transcript uses a different gene so this one was never stored. We'll have to make those now
    missing_gene_versions = gene_versions_used_by_transcripts - set(gene_version_ids_by_accession)
    if missing_gene_versions:
        print(f"Missing {len(missing_gene_versions)} gene versions")
        max_gene_version = GeneVersion.objects.all().aggregate(pk__max=Max("pk"))["pk__max"]
        fake_cdot_data = {
            "genes": {k: v for k, v in cdot_data["genes"].items() if k in missing_gene_versions},
            "transcripts": {},  # Don't insert any of these
        }
        GeneAnnotationCommand.import_cdot_data(genome_build, annotation_consortium, fake_cdot_data, cdot_version)
        new_gene_versions = GeneVersion.objects.filter(genome_build=genome_build,
                                                       gene__annotation_consortium=annotation_consortium,
                                                       pk__gt=max_gene_version)
        gene_version_ids_by_accession.update({gv.accession: gv.pk for gv in new_gene_versions})

    release_gene_version_list = []
    for gene_accession in cdot_data["genes"]:
        # Gene accession may not have
        # We only store gene versions that are used in the merged files (which is what's used to insert data)
        if gene_accession in gene_versions_used_by_transcripts:
            gene_version_id = gene_version_ids_by_accession[gene_accession]
            release_gene_version_list.append(ReleaseGeneVersion(release=release, gene_version_id=gene_version_id))

    if release_gene_version_list:
        print(f"Inserting {len(release_gene_version_list)} gene versions for release {release}")
        ReleaseGeneVersion.objects.bulk_create(release_gene_version_list)

    if release_transcript_version_list:
        print(f"Inserting {len(release_transcript_version_list)} transcript versions for release {release}")
        ReleaseTranscriptVersion.objects.bulk_create(release_transcript_version_list)

    print("Matching existing gene list symbols to this release...")
    gm = ReleaseGeneMatcher(release)
    gm.match_unmatched_in_hgnc_and_gene_lists()
    return release


def read_cdot_data(file_obj) -> tuple[dict, str]:
    """ A release needs random access / multiple passes over the data, so load it all in """
    cdot_data = json.load(file_obj)
    cdot_version = cdot_data.get("cdot_version")
    if cdot_version is None:
        raise ValueError("JSON does not contain 'cdot_version' key")
    print(f"JSON uses cdot version {cdot_version}")
    return cdot_data, cdot_version


class Command(BaseCommand):
    """
        Install the GeneAnnotationRelease matching the VEP a build was annotated with.

        A release has to hold the exact gene/transcript versions VEP used, so we resolve which cdot
        per-GFF file corresponds to the VariantAnnotationVersion's VEP version strings, download it
        from the cdot GitHub data release, import it, and link it to the VariantAnnotationVersion.

        @see https://github.com/SACGF/cdot
    """

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', choices=self.genome_builds, required=False,
                            help="Defaults to every build with annotation")
        parser.add_argument('--status', choices=list(VariantAnnotationVersion.Status), required=False,
                            help="Only act on VariantAnnotationVersions with this status. Defaults to every "
                                 "NEW and ACTIVE version, so a half configured build gets finished off")
        parser.add_argument('--json-file', required=False,
                            help="Create the release from a locally produced cdot JSON.gz instead of "
                                 "downloading the published one. Acts on the latest version only, so it "
                                 "requires --genome-build")
        parser.add_argument('--release-name', required=False,
                            help="Name for the GeneAnnotationRelease. Defaults to one derived from the "
                                 "VariantAnnotationVersion's VEP versions")
        parser.add_argument('--relink', action='store_true',
                            help="Re-point versions whose linked release disagrees with their VEP versions. "
                                 "This changes which transcript versions existing analyses resolve to")
        parser.add_argument('--force', action='store_true',
                            help="Re-download and rebuild the release even if one is already linked")
        parser.add_argument('--gene-annotation', action=argparse.BooleanOptionalAction, default=True,
                            help="Create any missing GeneAnnotationVersion and bring the AnnotationVersion "
                                 "into line afterwards. Default on")

    @property
    def genome_builds(self):
        return [gb.name for gb in GenomeBuild.builds_with_annotation()]

    def handle(self, *args, **options):
        if build_name := options.get("genome_build"):
            genome_builds = [GenomeBuild.get_name_or_alias(build_name)]
        elif options.get("json_file"):
            raise ValueError("--json-file requires --genome-build")
        else:
            genome_builds = list(GenomeBuild.builds_with_annotation())

        releases_by_build = {}
        for genome_build in genome_builds:
            print(f"=== {genome_build} ===")
            releases = {}  # keyed by pk so several versions sharing a release only annotate once
            for vav in self._variant_annotation_versions(genome_build, options):
                if release := self._ensure_release(vav, options):
                    releases[release.pk] = release
            if releases:
                releases_by_build[genome_build] = list(releases.values())

        if not releases_by_build:
            return

        if options["gene_annotation"]:
            for genome_build, releases in releases_by_build.items():
                self._ensure_gene_annotation(genome_build, releases)
            return

        missing = [r for releases in releases_by_build.values() for r in releases
                   if not GeneAnnotationVersion.objects.filter(gene_annotation_release=r).exists()]
        if missing:
            release_ids = ",".join(str(r.pk) for r in missing)
            print(f"These releases have no GeneAnnotationVersion: {', '.join(str(r) for r in missing)}. "
                  f"To create them, run:")
            print(f"    python3 manage.py gene_annotation --gene-annotation-release={release_ids}")
        for genome_build in releases_by_build:
            self._report_annotation_versions(genome_build)

    def _variant_annotation_versions(self, genome_build: GenomeBuild, options) -> list[VariantAnnotationVersion]:
        if options.get("json_file"):
            # A single local file makes a single release, so only the newest version can be meant
            vav = self._latest_variant_annotation_version(genome_build, options.get("status"))
            return [vav] if vav else []

        if status := options.get("status"):
            statuses = [status]
        else:
            # NEW is the version being built, ACTIVE the one in use - both need a release
            statuses = [VariantAnnotationVersion.Status.NEW, VariantAnnotationVersion.Status.ACTIVE]

        vavs = list(VariantAnnotationVersion.objects.filter(genome_build=genome_build, status__in=statuses)
                    .order_by("annotation_date"))
        if not vavs:
            print(f"No {'/'.join(statuses)} VariantAnnotationVersion - run "
                  f"'python3 manage.py create_new_variant_annotation_version --genome-build={genome_build}'")
        return vavs

    @staticmethod
    def _latest_variant_annotation_version(genome_build: GenomeBuild, status) -> VariantAnnotationVersion:
        if status:
            return VariantAnnotationVersion.latest(genome_build, status=status)
        # A release is normally installed for the version being built, but a build whose version has
        # already been promoted still needs one
        for fallback_status in (VariantAnnotationVersion.Status.NEW, VariantAnnotationVersion.Status.ACTIVE):
            if vav := VariantAnnotationVersion.latest(genome_build, status=fallback_status):
                return vav
        return None

    def _ensure_release(self, vav: VariantAnnotationVersion, options) -> GeneAnnotationRelease:
        """ Returns the release this version should be using, whether we just installed it or it was
            already linked - the caller still has to check it has a GeneAnnotationVersion """
        print(f"{vav.status} {vav}")
        if vav.gene_annotation_release and not options["force"]:
            return self._check_existing_release(vav, options)

        annotation_consortium = AnnotationConsortium(vav.annotation_consortium)
        if json_file := options.get("json_file"):
            release_name = self._release_name(vav, options)
            if release_name is None:
                return None
            with open_handle_gzip(json_file, "rb") as f:
                cdot_data, cdot_version = read_cdot_data(f)
        else:
            # An existing release often already matches - a build's annotation doesn't always change
            # between VEP versions
            if not options["force"]:
                if release := vav.link_gene_annotation_release():
                    print(f"  linked existing GeneAnnotationRelease '{release}' - no download needed")
                    return release

            release_token = vav.cdot_gene_release_token
            if release_token is None:
                print(f"  could not work out the gene set release from this VEP. "
                      f"refseq={vav.refseq!r} vep={vav.vep} vep_cache={vav.vep_cache} genebuild={vav.genebuild!r}")
                return None

            release_name = self._release_name(vav, options)
            if release_name is None:
                return None

            print(f"  VEP used {annotation_consortium.label} release '{release_token}'")
            asset, available = find_latest_release_asset(vav.genome_build.name, annotation_consortium.label,
                                                        release_token)
            if asset is None:
                available_releases = ", ".join(a.release for a in available) or "none"
                print(f"  the latest cdot data release has no file for "
                      f"{vav.genome_build}/{annotation_consortium.label} release '{release_token}'. "
                      f"It publishes: {available_releases}")
                return None

            self._ensure_transcripts_installed(vav.genome_build, annotation_consortium)
            print(f"  creating GeneAnnotationRelease '{release_name}' from {asset.filename}")
            with download_cdot_json(asset.url) as f:
                cdot_data, cdot_version = read_cdot_data(f)

        release = create_gene_annotation_release(vav.genome_build, annotation_consortium, release_name,
                                                 cdot_data, cdot_version)
        vav.gene_annotation_release = release
        vav.save()
        print(f"  linked GeneAnnotationRelease '{release}'")
        return release

    @staticmethod
    def _check_existing_release(vav: VariantAnnotationVersion, options) -> GeneAnnotationRelease:
        """ Already linked - keep it, unless it disagrees with the VEP versions and --relink was passed """
        release = vav.gene_annotation_release
        match = vav.match_gene_annotation_release()
        if match.release_ids and release.pk not in match.release_ids:
            matched = GeneAnnotationRelease.objects.filter(pk__in=match.release_ids)
            matched_str = ", ".join(str(r) for r in matched)
            if options["relink"] and (release_id := match.unique_release_id):
                vav.gene_annotation_release = GeneAnnotationRelease.objects.get(pk=release_id)
                vav.save()
                print(f"  re-pointed from '{release}' to '{vav.gene_annotation_release}' "
                      f"(matched on {match.description})")
                return vav.gene_annotation_release
            print(f"  linked to '{release}' but its VEP versions match {matched_str} "
                  f"(on {match.description}) - pass --relink to re-point it")
            return release

        print(f"  has GeneAnnotationRelease '{release}'")
        return release

    @staticmethod
    def _release_name(vav: VariantAnnotationVersion, options) -> str:
        release_name = options.get("release_name") or vav.suggested_gene_annotation_release_name
        if release_name == VariantAnnotationVersion.UNKNOWN_GENE_ANNOTATION_RELEASE_NAME:
            print(f"  could not derive a release name from this VEP (refseq={vav.refseq!r} "
                  f"vep={vav.vep} genebuild={vav.genebuild!r}) - pass --release-name")
            return None
        return release_name

    def _ensure_gene_annotation(self, genome_build: GenomeBuild, releases: list[GeneAnnotationRelease]):
        """ A release is only usable once it has a GeneAnnotationVersion and the AnnotationVersion points
            at it - otherwise AnnotationVersion.validate() raises and variant pages break """
        for release in releases:
            if GeneAnnotationVersion.objects.filter(gene_annotation_release=release).exists():
                continue
            print(f"Creating GeneAnnotationVersion for {release}")
            call_command("gene_annotation", gene_annotation_release=release.pk)

        print(f"Bringing {genome_build} AnnotationVersion into line")
        AnnotationVersion.new_sub_version(genome_build)
        self._report_annotation_versions(genome_build)

    @staticmethod
    def _report_annotation_versions(genome_build: GenomeBuild):
        for status in (VariantAnnotationVersion.Status.NEW, VariantAnnotationVersion.Status.ACTIVE):
            av = AnnotationVersion.latest(genome_build, validate=False, status=status)
            if av is None:
                continue
            try:
                av.validate()
                print(f"  {status} AnnotationVersion {av.pk}: ok")
            except InvalidAnnotationVersionError as e:
                print(f"  {status} AnnotationVersion {av.pk}: {e}")

    @staticmethod
    def _ensure_transcripts_installed(genome_build: GenomeBuild, annotation_consortium: AnnotationConsortium):
        """ A release only links to existing gene/transcript versions, so the combo file (which holds the
            latest copy of every transcript, including those in the per-GFF files) has to be in first """
        needs_update, ours, latest = cdot_data_needs_update([genome_build.name], [annotation_consortium.label])
        if needs_update:
            print(f"Our transcripts were built by cdot {ours}, latest is {latest} - importing combo file")
            import_latest_combo_file(genome_build, annotation_consortium)
        else:
            print(f"Transcripts already at latest cdot ({latest})")
