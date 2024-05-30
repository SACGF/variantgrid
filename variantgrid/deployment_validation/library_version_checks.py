from importlib import metadata


def check_library_versions() -> dict:
    """ Check library versions to make sure bug fixes have been applied """

    def _test_biocommons_hgvs():
        from snpdb.models import GenomeBuild
        from genes.hgvs import HGVSMatcher
        matcher = HGVSMatcher(GenomeBuild.grch38())
        matcher.get_variant_coordinate("NC_000006.12:g.49949407_49949408=")

    minimum_versions = {
        "cdot": (0, 2, 21),
        "hgvs": _test_biocommons_hgvs,
        "pyhgvs": (0, 12, 4),
    }

    library_version_valid = {}
    for name, version_required in minimum_versions.items():
        try:
            if callable(version_required):
                version_required()
            else:
                version_str = metadata.version(name)
                version = tuple(int(i) for i in version_str.split("."))
                assert version >= version_required, "Library %s (%s) requires version >= %s" % (name, version, version_required)
            valid = True
        except:
            valid = False
        library_version_valid[name] = {
            "valid": valid,
            "fix": "Upgrade the library using the version in requirements.txt",
        }

    return library_version_valid
