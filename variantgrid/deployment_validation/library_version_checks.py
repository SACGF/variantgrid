from importlib import metadata

from genes.hgvs.hgvs_converter import HGVSConverterType
from snpdb.models import VariantCoordinate


def check_library_versions() -> dict:
    """ Check library versions to make sure bug fixes have been applied """

    def _test_biocommons_hgvs():
        from snpdb.models import GenomeBuild
        from genes.hgvs import HGVSMatcher
        matcher = HGVSMatcher(GenomeBuild.grch38(), hgvs_converter_type=HGVSConverterType.BIOCOMMONS_HGVS)
        # Check it can handle reference variants
        matcher.get_variant_coordinate("NC_000006.12:g.49949407_49949408=")

        # Check it can handle contig names as chrom names
        vc = VariantCoordinate(chrom="NC_000006.12", position=386486, ref="A", alt="<DUP>", svlen=5000)
        matcher.variant_coordinate_to_g_hgvs(vc)

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
            library_version_valid[name] = {
                "valid": True,
                "fix": f"All good",
            }
        except Exception as ex:
            library_version_valid[name] = {
                "valid": False,
                "fix": f"Upgrade the library using the version in requirements.txt - error {ex}",
            }

    return library_version_valid
