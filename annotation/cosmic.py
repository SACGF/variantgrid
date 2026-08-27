# See https://clinicaltables.nlm.nih.gov/apidoc/cosmic/v4/doc.html
import requests

from snpdb.models import GenomeBuild


class CosmicAPI:
    """ The API only has GRCh37/GRCh38 data - asking for any other build returns a 500 """
    GRCH_VERSIONS = {
        "GRCh37": "37",
        "GRCh38": "38",
    }

    @classmethod
    def supports_genome_build(cls, genome_build: GenomeBuild) -> bool:
        return genome_build.name in cls.GRCH_VERSIONS

    def __init__(self, mutation_id: str, genome_build: GenomeBuild):
        grch_version = self.GRCH_VERSIONS.get(genome_build.name)
        if grch_version is None:
            raise ValueError(f"COSMIC API does not support genome build '{genome_build}'")

        fields = ["AccessionNumber", "MutationCDS"]
        params = {
            "terms": mutation_id,
            "df": ",".join(fields),
            "grchv": grch_version,
        }
        url = "https://clinicaltables.nlm.nih.gov/api/cosmic/v4/search"
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        self.data = response.json()

    def get_hgvs_list(self) -> list[str]:
        hgvs = []
        for accession, mutation_cds in self.data[3]:
            hgvs.append(f"{accession}:{mutation_cds}")
        return hgvs
