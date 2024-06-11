# See https://clinicaltables.nlm.nih.gov/apidoc/cosmic/v4/doc.html
import requests

from snpdb.models import GenomeBuild


class CosmicAPI:
    def __init__(self, mutation_id: str, genome_build: GenomeBuild):
        fields = ["AccessionNumber", "MutationCDS"]
        build_numer = genome_build.name.lstrip("GRCh")
        params = {
            "terms": mutation_id,
            "df": ",".join(fields),
            "grchv": build_numer,
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
