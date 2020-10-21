from annotation.models import MonarchDiseaseOntology
from annotation.models.models_mim_hpo import HumanPhenotypeOntology, MIMMorbid


class ExternalDatabaseRefLookup:

    def lookup(self, database: str, idx: str):
        """
        Looks up information about an external database, currently will block.
        In future might call multiple times in threads.
        @param database The string representation of the database we're looking up, see regexes's "db" value for examples
        @param idx The public id of the record we're looking up, most likely.
        @return HTML summary of the reference or None if no reference could be found
        """

        if '.' not in idx:
            try:
                # strip leading 0s
                idx = str(int(idx))
            except:
                pass

        if database == 'HP':
            if hpo := HumanPhenotypeOntology.objects.filter(pk=idx).first():
                return hpo.name

        elif database == 'OMIM':
            if omim := MIMMorbid.objects.filter(pk=idx).first():
                return omim.description

        elif database == 'MONDO':
            if mondo := MonarchDiseaseOntology.objects.filter(pk=idx).first():
                return mondo.name

        return None


externalDatabaseRefLookupInstance = ExternalDatabaseRefLookup()
