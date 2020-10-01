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
            hpo = HumanPhenotypeOntology.objects.filter(pk=idx).first()  # @UndefinedVariable
            if hpo:
                return hpo.name

        elif database == 'OMIM':
            omim = MIMMorbid.objects.filter(pk=idx).first()
            if omim:
                return omim.description

        return None

externalDatabaseRefLookupInstance = ExternalDatabaseRefLookup()
