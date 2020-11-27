from library.utils import invert_dict


class NucleicAcid:
    DNA = 'D'
    RNA = 'R'
    CHOICES = [
        (DNA, 'DNA'),
        (RNA, 'RNA'),
    ]


class Mutation:
    GERMLINE = 'G'
    SOMATIC = 'S'
    CHOICES = [
        (GERMLINE, 'Germline'),
        (SOMATIC, 'Somatic')
    ]


class SimpleZygosity:
    ANY_GERMLINE = 'A'
    ANY_ZYGOSITY = 'Z'
    HET = 'E'
    HOM_ALT = 'O'
    REF = 'R'

    CHOICES = [
        (ANY_GERMLINE, 'Het or Hom Alt'),
        (ANY_ZYGOSITY, 'Any zygosity call'),
        (HOM_ALT, "Hom alt"),
        (HET, "Het"),
        (REF, "Ref"),
    ]


class Zygosity:
    """ Observed Variant Zygosity """
    UNKNOWN_ZYGOSITY = 'U'  # No genotype call
    HOM_REF = 'R'
    HET = 'E'
    HOM_ALT = 'O'
    # TODO: Missing shouldn't be '.' as that's confusing with './.' in VCF which is UNKNOWN_ZYGOSITY
    MISSING = '.'  # Sample has reference (ie is missing variant) in a multi-sample VCF

    CHOICES = [
        (HOM_REF, "HOM_REF"),
        (HET, "HET"),
        (HOM_ALT, "HOM_ALT"),
        (UNKNOWN_ZYGOSITY, '?'),
    ]
    REVERSE_CHOICES_LOOKUP = {v: k for k, v in CHOICES + [(MISSING, ".")]}

    GENOTYPES = {
        HOM_REF: '1/1',  # Will only be stored on ref variant
        HOM_ALT: '1/1',
        HET: '0/1',
        UNKNOWN_ZYGOSITY: './.',
        MISSING: '0/0',
    }

    ALL_ZYGOSITIES_SET = set(GENOTYPES)
    VARIANT = {HOM_REF, HET, HOM_ALT}

    @staticmethod
    def display(zygosity):
        zygosity_display = dict(Zygosity.CHOICES)
        return zygosity_display.get(zygosity, zygosity)

    @staticmethod
    def get_regex_match(zygosities_set: set, require_zygosity=True):
        if zygosities_set:
            if not require_zygosity:
                zygosities_set.update({Zygosity.UNKNOWN_ZYGOSITY})

            if zygosities_set == Zygosity.ALL_ZYGOSITIES_SET:
                return "."

            escaped_zygosities_set = set()
            for z in zygosities_set:
                if z == Zygosity.MISSING:
                    escaped_zygosities_set.add(r'\.')
                else:
                    escaped_zygosities_set.add(z)

            if len(zygosities_set) == 1:
                (zygosity_char_match,) = escaped_zygosities_set
            else:
                zygosity_char_match = '[%s]' % ''.join(escaped_zygosities_set)
        else:
            zygosity_char_match = '.'  # All
        return zygosity_char_match

    @staticmethod
    def get_genotype(zygosity):
        return Zygosity.GENOTYPES.get(zygosity, './.')

    @staticmethod
    def get_genotype_from_expanded_zygosity(expanded_zygosity):
        """ expanded_zygosity = (HOM_REF, HOM_ALT, HET) """

        return Zygosity.get_genotype(Zygosity.REVERSE_CHOICES_LOOKUP[expanded_zygosity])


class TrioSample:
    MOTHER = 'M'
    FATHER = 'F'
    PROBAND = 'P'
    CHOICES = [(MOTHER, 'Mother'),
               (FATHER, 'Father'),
               (PROBAND, 'Proband')]


class Sex:
    UNKNOWN = 'U'
    MALE = 'M'
    FEMALE = 'F'
    FILLED_IN_CHOICES = [MALE, FEMALE]
    CHOICES = [
        (UNKNOWN, 'unknown'),
        (MALE, 'male'),
        (FEMALE, 'female')
    ]

    @staticmethod
    def string_to_sex(s):
        sd = invert_dict(dict(Sex.CHOICES))
        return sd.get(s.lower(), Sex.UNKNOWN)


class GnomADPopulation:
    """ http://gnomad.broadinstitute.org/faq """
    AFRICAN_AFRICAN_AMERICAN = 'AFR'
    LATINO = 'AMR'
    ASHKENAZI_JEWISH = 'ASJ'
    EAST_ASIAN = 'EAS'
    FINNISH = 'FIN'
    NON_FINNISH_EUROPEAN = 'NFE'
    SOUTH_ASIAN = 'SAS'
    OTHER = 'OTH'

    CHOICES = [
        (AFRICAN_AFRICAN_AMERICAN, 'African/African American'),
        (ASHKENAZI_JEWISH, 'Ashkenazi Jewish'),
        (EAST_ASIAN, 'East Asian'),
        (FINNISH, 'Finnish'),
        (LATINO, 'Latino / Mixed Amerindian'),
        (NON_FINNISH_EUROPEAN, 'Non-Finnish European'),
        (OTHER, 'Other'),
        (SOUTH_ASIAN, 'South Asian'),
    ]


class PopulationGroup(GnomADPopulation):
    AUSTRALO_MELANESIAN = 'AM'
    CENTRAL_ASIAN = 'CA'
    POLYNESIAN = 'PO'
    SOUTH_EAST_ASIAN = 'SEA'
    MIDDLE_EAST_NORTH_AFRICAN = 'MEA'

    CHOICES = [
        (GnomADPopulation.AFRICAN_AFRICAN_AMERICAN, 'African/African American'),
        (GnomADPopulation.ASHKENAZI_JEWISH, 'Ashkenazi Jewish'),
        (AUSTRALO_MELANESIAN, 'Australo Melanesian'),
        (CENTRAL_ASIAN, 'Central Asian'),
        (GnomADPopulation.EAST_ASIAN, 'East Asian'),
        (GnomADPopulation.FINNISH, 'Finnish'),
        (GnomADPopulation.LATINO, 'Latino / Mixed Amerindian'),
        (MIDDLE_EAST_NORTH_AFRICAN, 'Middle East / North African'),
        (GnomADPopulation.NON_FINNISH_EUROPEAN, 'Non-Finnish European'),
        (POLYNESIAN, 'Polynesian'),
        (GnomADPopulation.SOUTH_ASIAN, 'South Asian'),
        (SOUTH_EAST_ASIAN, 'South East Asian'),
        (GnomADPopulation.OTHER, 'Other'),
    ]
