from django.db import models

from library.utils import invert_dict, Constant


class NucleicAcid(models.TextChoices):
    DNA = 'D', 'DNA'
    RNA = 'R', 'RNA'


class Mutation(models.TextChoices):
    GERMLINE = 'G', 'Germline'
    SOMATIC = 'S', 'Somatic'


class SimpleZygosity(models.TextChoices):
    ANY_GERMLINE = 'A', 'Het or Hom Alt'
    ANY_ZYGOSITY = 'Z', 'Any zygosity call'
    HOM_ALT = 'O', "Hom alt"
    HET = 'E', "Het"
    REF = 'R', "Ref"


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


class Sex(models.TextChoices):
    UNKNOWN = 'U', 'unknown'
    MALE = 'M', 'male'
    FEMALE = 'F', 'female'
    FILLED_IN_CHOICES = Constant((MALE, FEMALE))

    @staticmethod
    def string_to_sex(s):
        sd = invert_dict(dict(Sex.choices))
        return sd.get(s.lower(), Sex.UNKNOWN)


class GnomADPopulation(models.TextChoices):
    """ http://gnomad.broadinstitute.org/faq """
    AFRICAN_AFRICAN_AMERICAN = 'AFR', 'African/African American'
    ASHKENAZI_JEWISH = 'ASJ', 'Ashkenazi Jewish'
    EAST_ASIAN = 'EAS', 'East Asian'
    FINNISH = 'FIN', 'Finnish'
    LATINO = 'AMR', 'Latino / Mixed Amerindian'
    NON_FINNISH_EUROPEAN = 'NFE', 'Non-Finnish European'
    OTHER = 'OTH', 'Other'
    SOUTH_ASIAN = 'SAS', 'South Asian'


class PopulationGroup(models.TextChoices):
    """ Extends gnomAD for common Australian populations """
    AFRICAN_AFRICAN_AMERICAN = 'AFR', 'African/African American'
    ASHKENAZI_JEWISH = 'ASJ', 'Ashkenazi Jewish'
    AUSTRALO_MELANESIAN = 'AM', 'Australo Melanesian'
    CENTRAL_ASIAN = 'CA', 'Central Asian'
    EAST_ASIAN = 'EAS', 'East Asian'
    FINNISH = 'FIN', 'Finnish'
    LATINO = 'AMR', 'Latino / Mixed Amerindian'
    MIDDLE_EAST_NORTH_AFRICAN = 'MEA', 'Middle East / North African'
    NON_FINNISH_EUROPEAN = 'NFE', 'Non-Finnish European'
    POLYNESIAN = 'PO', 'Polynesian'
    SOUTH_ASIAN = 'SAS', 'South Asian'
    SOUTH_EAST_ASIAN = 'SEA', 'South East Asian'
    OTHER = 'OTH', 'Other'
