from django.db.models.query_utils import Q

from annotation.models.models_enums import VariantClass


class AbstractPathogenicity:
    """ CHOICES *must* be in ascending order of pathogenicity! """

    MINIMUM_FLAG_DAMAGE_LEVEL = None  # 1st choice that causes damage
    VARIANT_PATH = None  # QuerySet path to field - ie "variantannotation__impact"
    CHOICES = []

    @classmethod
    def get_most_damaging(cls, predictions_list):
        most_damaging = None
        predictions_set = set(predictions_list)
        for prediction, _ in cls.CHOICES:
            if prediction in predictions_set:
                most_damaging = prediction
        return most_damaging

    @classmethod
    def get_damage_or_greater_levels(cls, min_level=None):
        if min_level is None:
            min_level = cls.MINIMUM_FLAG_DAMAGE_LEVEL

        higher_levels = set()
        found = False
        for k, _ in cls.CHOICES:
            if k == min_level:
                found = True
            if found:
                higher_levels.add(k)
        return higher_levels

    @classmethod
    def _get_q_for_damage_levels(cls, damage_levels, allow_null=False) -> Q:
        q = Q(**{cls.VARIANT_PATH + "__in": damage_levels})
        if allow_null:
            q |= Q(**{cls.VARIANT_PATH + "__isnull": True})
        return q

    @classmethod
    def get_q(cls, min_level=None, allow_null=False) -> Q:
        damage_levels = cls.get_damage_or_greater_levels(min_level)
        return cls._get_q_for_damage_levels(damage_levels, allow_null)

    @classmethod
    def is_level_flagged(cls, level):
        """ level >= MINIMUM_FLAG_DAMAGE_LEVEL """
        if level is None:
            return False
        level_score = {l: i for i, l in enumerate(dict(cls.CHOICES))}
        return level_score[level] >= level_score[cls.MINIMUM_FLAG_DAMAGE_LEVEL]


class PathogenicityImpact(AbstractPathogenicity):
    MODIFIER = 'O'
    LOW = 'L'
    MODERATE = 'M'
    MODERATE_OTHER = '*'  # This is a custom filter not a VEP value, so won't be in annotation data
    HIGH = 'H'

    CHOICES = [
        (MODIFIER, "MODIFIER"),
        (LOW, "LOW"),
        (MODERATE, "MODERATE"),
        (MODERATE_OTHER, "MODERATE*"),
        (HIGH, 'HIGH'),
    ]
    MINIMUM_FLAG_DAMAGE_LEVEL = MODERATE
    VARIANT_PATH = "variantannotation__impact"

    @classmethod
    def get_damage_or_greater_levels(cls, min_level=None):
        levels = super().get_damage_or_greater_levels(min_level)
        levels.discard(cls.MODERATE_OTHER)
        return levels

    @classmethod
    def get_q(cls, min_level=None, allow_null=False):
        q = super().get_q(min_level=min_level, allow_null=allow_null)
        if min_level == cls.MODERATE_OTHER:
            # SNV=HIGH, non SNV HIGH or MODERATE
            q_non_snv = cls._get_q_for_damage_levels([cls.MODERATE], allow_null)
            q |= q_non_snv & ~Q(variantannotation__variant_class=VariantClass.SNV)
        return q


class SIFTPrediction(AbstractPathogenicity):
    TOLERATED = 'T'
    DAMAGING = 'D'

    CHOICES = [
        (TOLERATED, "Tolerated"),
        (DAMAGING, "Damaging"),
    ]
    MINIMUM_FLAG_DAMAGE_LEVEL = DAMAGING
    VARIANT_PATH = "variantannotation__sift"


class Polyphen2Prediction(AbstractPathogenicity):
    BENIGN = 'B'
    POSSIBLY_DAMAGING = 'P'
    PROBABLY_DAMAGING = 'D'

    CHOICES = [
        (BENIGN, "Benign"),
        (POSSIBLY_DAMAGING, "Possibly Damaging"),
        (PROBABLY_DAMAGING, "Probably Damaging"),
    ]
    MINIMUM_FLAG_DAMAGE_LEVEL = POSSIBLY_DAMAGING
    VARIANT_PATH = "variantannotation__polyphen2_hvar_pred_most_damaging"


class MutationTasterPrediction(AbstractPathogenicity):
    DISEASE_CAUSING_AUTOMATIC = 'A'
    DISEASE_CAUSING = 'D'
    POLYMORPHISM = 'N'
    POLYMORPHISM_AUTOMATIC = 'P'

    CHOICES = [  # Ordered least to most pathogenic
        (POLYMORPHISM_AUTOMATIC, "Polymorphism (automatic)"),
        (POLYMORPHISM, "Polymorphism"),
        (DISEASE_CAUSING, "Disease Causing"),
        (DISEASE_CAUSING_AUTOMATIC, "Disease causing (automatic)"),
    ]
    MINIMUM_FLAG_DAMAGE_LEVEL = DISEASE_CAUSING
    VARIANT_PATH = "variantannotation__mutation_taster_pred_most_damaging"


class FATHMMPrediction(AbstractPathogenicity):
    TOLERATED = 'T'
    DAMAGING = 'D'

    CHOICES = [
        (TOLERATED, "Tolerated"),
        (DAMAGING, "Damaging"),
    ]
    MINIMUM_FLAG_DAMAGE_LEVEL = DAMAGING
    VARIANT_PATH = "variantannotation__fathmm_pred_most_damaging"


class MutationAssessorPrediction(AbstractPathogenicity):
    NEUTRAL = 'N'
    LOW = 'L'
    MEDIUM = 'M'
    HIGH = 'H'

    CHOICES = [
        (NEUTRAL, 'Neutral'),
        (LOW, 'Low'),
        (MEDIUM, 'Medium'),
        (HIGH, 'High'),
    ]
    MINIMUM_FLAG_DAMAGE_LEVEL = MEDIUM
    VARIANT_PATH = "variantannotation__mutation_assessor_pred_most_damaging"
