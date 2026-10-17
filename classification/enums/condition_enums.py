"""
Enums for condition resolution: MultiCondition says how a classification's multiple condition terms combine.
"""
from django.db import models


class MultiCondition(models.TextChoices):
    """
    For when multiple terms exist, are we uncertain which one it is, or are we curating for multiple.
    NOT_DECIDED means a choice is still required
    """
    NOT_DECIDED = 'N', 'Not decided'
    UNCERTAIN = 'U', 'Uncertain'  # aka uncertain
    CO_OCCURRING = 'C', 'Co-occurring'  # aka combined

    @property
    def clinvar_label(self):
        return self.label
