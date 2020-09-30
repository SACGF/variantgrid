from django.contrib import admin

from pathtests import models

#admin.site.register(models.Case)
#admin.site.register(models.CaseClinician)
#admin.site.register(models.FollowLeadScientist)
admin.site.register(models.PathologyTest)
admin.site.register(models.PathologyTestVersion)
admin.site.register(models.PathologyTestSynonyms)
admin.site.register(models.PathologyTestGeneModificationRequest)
