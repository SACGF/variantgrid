from annotation.models.models_mim_hpo import HumanPhenotypeOntology, HPOSynonym, MIMMorbid, MIMMorbidAlias


def create_mim_hpo_test_data():
    MIM_AND_ALIAS = [{'description': 'MULTIPLE ENDOCRINE NEOPLASIA, TYPE I; MEN1', 'accession': 131100},
                     {'description': 'COLORECTAL CANCER; CRC', 'accession': 114500},
                     {'description': 'MAPLE SYRUP URINE DISEASE; MSUD', 'accession': 248600},
                     {'description': 'MEGALOBLASTIC ANEMIA 1', 'accession': 261100},
                     {'description': 'PLATELET DISORDER, FAMILIAL, WITH ASSOCIATED MYELOID MALIGNANCY; FPDMM',
                      "accession": 601399},
                     {"description": 'LEUKEMIA, ACUTE MYELOID', "accession": 79929}]

    MIM_ALIAS = [{'id': 53, 'mim_morbid_id': 261100, 'description': 'IMERSLUND-GRASBECK SYNDROME; IGS'}]

    HPO = [{'id': 1657, 'name': 'Prolonged QT interval'},
           {'id': 1508, 'name': 'Failure to thrive'},
           {'id': 2925, 'name': 'Thyroid-stimulating hormone excess'},
           {'id': 4762, 'name': 'Hypoplasia of right ventricle'},
           {'id': 200063, 'name': 'Colorectal polyposis'}]

    HPO_SYNONYMS = [{'id': 3224, 'scope': 'E', 'name': 'Long QT syndrome', 'hpo_id': 1657},
                    {'id': 15668, 'scope': 'E', 'name': 'Colorectal polyps', 'hpo_id': 200063}]

    for mim in MIM_AND_ALIAS:
        mim_morbid = MIMMorbid.objects.create(**mim)
        # Because we auto-complete MIMMorbidAlias (to get all the terms), we need to create a duplicate of MIMMorbid as an alias
        # MIMMorbid.pk == MIMMorbidAlias.pk for easy pk checking in check_expected_results_for_description()
        MIMMorbidAlias.objects.create(mim_morbid=mim_morbid, id=mim["accession"], description=mim["description"])

    for mim_alias in MIM_ALIAS:
        MIMMorbidAlias.objects.create(**mim_alias)

    for hpo in HPO:
        HumanPhenotypeOntology.objects.create(**hpo)  # @UndefinedVariable

    for hpo_synonym in HPO_SYNONYMS:
        HPOSynonym.objects.create(**hpo_synonym)
