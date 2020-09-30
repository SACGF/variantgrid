from io import StringIO
from datetime import datetime
from operator import attrgetter
from typing import Optional, List
from xml.etree.cElementTree import Element, SubElement, tostring

from django.contrib.sites.models import Site
from django.utils.timezone import now

from library.utils import group_by_key, cautious_attempt_html_to_text
from snpdb.models import Lab
from variantclassification.enums import SpecialEKeys
from variantclassification.models import EvidenceKeyMap, VCBlobKeys
from variantclassification.models.variant_classification import VariantClassificationModification
from variantclassification.regexes import db_ref_regexes
from variantclassification.views.clinvar.clinvar import SubmissionSetType, SubmitterType, PersonType, NameType, \
    OrganizationType, ContactType, OrganizationCategoryList, InstitutionType, SubmissionType, RecordStatusType, \
    ReleaseStatusType, ClinvarSubmissionIDType, MeasureTraitType, AssertionType, AssertionTypeType, AssertionType4, \
    ClinicalSignificanceType, ReviewStatusType, ObservationSet, SampleType, OriginType, SpeciesType, GenderType, \
    MeasureSetType, MeasureSetTypeType, MeasureType, MeasureTypeType, AttributeSetType, MeasureAttributeTypeType, \
    MeasureRelationshipAttributetype, MeasureRelationshipType, MeasureRelationshipTypeType, MeasureRelationshiptype, \
    ElementValueType, ElementValueTypeType, ElementValuetype, Measuresettypelist, Measuretype, SetElementSetType, \
    TraitSetType, TraitSetTypeType, TraitRelationshiptype, TraitTypeType, TraitType, AttributeType, MethodType, \
    MethodTypeType, Methodtypelist, ObservedDataType, MeasureAttributetype, XrefType, CitationType, IDType, \
    ZygosityType, ObsAttributeTypeType, ObsAttributetype
from variantclassification.views.variant_classification_export_utils import ExportFormatter


class LabGrouping:
    def __init__(self, lab: Lab, records: List[VariantClassificationModification]):
        self.lab = lab
        self.records = records

CLINVAR_CS_DICT = {
    "B": "Benign",
    "LB": "Likely Bengign",
    "VUS": "Uncertain significance",
    "LP": "Likely Pathogenic",
    "P": "Pathogenic"
}

class ExportFormatterClinvar(ExportFormatter):

    @property
    def version(self):
        return '0.1'

    def __init__(self, filename_override:str = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filename_override = filename_override

    def row_iterator(self):
        by_lab_qs = self.qs.order_by('variant_classification__lab_id')
        for lab, records in group_by_key(by_lab_qs, attrgetter('variant_classification.lab')):
            yield LabGrouping(lab, records)

    def header(self) -> Optional[str]:
        return None

    def row(self, vcm: LabGrouping) -> str:
        """
        Clinvar Export
        """
        lab = vcm.lab
        records = vcm.records

        # TODO should this be last modified of the record
        export_date = datetime.today().strftime('%Y-%m-%d')
        current_site: Site = Site.objects.get_current()

        submitter_of_record = SubmitterType(
            Person=PersonType(
                Name=NameType(
                    First=self.user.first_name,
                    Last=self.user.last_name,
                    NCBIOrganizationID=555
                )
            ),
            # 'lab', 'LSDB', 'clinic', 'resource', 'consortium', 'patient registry', 'other'

            Organization=OrganizationType(
                Name=current_site.name,
                OrganizationCategory=OrganizationCategoryList.RESOURCE.value,
                LabContact=[ContactType(
                    Email=self.user.email
                )],
                URL=current_site.domain
            ),
        )
        # use the lab head(s) as contacts for the submitting lab
        lab_head_contacts = [ContactType(Email=lh.user.email) for lh in lab.labhead_set.all()]

        submitter = SubmitterType(
            Organization=OrganizationType(
                Name=lab.name,
                Institution=InstitutionType(valueOf_=lab.organization.name),
                OrganizationCategory=OrganizationCategoryList.LAB.value,
                LabContact=lab_head_contacts
            )
        )

        submissions = []
        for vcm in records:
            vc = vcm.variant_classification
            cs = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            cs_label = CLINVAR_CS_DICT.get(cs)
            if not cs_label:
                cs_label = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(cs)

            pubmed_citations = [CitationType(
                ID=IDType(Source='PubMed', valueOf_=ref.get('idx'))
            ) for ref in vcm.db_refs if ref.get('db') == 'PubMed']

            variant_x_refs = None
            if vcm.get(SpecialEKeys.DB_SNP):
                variant_x_refs = [XrefType(
                    db="dbSNP", id=vcm.get(SpecialEKeys.DB_SNP), type_="rsNumber"
                )]
            else:
                #  too dangerous just to find any and all rsXXXXXX in the classification
                #  but only SA Pathology uses the dbSNP ID field
                # variant_x_refs = [XrefType(
                #    db="dbSNP", id=ref.get('idx'), type_="rsNumber"
                # ) for ref in vcm.db_refs if ref.get('db') == 'SNP']
                pass

            # TODO If several OMIMs are listed, should they be considered alternatives to each other,
            # or multiple traits, should we say "Alternate" if we see "uncertain" and otherwise
            # have them as multile traits?
            trait = TraitType(
                TraitType=TraitTypeType(val_type="name", valueOf_=TraitRelationshiptype.DISEASE),
            )

            zygosity = None
            zygosity_val = vcm.get(SpecialEKeys.ZYGOSITY)
            if isinstance(zygosity_val, list):
                zygosity_val = zygosity_val[0]
            zygosity_obs = None
            if zygosity_val:
                if zygosity_val == 'homozygous':
                    zygosity = ZygosityType.HOMOZYGOTE.value
                elif zygosity_val == 'compound_heterozygous':
                    zygosity = ZygosityType.COMPOUNDHETEROZYGOTE.value
                elif zygosity_val == 'hemizygous':
                    zygosity = ZygosityType.HEMIZYGOTE.value
                elif zygosity_val == 'heteroplasmic':
                    zygosity = ZygosityType.SINGLEHETEROZYGOTE.value
                elif zygosity == 'mosaic':
                    zygosity_obs = ObsAttributeTypeType(
                        val_type=ObsAttributetype.NUMBER_MOSAIC.name(),
                        valueOf_=1
                    )
                    pass

            condition_xrefs = vcm.evidence.get(SpecialEKeys.CONDITION, {}).get(VCBlobKeys.DB_REFS.value)
            if condition_xrefs:
                val_type = ElementValuetype.PREFERRED.value
                for condition_xref in condition_xrefs:
                    db, idx, url, summary = condition_xref.get('db'), condition_xref.get('idx'), condition_xref.get('url'), condition_xref.get('summary')
                    xref_type = 'MIM' if db == 'OMIM' else None
                    if summary:
                        name = summary
                        symbol = None
                        summary_parts = summary.split(';')
                        if len(summary_parts) == 2:
                            name = summary_parts[0].strip()
                            symbol = summary_parts[1].strip()
                        trait.get_Name().append(
                            SetElementSetType(
                                ElementValue=ElementValueTypeType(val_type=val_type, valueOf_=name),
                                XRef=[XrefType(db=db, id=idx, URL=url, type_=xref_type)]
                            )
                        )
                        if symbol:
                            trait.get_Symbol().append(
                                SetElementSetType(
                                    ElementValue=ElementValueTypeType(val_type=val_type, valueOf_=symbol),
                                    XRef=[XrefType(db=db, id=idx, URL=url, type_=xref_type)]
                                )
                            )
                        val_type = ElementValuetype.ALTERNATE.value

            if not trait.get_Name():
                trait.get_Name().append(SetElementSetType(
                    ElementValue=ElementValueTypeType(val_type=ElementValuetype.PREFERRED.value,
                                                      valueOf_=vcm.get(SpecialEKeys.CONDITION))
                ))

            submitter_date = vcm.get(SpecialEKeys.CURATION_DATE) or vcm.created.strftime('%Y-%m-%d')  # TODO check CURATION_DATE is valid
            current_date = now().strftime('%Y-%m-%d')
            submission_text = '\n'.join([text for text in [
                vcm.get(SpecialEKeys.INTERPRETATION_SUMMARY),
                vcm.criteria_strength_summary(only_acmg=True)
            ] if text])
            submission = SubmissionType(
                RecordStatus=RecordStatusType.NOVEL.value,  # TODO should this be DELETE/UPDATE on subsequent runs
                ReleaseStatus=ReleaseStatusType.PUBLIC.value,
                ClinvarSubmissionID=ClinvarSubmissionIDType(localKey=str(vc.id), submitter_date=current_date),  # TODO can localKey be a long string that's the lab ID?
                MeasureTrait=MeasureTraitType(
                    Assertion=AssertionType4(
                        AssertionType=AssertionTypeType(val_tpe="name", valueOf_=AssertionType.VARIATIONTODISEASE.value)  # TODO is this always variation to disease?
                    ),
                    ClinicalSignificance=ClinicalSignificanceType(
                        ReviewStatus=ReviewStatusType.CRITERIAPROVIDEDSINGLESUBMITTER.value,
                        DateLastEvaluated=submitter_date,
                        Description=[cs_label]
                    ),
                    ObservedIn=[ObservationSet(
                        Sample=SampleType(
                            Origin=OriginType.GERMLINE.value,
                            Species=SpeciesType(TaxonomyId=9606, valueOf_="human"),
                            Gender=vcm.get(SpecialEKeys.SEX)  # current values for sex match ClinVar's Gender
                        ),
                        Method=[MethodType(
                            MethodType=MethodTypeType(valueOf_=Methodtypelist.CURATION.value)  # TODO is this safe to always assume?
                        )],
                        ObservedData=[ObservedDataType(
                            Attribute=AttributeType(Type=MeasureAttributetype.DESCRIPTION.value, valueOf_=cautious_attempt_html_to_text(
                                text=submission_text
                            )),
                            Zygosity=zygosity,
                            Citation=pubmed_citations
                        )]
                    )]
                ),
                MeasureSet=MeasureSetType(
                    MeasureSetType=MeasureSetTypeType(val_type="name", valueOf_=Measuresettypelist.VARIANT.value),
                    Measure=[MeasureType(
                        MeasureType=MeasureTypeType(val_type="name", valueOf_=Measuretype.VARIATION.value),
                        AttributeSet=[
                            AttributeSetType(
                                MeasureTraitAttributeType=MeasureAttributeTypeType(
                                    val_type="name",
                                    valueOf_="HGVS"
                                ),
                                Attribute=AttributeType(valueOf_=vcm.get(SpecialEKeys.C_HGVS))  # TODO do we need to remove gene symbol from c.hgvs?
                            )
                        ],
                        XRef=variant_x_refs,
                        MeasureRelationship=[MeasureRelationshipType(
                            MeasureRelationshipType=MeasureRelationshipTypeType(val_type="name", valueOf_=MeasureRelationshiptype.VARIANTINGENE.value),
                            Symbol=[SetElementSetType(
                                ElementValueType=ElementValueTypeType(val_type="name", valueOf_=ElementValuetype.PREFERRED.value),
                                ElementValue=ElementValueType(valueOf_=vcm.get(SpecialEKeys.GENE_SYMBOL))
                            )]
                        )]
                    )]
                ),
                TraitSet=TraitSetType(
                    TraitSetType=TraitSetTypeType(val_type="name", valueOf_=TraitRelationshiptype.DISEASE),
                    Trait=[trait]
                )
            )
            submissions.append(submission)

        submission_set = SubmissionSetType(
            Date=export_date,
            SubmitterOfRecord=submitter_of_record,
            Submitter=[submitter],
            ClinvarSubmission=submissions
        )
        buffer = StringIO()
        submission_set.export(outfile=buffer, level=0)
        xml_str = buffer.getvalue()
        buffer.close()

        return xml_str

    def content_type(self) -> str:
        return 'text/xml'

    def filename(self) -> str:
        if self.filename_override:
            return self.filename_override
        return self.generate_filename(include_genome_build=False, extension='xml', suffix='clinvar')
