# Classifications

A classification is the review of a single variant to determine its likelihood in causing a disease.
It is the 2nd last step in a variant analysis followed by generating a report.

Note this documentation often refers to the Shariant project as that's the environment that utilises the advanced functions around Classifications the most, but everything here can be utilised by any VariantGrid environment if the features were to be enabled.

For the Shariant project, classifications is the key higher level app as classifications are imported by many labs, compared to each other, shared, uploaded to ClinVar.

## Main Classes / Concepts

## Classification
Relates to what would have been one physical form of a classification.
Keeps a record of the latest evidence in JSON format.

Each piece of evidence is related to an evidence key, can have a value, a free text note, and an explain (typically auto-generated text to describe if a lab uses that field in a non-standard way).

### ClassificationModification
Classifications are constantly updated, and each set of changes will make a new ClassificationModification linked to the parent Classification. Some Modifications will be "Published" at a certain share level, published versions are the only ones that can be seen by users that don't belong to the lab that created the classification.
When a classification is published at a certain share level, it can only be published at that share level or higher from then on.

### ImportedAlleleInfo
When importing classifications we use a combination of its genome build, c.HGVS (or g.HGVS) to resolve it to an Allele. We do this by first linking the classification to an ImportedAlleleInfo for that combination of data, and resolving that (if it hasn't already been resolved). This object also links to ResolvedVariantInfo which has cached representations of 37 and 38 c.HGVS values.

ImportedAlleleInfo also keeps track of if there were any resolution warnings or errors, and classifications associated with errored ImportedAlleleInfos will not be exported out of the system.

### Withdrawing
Withdrawing a classification is simply a soft delete.

### Evidence Key
There's hundreds of potential bits of information that up the evidence in a classification.
Evidence keys are used to generate the web form and validate uploaded data.

The evidence key determines which section of the form the field appears in, order, help text, and the type of field (integer, text, multi-line text, select, criteria) as well as some meta-data (such as ClinVar values, buckets for discordances).

The Evidence Key configuration 

#### Criteria

Criteria is a value type for an evidence key that is given special prominence. For the classic examples see ACMG Curation Guidelines re PM1, PM2, etc.
Criterias will be given values from Benign Stand Alone, Benign Strong, Benign Moderate, Benign Supporting, Neutral, Pathogenic Supporting etc and typically labs use the accumulation of these to determine the overall clinical significance.

#### Clinical Significance

This is the overall value of a classification that goes from Benign, Likely Benign, VUS (Variant of Uncertain Significance), Likely Pathogenic, Pathogenic. In addition, a lot of environments are configured to have different levels of VUS, VUS_A, VUS_B, VUS_C with VUS_B being in the middle and VUS_A being more pathogenic.

## Discordance

When two or more classifications for the same Allele have their clinical significances fall into different "buckets" it is considered a Discordance (Discordances can be enabled/disabled in the settings file).
When this is first detected a DiscordanceReport is created. At the time that all (non-withdrawn) classifications are in the same bucket, the DiscordanceReport is closed.

### Pending Concordance

When reviewing a DiscordanceReport, users are able to agree to change their clinical significances. This raises pending change flags, with the idea that the users will then change the data in their own curation systems and upload those changes sometime in the future.

## Condition Text Matching

Much of the data imported into Shariant is provided with textual descriptions of the condition, where for our own purposes and for uploading to ClinVar, a standard ontology term is greatly preferred.

This is done by looking at the text on the condition, and building up a hierarchy of
* "Lab / Context Text"
* ...."Gene Symbol"
* ........"Mode of Inheritance"
There is some auto matching that can be done, but users are then able to match free text to standard conditions (at any one of these levels), and if that free text is seen again at a level that has already been set, the standard condition is assigned to the classification.

## Importing

Classification data in Shariant is generally obtained by taking exports from the contributing labs curation systems, transforming that data to the Shariant format using the OmniImporter project.
There is also more live data through APIs, which is done with the Syncing app which is discussed later.

### OmniImporter

As even labs that use the same curation system provide data in very different formats, the most efficient method of converting the data is for the labs to provide us their un-converted files, the OmniImporter project convert it to VariantGrid JSON format (as well as tracking).

The OmniImporter has a class ClassificationRecord which is generated from EvidenceKeys, allowing type safety to be enforced during the conversion process.

In the Shariant setup, labs will upload a file to Shariant (which will then save it to Amazon's s3 bucket). When the file is processed, Shariant will download the file and place it in a directory accessible by the OmniImporter, and then invoke the OmniImporter with parameters about the files location, what share level the data should be imported as, the organisation and lab doing the importing.

The OmniImporter picks the appropriate converter based on the org and lab, makes a HTTP connection back to Shariant, and sends the converted data.

### Syncing

Syncing is it's own Django app for automated upload and download of data.
Examples include SA Path's VariantGrid uploading/downloading to/from Shariant - utilising Shariant's import process, as well as its export process. In this mode each classification to be uploaded is individually tracked to see if the last uploaded version is up to date.

Also introduced in 2023 was syncing to Alissa (a commonly used curation system), where all the relevant data from Alissa has to be downloaded each time, and then passed onto OmniImporter. When uploading to Alissa check the last successful since date of upload, then using the export functionality in the MVL JSON format.

## Export

VariantGrid can export its classification data in a myriad of formats, with the intent that this data can then be imported into other curation systems.
Export formats include Alissa's MVL, CSV (which we can review for general validation), JSON (can be used to import into other VariantGrid systems), VCF, REDCap, ClinVar Compare (which shows where a lab's data differs from what's present in ClinVar)

## ClinVar Export

Exporting to ClinVar was a little too complicated to perform via the Syncing app. This is due to a combination of factors such as labs wanting to review data before it is exported, ClinVar only wanting one record per clinvarkey/allele/condition, being able to list all validation errors as to why some records are too incomplete to upload to ClinVar etc.

The process is:
Review classifications to find the latest record for each clinvar key, allele, condition combination.
* Convert that classification data to ClinVar JSON.
* Upload the ClinVar JSON in a batch, versioning off the data.
* Poll ClinVar to see if that data has been processed yet
* Process the response file to assign ClinVar IDs (known as SCVs) to the uploaded data so we can distinguish between updates and creation of new records.

### ClinVar Key

When exporting to ClinVar, a private ClinVar Key is required which also identifies the organisation that is providing the data. Within VariantGrid this is called a "ClinVar Key". A single key can be optionally shared between labs (such as between SA Path labs) or be unique to a lab.

### ClinVarExportConverter

This class takes a ClassificationModification and converts it the JSON used by VariantGrid. It uses the ValidatedJson library which is able to embed warnings/errors in the JSON for display purposes, but also to export clean JSON to ClinVar.
