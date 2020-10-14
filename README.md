# VariantGrid

VariantGrid is a database and web application for storing, analysing and classifying variants.

Upload VCFs to:

* Automatically annotate with Ensembl VEP
* See samples from all VCFs that share a variant
* Analyse and filter samples, trios or cohorts via real time drag & drop interactive analyses
* Classify variants using customisable ACMG form
* GRCh37 (hg19), GRCh38, Ensembl and RefSeq in the same database (some automatic conversion)
* Manage and curate data, including patient phenotypes

Research use is free, while commercial use requires a [licence](https://github.com/SACGF/variantgrid/blob/master/LICENCE.md), before code becomes fully free/open source in 4 years. This is our attempt to obtain sustainable funding to for future development. See [Licence discussion](https://github.com/SACGF/variantgrid/wiki/Licence-discussion)

![VG screenshots](https://user-images.githubusercontent.com/763201/95926363-a0c6b580-0e03-11eb-9a5b-1d48e0e46722.png)

## Use cases

Private servers:

* Clinical/diagnostic use by [SA Pathology](https://www.sapathology.sa.gov.au), the South Australian public pathology service. 
* Research exomes by labs at the [Centre for Cancer Biology](https://www.centreforcancerbiology.org.au)

Public data sharing sites:

* https://shariant.org.au - Australian Genomics variant classification sharing server
* https://runx1db.runx1-fpd.org - Rare blood disease data sharing

## Other resources

This wiki is for a technical audience looking to install/modify VariantGrid code. For users of the software, please see:

* [User guide at read the docs](https://variantgrid.readthedocs.io/en/latest/)
* [VariantGrid.com](https://variantgrid.com) main site, and research cloud server

## Installation

* Contact us and we can setup and manage a cloud instance for you
* [Install](https://github.com/SACGF/variantgrid/wiki/Install) from source/scratch
* [Clone a VM](https://github.com/SACGF/variantgrid/wiki/Clone-a-VM) (AWS, VirtualBox, NECTAR, NCBI instances)

## System design

Written in Python3, using Django and PostgreSQL. See [wiki for technical details](https://github.com/SACGF/variantgrid/wiki)
