## VariantGrid

A [source available](LICENCE.md) genetic database and web analysis platform written in Python/Django and Postgresql.

Developed by SA Pathology (South Australian Health) for research or diagnostic labs to store, analyse and report on 10-100k exomes (VCF).

Currently used by SA Pathology, UniSA Centre for Cancer Biology research labs, RUNX1db, the RUNX1 foundation data sharing site, and Shariant, the Australian Genomics variant classification sharing server.

## FEATURES

* Open and extendible, API to integrate with other systems
* Create custom analyses in real time via drag and drop interface
* Imports VCF and automatically annotates variants with VEP + a range of plugins.  
* Analyse singletons, trios, cohorts. See all variants in a gene, or everyone who has a variant
* Classify and report on clinically relevant variants
* GRCh37 (hg19), GRCh38, Ensembl and RefSeq in the same database (some automatic conversion)
* Many QC features and graphs

## INSTALL

Use an existing public server (research use only)
Clone a VM (AWS, Nectar or NCI)
Step by step guide to install from scratch
Buy a support contract for an intranet or cloud server

## TECH

Python 3.8, Django 3.1 and Postgresql. Async celery tasks.

See technical guide in the wik.
