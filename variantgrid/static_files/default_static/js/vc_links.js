let VCLink = (function() {

    let VCLink = function(params) {
        this.text = params.text;
        this.href = params.href;
        this.missing = params.missing;
        this.title = params.title;
        this.build = params.build;
    };

    VCLink.prototype = {
        isMissing() {
            return !!this.missing;
        },

        asAnchor(style) {
            let title = this.title || this.text;
            if (this.build) {
                title = `${title} (${this.build})`;
            }
            let cssClass = 'hover-link external-link';
            if (this.isMissing()) {
                title = (title ? `${title} ` : '') + ` : Missing ${this.missing}`;
                cssClass += " general";
            }
            let attribs = {href: this.href, title: title, class: cssClass};
            if (!this.href.startsWith("javascript:")) {
                attribs.target = '_blank'; // open non-JS in new window
            }
            if (this.build) {
                let buildTitle = this.buildTitle ? this.buildTitle :
                attribs.html = [
                    $('<span>', {text: this.text + ' '}),
                    $('<span>', {class: 'build-indicator', text: this.build})
                ];
            } else {
                attribs.html = [this.text];
            }
            if (style === "bootstrap") {
                attribs["class"] = "list-group-item list-group-item-action hover-link external-link";
                if (this.isMissing()) {
                    attribs["class"] += " list-group-item-light";
                }
            }

            return $('<a>', attribs);
        }
    };
    return VCLink;

})();

let VCLinks = (function() {

    let VCLinks = function(eKeys) {
        this.eKeys = eKeys;
        this.data = {};
        this.variant_coordinate_parts = null;
    };

    VCLinks.ALL_KEYS = [
        SpecialEKeys.C_HGVS,
        SpecialEKeys.VARIANT_COORDINATE,
        SpecialEKeys.GENE_SYMBOL,
        SpecialEKeys.REFSEQ_TRANSCRIPT_ID,
        SpecialEKeys.P_HGVS,

        SpecialEKeys.CLINGEN_ALLELE_ID,
        SpecialEKeys.GENOME_BUILD,
        SpecialEKeys.UNIPROT_ID,
        SpecialEKeys.CLINVAR_VARIANTION_ID,
        SpecialEKeys.HGNC_ID,
        SpecialEKeys.GENE_OMIM_ID,
        SpecialEKeys.UNIPROT_ID,
        SpecialEKeys.SEARCH_TERMS,
    ];

    VCLinks.prototype = {

        generateLinks(data) {
            this.data = data;

            this.variant_coordinate_parts = null;
            let variant_coordinate = this.data[SpecialEKeys.VARIANT_COORDINATE];
            if (variant_coordinate) {
                this.variant_coordinate_parts = variant_coordinate.match(/(.+?):([0-9]+)(?:[ ])?([ATCG].+)/);
            }

            let links = [];
            let igvLink = this.generateIgv();
            if (igvLink) {
                // always put igvLink first as it's special
                links.push(igvLink);
            }

            links.push(this.generateBeacon());
            links.push(this.makeLink('ClinGen Allele Reg.', 'http://reg.clinicalgenome.org', '/redmine/projects/registry/genboree_registry/by_caid?caid=@@', SpecialEKeys.CLINGEN_ALLELE_ID, 'Clingen Allele Registry'));
            links.push(this.generateClingenKb());
            links.push(this.makeLink('Clinvar Variant', 'https://www.ncbi.nlm.nih.gov', '/clinvar/variation/@@', SpecialEKeys.CLINVAR_VARIANTION_ID));
            links.push(this.makeLink('Cosmic (Gene)', 'https://cancer.sanger.ac.uk/cosmic', '/gene/analysis?ln=@@', SpecialEKeys.GENE_SYMBOL));
            links.push(this.makeLink('GHR (Gene)', 'https://ghr.nlm.nih.gov', '/gene/@@', SpecialEKeys.GENE_SYMBOL, 'Genetics Home Reference'));
            // commenting out genomizer due to build confusion
            //links.push(this.generateGenomizer());
            links.push(this.generateGnomad());
            links.push(this.makeLink('Google', 'https://google.com/search', '?q=@@', SpecialEKeys.SEARCH_TERMS));
            links.push(this.makeLink('GTEx', 'https://gtexportal.org/home/', 'gene/@@', SpecialEKeys.GENE_SYMBOL, 'Genotype-Tissue Expression'));
            links.push(this.makeLink('Monarch Phen.', 'https://monarchinitiative.org', '/gene/HGNC:@@', SpecialEKeys.HGNC_ID, 'Monarch Phenotype'));
            links.push(this.makeLink('NCBI (Gene)', 'https://www.ncbi.nlm.nih.gov/gene/', '?term=@@', SpecialEKeys.GENE_SYMBOL));
            links.push(this.makeLink('OMIM (Gene)', 'https://www.omim.org', '/entry/@@', SpecialEKeys.GENE_OMIM_ID));
            links.push(this.makeLink('PDB', 'https://www.rcsb.org/pdb', '/protein/@@', SpecialEKeys.UNIPROT_ID, 'Protein Data Bank'));
            links.push(this.generateUcsc());
            links.push(this.makeLink('Uniprot ID', 'https://www.uniprot.org/uniprot/', '@@', SpecialEKeys.UNIPROT_ID));
            links.push(this.generateVarsome());

            return links;
        },

        buildName() {
            let genome_build = this.data[SpecialEKeys.GENOME_BUILD];
            if (genome_build) {
                return genome_build.indexOf('38') !== -1 ? 'GRCh38' : 'GRCh37';
            } else {
                return '';
            }
        },

        generateClingenKb() {
            let val = this.data[SpecialEKeys.HGNC_ID];
            if (val) {
                if (!val.toUpperCase().startsWith('HGNC')) {
                    val = 'HGNC:' + val;
                }
                let link = `https://search.clinicalgenome.org/kb/genes/${val}`;
                return new VCLink({text: 'ClinGen KB', href: link});
            } else {
                return new VCLink({
                    text: 'ClinGen KB',
                    href: 'https://search.clinicalgenome.org',
                    missing: this.eKeys.key(SpecialEKeys.HGNC_ID).label
                });
            }
        },

        generateBeacon() {
            // https://beacon-network.org/#/search?pos=25288616&chrom=20&allele=GGCTCTTA&ref=G&rs=GRCh37
            if (this.variant_coordinate_parts) {
                let parts = this.variant_coordinate_parts;
                let chr = parts[1];
                let pos = parts[2];
                let ref = null;
                let alt = null;
                let refAlt = parts[3].match(/([ATCG]+)>([ATCG]+)/);
                if (refAlt) {
                    ref = refAlt[1];
                    alt = refAlt[2];
                } else if (refAlt = parts[3].match(/([ATCG]+)=/)) {
                    ref = refAlt[1];
                    alt = ref;
                }
                let build = null;
                let genome_build = this.data[SpecialEKeys.GENOME_BUILD];
                if (genome_build) {
                    build = genome_build.indexOf('38') != -1 ? 'GRCh38' : 'GRCh37';
                }
                if (ref && alt && build) {
                    let link = `https://beacon-network.org/#/search?pos=${encodeURIComponent(pos)}&chrom=${encodeURIComponent(chr)}&allele=${encodeURIComponent(alt)}&ref=${encodeURIComponent(ref)}&rs=${encodeURIComponent(build)}`;
                    return new VCLink({
                        text:'Beacon',
                        href:link,
                        build: this.buildName()
                    });
                }
            }
            return new VCLink({
                text:'Beacon',
                href: `https://beacon-network.org`,
                missing: 'or complex ' + this.eKeys.key(SpecialEKeys.VARIANT_COORDINATE).label
            });
        },

        generateGnomad() {
            if (this.variant_coordinate_parts && this.data[SpecialEKeys.GENOME_BUILD]) {
                let parts = this.variant_coordinate_parts;
                let genome_build = this.data[SpecialEKeys.GENOME_BUILD];
                let safeUrl = parts[1] + '-' + parts[2] + '-' + (parts[3] || '').replace('>', '-');
                let dataset = genome_build.indexOf('38') != -1 ? 'gnomad_r3' : 'gnomad_r2_1';
                return new VCLink({
                    text:'gnomAD',
                    href:`https://gnomad.broadinstitute.org/variant/${safeUrl}?dataset=${dataset}`,
                    build: this.buildName()
            });
            }
            return new VCLink({
                text: 'gnomAD',
                href: `https://gnomad.broadinstitute.org`,
                missing: 'or complex ' +this.eKeys.key(SpecialEKeys.VARIANT_COORDINATE).label
            });
        },

        generateUcsc() {
            let genome_build = this.data[SpecialEKeys.GENOME_BUILD];
            if (this.variant_coordinate_parts && genome_build) {
                try {
                    let parts = this.variant_coordinate_parts;
                    let coordinate = parseInt(parts[2]);
                    let use_build = 'hg19';
                    if (genome_build && genome_build.indexOf('38') != -1) {
                        use_build = 'hg38';
                    }
                    let range = `chr${parts[1]}%3A${coordinate - 20}-${coordinate + 20}`;
                    let url = `https://genome.ucsc.edu/cgi-bin/hgTracks?db=${use_build}&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=${range}`;
                    return new VCLink({
                        text:'UCSC',
                        href:url,
                        build: this.buildName()
                    })
                } catch (e) {
                }
            }
            return new VCLink({
                text: 'UCSC',
                href: `https://genome.ucsc.edu`,
                missing: 'or complex ' + this.eKeys.key(SpecialEKeys.VARIANT_COORDINATE).label
            });
        },

        generateVarsome() {
            let genome_build = this.data[SpecialEKeys.GENOME_BUILD];
            if (genome_build) {
                let use_build = 'hg19';
                if (genome_build && genome_build.indexOf('38') != -1) {
                    use_build = 'hg38';
                }
                let searchFor = this.data[SpecialEKeys.VARIANT_COORDINATE];
                if (searchFor == null) {
                    searchFor = this.data[SpecialEKeys.C_HGVS];
                }
                if (searchFor) {
                    return new VCLink({
                        text: 'varsome',
                        href: `https://varsome.com/variant/${use_build}/${encodeURIComponent(searchFor)}`,
                        build: this.buildName()
                    });
                }
            }
            return new VCLink({
                text: 'varsome',
                href: `https://varsome.com`,
                missing: this.eKeys.key(SpecialEKeys.C_HGVS).label
            });
        },

        generateIgv() {
            if (typeof createIgvUrl == 'undefined' || !this.variant_coordinate_parts) {
                return null;
            }
            let locus = this.variant_coordinate_parts[1] + ':' + this.variant_coordinate_parts[2];
            return new VCLink({text:'IGV', href:createIgvUrl(locus), title:`Open ${locus} in IGV`});
        },

        makeLink(text, link, param, key, title) {
            let val = this.data[key];
            if (val) {
                link = link + param.replace('@@', encodeURIComponent(val));
                return new VCLink({text:text, href:link, title: title});
            } else {
                return new VCLink({text:text, href:link, missing:this.eKeys.key(key).label, title: title});
            }
        }
    };

    return VCLinks;

})();