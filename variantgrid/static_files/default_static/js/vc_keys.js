let EKey = (function() {
   
    let EKey = function(keyData, index) {
        Object.assign(this, keyData);
        this.safeKey = keyData.key.replace(':','_');
        this.label = keyData.label || EKey.prettyKey(keyData.key);
        this.value_type = keyData.value_type || 'F';
        this.evidence_category = keyData.evidence_category || 'U';
        this.exclude_namespace = keyData.exclude_namespace || false;
        this.index = index;
    };
   
    EKey.prototype = {
    
        matchesFilter: function(val) {
            if (val.startsWith('#')) {
                return val.substring(1) == this.key.toLowerCase();
            }
            return this.key.toLowerCase().indexOf(val) != -1 ||
                this.label.toLowerCase().indexOf(val) != -1;
        },
    
        matchingOption: function(val) {
            var options = this.options || [];
            return options.find(o => {
                if ((val === '' || val === null) && (!('key' in o) || o.key === '')) {
                    return o;
                } else if (o.key === val) {
                    return o;
                }
            }) || null;
        },
        
        namespace: function() {
            let dividerIndex = this.key.indexOf(':');
            if (dividerIndex == -1) {
                return null;
            } else {
                return this.key.substring(0, dividerIndex);
            }
        },
    
        prettyValue: function(val) {
            let isBlank = false;
            let isValidValue = true;
            let freeText = false;
            if (typeof(val) === 'undefined') {
                val = null;
            }

            // note only M should ever have Multiplevalues but considering data types have been changed on occasion
            // treat all drop downs like they could receive an array (it can be up to server validation to mark issues
            // in that case)
            if (this.value_type === 'S' || this.value_type === 'C' || this.value_type === 'M') {
                if (val === null) {
                    // check for special blank option
                    isBlank = true;
                    let matchingOpt = this.matchingOption(val);
                    if (matchingOpt) {
                        val = matchingOpt.label || EKey.prettyKey(matchingOpt.key);
                    }
                    return {val: val, isBlank: true, isValidValue: true};
                }

                let formattedValues = [];
                if (!Array.isArray(val)) {
                    val = [val];
                }
                for (let v of val) {
                    let matchingOpt = this.matchingOption(v);
                    if (matchingOpt) {
                        v = matchingOpt.label || EKey.prettyKey(matchingOpt.key);
                    } else {
                        isValidValue = !!this.allow_custom_values;
                    }
                    formattedValues.push(v);
                }
                val = formattedValues.join(', ');

            } else if (val === null) {
                return {val: null, isBlank: true};
            } else if (this.value_type === 'B') {
                if (val !== true && val !== false) {
                    isValidValue = false;
                }
            }
                
            if (val === true) {
                val = "True";
            } else if (val === false) {
                val = "False";
            }
            
            return {
                val: val,
                isBlank: isBlank,
                isValidValue: isValidValue
            };
        },
        
        value: function(val) {
            return this.prettyValue(val).val;
        },
        
        defaultValue: function() {
            let default_option = this.options.find(o => o['default']);
            if (default_option) {
                return default_option.key;
            } else {
                return null;
            }
        },

        asSafeHtml: function(val) {
            let html = $('<span>' + val + '</span>');
            let allTags = html.find('*');
            for (let tag of allTags) {
                if (!EKey.HTML_WHITE_LIST.has(tag.tagName)) {
                    throw Error(`unsafe HTML tag ${tag.tagName}`);
                }
                if (tag.attributes.length) {
                    throw Error(`tag ${tag.tagName} contains attributes`);
                }
            }
            return html;
        },

        safeHtml: function(val, elem) {
            try {
                elem.html(this.asSafeHtml(val));
            } catch (e) {
                console.log(e);
                elem.text(val);
            }
        },
    
        formatValue: function(value, element, report) {
            if (typeof(element) === 'undefined') {
                element = $('<div>', {});
            }
            let {val, isBlank, isValidValue} = this.prettyValue(value);
            // free text or text area can have regexes
            if (val === '' || val === null || typeof(val) === 'undefined') {
                element.text('');
                element.addClass('blank');
            } else if (this.value_type === 'C' && (value === 'NM' || value === 'NA')) {
                element.text(val);
                element.addClass('not-met');
            } else if (this.value_type === 'C' && value !== this.defaultValue()) {
                element.text(val);
                element.addClass('override');
                if (!element.attr('title')) {
                    element.attr('title', 'Using an override strength for this ACMG criterion');
                }
            } else if (this.value_type === 'T') {
                this.safeHtml(val, element);
                element.removeClass('blank');
            } else if (this.value_type === 'F') {
                element.text(val);
                element.removeClass('blank');
            } else if (this.value_type === 'N' || this.value_type === 'P') {
                if (report) {
                    let numValue = Number(val);
                    if (val && !isNaN(numValue)) {
                        numValue = toPercent(numValue, this.value_type === 'P' ? 1 : null);
                        element.text(numValue);
                        element.attr('title', toFixedString(val));
                    } else {
                        element.text(val);
                    }
                } else {
                    val = toFixedString(val);
                    element.text(val);
                }
            } else {
                // TODO - add formatting if not isValidValue
                element.text(val);
                if (isBlank) {
                    element.addClass('blank');
                } else {
                    element.removeClass('blank');
                }
                if (value === false) {
                    element.addClass('not-met');
                }
            }
            if (!isValidValue) {
                element.addClass('error-value');
            } else {
                element.removeClass('error-value');
            }
            return element;
        }
    };
   
    return EKey;
    
})();

EKey.HTML_WHITE_LIST = new Set(['BR','B','I','OL','UL','LI','U']);

EKey.critValues = {
    "NM": "Not Met",
    "NA": "Not Applicable",
    "BA": "Benign Standalone",
    "BS": "Benign Strong",
    "BP": "Benign Supporting",
    "N": "Neutral",
    "PP": "Pathogenic Supporting",
    "PM": "Pathogenic Moderate",
    "PS": "Pathogenic Strong",
    "PVS": "Pathogenic Very Strong"
};

EKey.families = {
    V: 'Variant',
    H: 'Gene',
    HP: 'Patient',
    HT: 'Test',
    P: 'Population data',
    CP: 'Computational and predictive data',
    F: 'Functional data',
    S: 'Segregation data',
    D: 'De novo data',
    A: 'Allelic data',
    DB: 'Other database',
    O: 'Other data',
    L: 'Literature',
    HI: 'Interpretation',
    SO: 'Sign Off',
    U: 'Unknown'
};

EKey.prettyKey = function(key) {
    if (!key) {
        return '';
    }
    if (/^([A-Z|a-z])[_-]/.test(key)) {
        key = key.substring(0,1) + '-' + key.substring(2);
    }
    return key.substring(0,1).toUpperCase() + key.substring(1).replace(/_/g,' ')
};

var EKeys = (function() {
    
    let EKeys = function(loadedKeys) {
        this._map = new Map();
        this.hasUnknown = false;
        loadedKeys.forEach((data, index) => {
            this._set(new EKey(data, index));
        });
    };
    
    EKeys.prototype = {
        
        configCopy(labConfig) {
            labConfig = labConfig || {};
            let namespaces = new Set(labConfig['namespaces']);
            namespaces.add(null);
            delete labConfig['namespaces'];

            return new EKeys(Array.from(this._map.values()).map(ekey => {
                let config = labConfig[ekey.key];
                if (config == true || config == false) {
                    config = {hide: config};
                }
                let exclude_namespace = !namespaces.has(ekey.namespace());
                if (exclude_namespace) {
                    config = config || {};
                    config.exclude_namespace = true;
                    config.mandatory = false;   
                }
                if (config) {
                    data = Object.assign({}, ekey, config);
                    let newkey = new EKey(data, ekey.index);
                    return newkey;
                } else {
                    return ekey;
                }
            }));
        },
    
        key(key) {
            var eKey = this._map.get(key);
            if (!eKey) {
                eKey = new EKey({key: key}, this._map.length);
                this._set(eKey);
            }
            return eKey;
        },
        
        _set(eKey) {
            this._map.set(eKey.key, eKey);
            this.hasUnknown = this.hasUnknown || eKey.evidence_category == 'U';
        },
        
        forEach(funkey) {
            return this._map.forEach(funkey);
        },
        
        criteria() {
            var c = [];
            this._map.forEach(e => {
                if (e.value_type === 'C' && !e.exclude_namespace) {
                    c.push(e);
                }
            });
            return c;
        },

        criteriaIgnoreNamespace() {
            var c = [];
            this._map.forEach(e => {
                if (e.value_type === 'C') {
                    c.push(e);
                }
            });
            return c;
        },
        
        renderClinicalSignificance(cs) {
            let {val, isBlank, isValidValue} = this.key('clinical_significance').prettyValue(cs);
            if (isBlank) {
                cs = 'unclassified';
                val = 'Unclassified';
            }
            return $('<span>', {class: `cs cs-${cs.toLowerCase()}`, text: val});
        },

        keySelectOptions() {
            let alphabeticKeys = [];
            this.forEach(e => {
                alphabeticKeys.push(e);
            });
            alphabeticKeys.sort((a, b) => {
                return a.label.localeCompare(b.label);
            });
            let options = alphabeticKeys.map(e => {
                return $('<option>', {value: e.key, text:e.label});
            });
            options.splice(0, 0, $('<option>'));
            return options;
        }
    };
    
    return EKeys;
})();

// cache keys for a minute
EKeys.cacheTimeMs = 60000;

EKeys.load = function() {
    if (EKeys.cachedKeys) {
        return Promise.resolve(EKeys.cachedKeys);
    }

    cached = window.localStorage.getItem('evidenceKeys');
    use_cache = null;
    if (cached) {
        try {
            cached = JSON.parse(cached);
            ageMs = Date.now() - cached.date;
            if (ageMs <= EKeys.cacheTimeMs) {
                use_cache = cached;
            }
        } catch (e) {
            console.log(e);
        }
    }
    if (use_cache) {
        EKeys.cachedKeys = new EKeys(use_cache.data);
        return Promise.resolve(
            EKeys.cachedKeys
        );
    } else {
        return $.getJSON( Urls.evidence_keys_api() ).then(keys => {
            window.localStorage.setItem('evidenceKeys', JSON.stringify({
                date: Date.now(),
                data: keys
            }));
            EKeys.cachedKeys = new EKeys(keys);
            return EKeys.cachedKeys;
        });
    }
};

EKeys.levelToIndex = {
    lab: 0,
    institution: 1,
    logged_in_users: 2,
    'public': 3
};

EKeys.shareLevelInfo = function(share_level, record, defaultToInstitution) {
    var base = '/static/icons/share_level/';
    if (!share_level || share_level == 'user') {
        return { icon: base + 'draft.png', title: 'Last Edited'};
    }
    
    let contentGenerator = () => {
        switch (share_level) {
            case 'latest': return { src: base + 'draft.png', title: 'Last Edited'};
            case 'lab':
                return {
                    icon : base + 'lab.png',
                    title : record && record.lab_name ? record.lab_name : 'Lab'
                };
            case 'institution':
                return {
                    icon: base + 'institution.png',
                    title : record && record.institution_name ? record.institution_name : VcSettings.INSTITUTION
                };
            case 'logged_in_users':
                return {
                    icon: base + 'logged_in_users.png',
                    title: VcSettings.LOGGED_IN_USERS
                };
            case 'public':
                return {
                    icon: base + 'public.png',
                    title: VcSettings.PUBLIC
                };
            default: return {icon: base + share_level, title: 'Shared with ' + share_level};
        }
    };
    let content = contentGenerator();
    if (record) {
        let consider_level = EKeys.levelToIndex[share_level] || 0;
        let current_level = EKeys.levelToIndex[record.publish_level] || (defaultToInstitution ? 1 : 0);
        content.current = current_level == consider_level;
        content.included = consider_level < current_level;
    }
    return content;
};

EKeys.fixDescription = function(htmlText) {
    htmlText = htmlText.trim();
    htmlText = htmlText.replaceAll(/(<br>|<br\/>|<br \/>)/gm, '\n');
    htmlText = htmlText.replaceAll(/^\s+$/gm, '');
    htmlText = htmlText.replaceAll(/\n{3,}/gm, '\n\n');
    let html = $('<span>', {html: htmlText});
    html.find('p').each(function() {
        let $this = $(this);
        let pHtml = $this.html().trim();
        if (pHtml.length === 0) {
            $this.remove();
        } else {
            $this.html(pHtml);
        }
    });
    return html;
}


SpecialEKeys = {};
SpecialEKeys.ASSERTION_METHOD = 'assertion_method';
SpecialEKeys.GENOME_BUILD = 'genome_build';
SpecialEKeys.REFSEQ_TRANSCRIPT_ID = 'refseq_transcript_id';
SpecialEKeys.CONDITION = 'condition';
SpecialEKeys.CLINICAL_SIGNIFICANCE = 'clinical_significance';
SpecialEKeys.GENE_SYMBOL = 'gene_symbol';
SpecialEKeys.GENE_OMIM_ID = 'gene_omim_id';
SpecialEKeys.UNIPROT_ID = 'uniprot_id';
SpecialEKeys.CLINVAR_VARIANTION_ID = 'clinvar_variation_id';
SpecialEKeys.VARIANT_COORDINATE = 'variant_coordinate';
SpecialEKeys.HGNC_ID = 'hgnc_id';
SpecialEKeys.CLINGEN_ALLELE_ID = "clingen_allele_id";
SpecialEKeys.C_HGVS = 'c_hgvs';
SpecialEKeys.P_HGVS = 'p_hgvs';
SpecialEKeys.G_HGVS = 'g_hgvs';
SpecialEKeys.SEARCH_TERMS = 'search_terms';