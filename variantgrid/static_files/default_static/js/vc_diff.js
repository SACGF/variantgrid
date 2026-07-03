//
// DIFF
//
// Render a table of rows that are different by VCForm fields

const Diff = (function() {

    const Diff = function(records, diffDom, eKeys, flagParams) {
        
        this.flagParams = flagParams;
        this.versions = records.map(r => {
            const merged = Object.assign({}, r, r.data); 
            delete merged['data'];
            return merged;
        });
        let comparingDiffRecords = false;
        let firstId = null;
        for (const record of records) {
            if (!firstId) {
                firstId = record.id;
            } else if (record.id !== firstId) {
                comparingDiffRecords = true;
                break;
            }
        }

        const groupMap = {};
        const criteriaMap = {};
        for (const familyKey of Object.keys(EKey.families)) {
            groupMap[familyKey] = {"label": EKey.families[familyKey], eKeys: []};
            criteriaMap[familyKey] = {eKeys: []};
        }
        
        eKeys.forEach(k => {
            const key = k.key;
            const isCriteria = k.value_type === 'C';
            const category = k.evidence_category;

            const dict = (isCriteria ? criteriaMap : groupMap)[category];
            if (dict) {
                dict.eKeys.push(k);
            }
        });

        const groups = Object.keys(groupMap).map(g => {
            const combined = Object.assign({}, groupMap[g]);
            combined.eKeys = criteriaMap[g].eKeys.concat(combined.eKeys);
            return combined;
        }).filter(g => g.eKeys.length);

        this.diffDom = diffDom;
        this.eKeys = eKeys;
        this.groups = groups;

        this.showAllAcmg = comparingDiffRecords;

        this.excluded = {};
    };

    Diff.LITERATURE_DBS = {
        'Bookshelf ID': true,
        'PMCID': true,
        'PMID': true
    };

    Diff.hash = function(val) {
        val = Diff.emptyToNull(val);
        if (val === null || val === 'NM' || val === 'NA') {
            // treat Not Met the same as null
            return null;
        } else if (Array.isArray(val)) {
            val.sort();
        } else {
            val = [val];
        }
        return val.map(e => `${e}`).join('&%');
    };
    
    Diff.emptyToNull = function(val) {
        if (val === null || typeof(val) === 'undefined') {
            return null;
        }
        if (Array.isArray(val) && val.length === 0) {
            return null;
        }
        if (typeof(val) === 'string' && val.trim().length === 0) {
            return null;
        }
        return val;
    };
    
    Diff.compareDbRefs = function(db1, db2) {
        if (db1.db && db2.db) {
            const sortDb = db1.db.localeCompare(db2.db);
            if (sortDb) {
                return sortDb;
            }
        }
        try {
            return parseFloat(db1.idx) - parseFloat(db2.idx);
        } catch (e) { }
        return db1.id.localeCompare(db2.id);
    },
    
    Diff.prototype = {
       
        isIncluded(v) {
            let uniqueId = `${v.id}`;
            if (v.version) {
                uniqueId += `.${v.version}`;
            }
            return !this.excluded[uniqueId];
        },

        toggledCriteriaHiding() {
            this.showAllAcmg = !!$('#show_all_acmg').prop('checked');
            this.applyCriteriaHiding();
        },

        applyCriteriaHiding() {
            if (this.showAllAcmg) {
                $('.diff tr.empty.criteria').show();
            } else {
                $('.diff tr.empty.criteria').hide();
            }
        },
       
        configureColumns() {
            const modalContent = createModalShell('configureColumnsModal', 'Choose which columns to compare');
            const content = modalContent.find('.modal-body');

            const compareBox = $('<div>').appendTo(content);
            this.versions.forEach((v, index) => {
                let uniqueId = `${v.id}`;
                if (v.version) {
                    uniqueId += `.${v.version}`;
                }
                
                const checkbox = $('<input>', {type:'checkbox', class:"form-check-input", id:`ch-${uniqueId}`, click:event => {
                    const checked = $(event.target).prop('checked');
                    this.excluded[uniqueId] = !checked;
                    this.render();
                }});
                if (!this.excluded[uniqueId]) {
                    checkbox.prop('checked', true);
                }
                const label = $('<label>', {class:'form', for:`ch-${uniqueId}`,  text: `${v.org_name} / ${v.lab_name} / ${v.cr_lab_id}`});

                compareBox.append($('<div>', {
                    class: 'form-group form-check',
                    html: [checkbox, label]
                }));
            });
            modalContent.on('hidden.bs.modal', function() {
                modalContent.modal('dispose');
                modalContent.remove();
            });
            modalContent.modal('show');
        },
               
        dbRefCompare(versionDbRefs, source) {
            // Move this out into its own function
            const allDbRefs = {};
            versionDbRefs.forEach((dbRefs, index) => {
                if (dbRefs) {
                    dbRefs.forEach(dbRef => {
                        let meta = allDbRefs[dbRef.id];
                        if (!meta) {
                            meta = { 'dbRef': dbRef, whoHas: {} };
                            allDbRefs[dbRef.id] = meta;
                            meta.source = source;
                        }
                        meta.whoHas[index] = true;
                        if (Object.keys(meta.whoHas).length === versionDbRefs.length) {
                            meta.everyoneHas = true;
                        }
                    });
                }
            });
            const orderedDbRefs = Object.keys(allDbRefs).map(k => allDbRefs[k]);
            return orderedDbRefs.sort((meta1, meta2) => Diff.compareDbRefs(meta1.dbRef, meta2.dbRef));
        },

        render: function() {
            const e = $(this.diffDom);
            const includedVersions = this.versions.filter(v => this.isIncluded(v));
           
            const table = $('<table>', {class: 'diff table'});
            const thead = $('<thead>').appendTo(table);
            const tbody = $('<tbody>').appendTo(table);
            const headerRow = $('<tr>', {class: 'versions no-compare'}).appendTo(thead);

            // top left cell
            $('<th>', {html: [
                $('<div>', {html: [
                    $('<a>', {text: 'Configure columns', class:'hover-link', click: ()=>this.configureColumns()})
                ]}),
                $('<br/>'),
                $('<div>', {class:'form-inline', html:
                    $('<div>', {class: 'form-check', html: [
                        $('<label>', {class:'form-check-label', html:[
                            $('<input>', {class:'form-check-input', type: 'checkbox', id: 'show_all_acmg', click: ()=>this.toggledCriteriaHiding(), 'checked': this.showAllAcmg}),
                            ' Show unmet ACMG criteria'
                        ]})
                    ]})
                })
            ]}).appendTo(headerRow);


            let renderDbRefCompare = function(dbRefCompared) {
                dbRefCompared.forEach(dbRefC => {
                    const dbRef = dbRefC.dbRef;
                    const whoHas = dbRefC.whoHas;

                    let rowClass;
                    let icon;
                    let heading;
                    if (dbRef.label) {
                        rowClass = 'option-row';
                        heading = dbRef.label;
                    } else {
                        rowClass = 'citation-row';
                        heading = CitationsManager.defaultManager.citationDomFor(dbRef, true);
                    }

                    const row = $('<tr>', {class: `${rowClass}`}).appendTo(tbody);

                    if (!dbRefC.everyoneHas) {
                        row.addClass('blanks');
                        icon = $('<i class="fas fa-not-equal" style="color:#b88"></i>');
                    } else {
                        icon = $('<i class="fas fa-check-square" style="color:#6d6"></i>');
                    }

                    const th = $('<th>', {html: [icon, heading]}).appendTo(row);

                    let hasCount = 0;
                    let hasNotCount = 0;
                    includedVersions.filter(v => this.isIncluded(v)).forEach((v, index) => {
                        let uniqueId = `${v.id}`;
                        if (v.version) {
                            uniqueId += `.${v.version}`;
                        }
                        if (whoHas[index]) {
                            hasCount++;
                            let title = null;
                            const whoHasData = whoHas[index];
                            if (Array.isArray(whoHasData)) {
                                title = 'Referenced in ' + whoHasData.join(', ');
                            }
                            $('<td>', {
                                title: title,
                                html: $('<i class="far fa-check-square"></i>')
                            }).appendTo(row);
                           //$('<td>', {text: 'included', class: primeInclusion ? '' : 'addition'}).appendTo(row);
                       } else {
                            hasNotCount++;
                            $('<td>', {html: '<i class="far fa-square" style="opacity:0.2""></i>'}).appendTo(row);
                        }
                    });
                    const tooltip = [];
                    const referenceType = dbRefC.source === 'multiselect' ? 'value' : 'reference';
                    if (hasCount) {
                        const verb = hasCount === 1 ? 'record has this' : 'records have this';
                        tooltip.push(`${hasCount} <span style="color:#444;font-size:smaller">x</span> ${verb} ${referenceType}`);
                    }
                    if (hasNotCount) {
                        const verb = hasNotCount === 1 ? 'record does not have this' : 'records do not have this';
                        tooltip.push(`${hasNotCount} <span style="color:#444;font-size:smaller">x</span> ${verb} ${referenceType}`);
                    }
                    th.attr('title', dbRef.label || dbRef.id);
                    th.attr('data-content', tooltip.join('<br>'));
                });
            };
            renderDbRefCompare = renderDbRefCompare.bind(this);
            //new Citations().refresh();

            //var percent = Math.floor(100 / Math.max(1, this.versions.length));
            includedVersions.forEach((v, index) => {

                const headers = [];

                let url = Urls.view_classification(v.id);
                if (v.version) {
                    url += `.${v.version}`;
                }

                const titleDom = $('<a>', {class:'hover-link text-center d-flex flex-column flex-align-center', href: url});
                headers.push(titleDom);

                const titlePart = $('<div>', {text:  v.org_name + ' / ' + v.lab_name + ' / ' + v.cr_lab_id}).appendTo(titleDom);
                if (v.first_seen) {
                    const first_seen_date = moment(v.first_seen * 1000).format('DD/MMM/YYYY HH:mm');
                    const last_seen_date = moment(v.version * 1000).format('DD/MMM/YYYY HH:mm');
                    if (first_seen_date !== last_seen_date) {
                        const title = `First uploaded at ${first_seen_date},<br/>Non-significant change at ${last_seen_date}`;
                        $('<div>', {class: 'timestamp', text: first_seen_date + '*', 'data-content': title}).appendTo(titleDom);
                    } else {
                        $('<div>', {class: 'timestamp', text: moment(v.version * 1000).format('DD/MMM/YYYY HH:mm')}).appendTo(titleDom);
                    }
                } else if (v.version) {
                    $('<div>', {class: 'timestamp', text: moment(v.version * 1000).format('DD/MMM/YYYY HH:mm')}).appendTo(titleDom);
                } else {
                    $('<div>', {class: 'timestamp', text: 'Working version', title: 'Working Version'}).appendTo(titleDom);
                }
                let content = null;

                const clinSigText = [];

                const clin_sig_key = this.eKeys.key(SpecialEKeys.CLINICAL_SIGNIFICANCE);
                const clin_sig = clin_sig_key.prettyValue((v.clinical_significance || {}).value);
                if (clin_sig.val) {
                    clinSigText.push(clin_sig.val);
                }

                const somatic_clin_sig_value = v[SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE];
                if (somatic_clin_sig_value && somatic_clin_sig_value.value) {
                    const somatic_clin_sig_key = this.eKeys.key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE);
                    const somatic_clin_sig_pretty = somatic_clin_sig_key.prettyValue(somatic_clin_sig_value.value);
                    clinSigText.push(somatic_clin_sig_pretty.val);
                }
                headers.push(
                    $('<div>', {class:'text-center my-1', text: clinSigText.join(" - ")})
                );

                /*
                let conditionRow = $('<div>', {class:'my-1'});
                if (v.resolved_condition) {
                    VCForm.format_condition(v.resolved_condition).appendTo(conditionRow);
                } else {
                    conditionRow.text((v.condition || {}).value);
                }
                */
                const flagRow = $('<div>', {class:'d-flex mt-2 align-items-center', style:'min-height:20px'}).appendTo(titleDom);
                $('<div>', {class:'flex-grow'}).appendTo(flagRow);

                if (!v.version_is_published) {
                    content = {
                        icon: '/static/icons/share_level/draft.png',
                        title: 'Unpublished working version'
                    };
                } else {
                    content = EKeys.shareLevelInfo(v.version_publish_level, v);
                    content.title = 'Shared with ' + content.title;
                }
                $('<img>', {
                    src: content.icon,
                    title: content.title,
                    width: '16px',
                    height: '16px',
                    display: 'inline-block'
                }).appendTo(flagRow);

                $('<div>', {'class': 'ml-2 text-center', 'data-flags': v.flag_collection, text: ''}).appendTo(flagRow);
                $('<div>', {class:'flex-grow'}).appendTo(flagRow);

                headers.push(flagRow);

                const th = $('<th>', {html: headers}).appendTo(headerRow);
            });

            const rowLabRecordId = $('<tr>', {class: 'group no-compare'}).appendTo(table);
            $('<th>', {text: "Lab ID", style: 'font-weight:normal'}).appendTo(rowLabRecordId);
            includedVersions.forEach(v => {
                $('<td>', {text: v.lab_record_id}).appendTo(rowLabRecordId);
            });

            let all_db_refs = {};

            // loop through groups of keys
            this.groups.forEach(group => {
                const row = $('<tr>', {class: 'group no-compare'}).appendTo(table);
                $('<th>', {text:group["label"]}).appendTo(row);
                $('<td>', {class:'filler', colspan: includedVersions.length}).appendTo(row);

                group.eKeys.forEach(eKey => {
                    const key = eKey.key;

                    includedVersions.forEach((v, index) => {
                        const db_refs = (v[key] || {}).db_refs;
                        if (db_refs) {
                            CitationsManager.normalizeInPlace(db_refs);

                            for (const db_ref of db_refs) {
                                let ref_details = all_db_refs[db_ref.id];

                                if (!ref_details) {
                                    ref_details = {dbRef: db_ref, whoHas: {}};
                                    all_db_refs[db_ref.id] = ref_details;
                                }
                                ref_details.source = 'reference';
                                ref_details.whoHas[index] = (ref_details.whoHas[index] || []).concat([eKey.label]);

                                if (Object.keys(ref_details.whoHas).length === includedVersions.length) {
                                    ref_details.everyoneHas = true;
                                }
                            }
                        }
                    });

                    for (const show of ['value', 'resolved', 'note']) {
                        let labelText = eKey.label;
                        if (show !== 'value') {
                            labelText += ` ${show}`;
                        }
                        const label = $('<span>', {text: labelText});

                        const row = $('<tr>', {class: `${show}`}).appendTo(table);

                        const th = $('<th>', {html: label}).appendTo(row);

                        const uniqueValues = {};
                        const blankValues = {};
                        let hasBlank = 0;
                        let hasMultiValues = false;
                        const isCriteria = show === 'value' && eKey.value_type === 'C' && eKey.namespace() === 'acmg';
                        if (eKey.value_type === 'T') {
                            row.addClass('textarea');
                        }
                        if (isCriteria) {
                            row.addClass('criteria');
                        }

                        includedVersions.forEach((v, index) => {

                            const blob = v[key] || {};
                            const {note, explain, resolved, value, hidden, db_refs} = blob;

                            const cell = $('<td>').appendTo(row);
                            let val = null;

                            if (hidden) {
                                val = null;
                                cell.text('hidden');
                                cell.addClass('hidden-value');
                            } else if (show === 'value') {
                                val = Diff.emptyToNull(value);
                                if (val !== null) {
                                    eKey.formatValue(val, cell, true);
                                    hasMultiValues = hasMultiValues || Array.isArray(value) && value.length > 1;
                                } else if (isCriteria) {
                                    cell.text('Not Set');
                                    cell.addClass('not-met');
                                }
                            } else if (show === 'note') {
                                val = Diff.emptyToNull(note);
                                cell.text(val);
                            } else if (show === 'resolved') {
                                if (resolved) {
                                    val = resolved.display_text;
                                    cell.html(VCTable.condition(resolved));
                                } else {
                                    val = null;
                                }
                            }

                            const valueKey = Diff.hash(val);
                            if (valueKey) {
                                cell.addClass(`not-blank`);
                                if (!uniqueValues[valueKey]) {
                                    uniqueValues[valueKey] = {value: `${eKey.prettyValue(val).val || ''}`, count: 1};
                                } else {
                                    uniqueValues[valueKey].count++;
                                }
                            } else {
                                if (isCriteria) {
                                    if (!blankValues[val]) {
                                        blankValues[val] = {value: `${eKey.prettyValue(val).val || 'Not Set'}`, count: 1};
                                    } else {
                                        blankValues[val].count++;
                                    }
                                }
                                cell.addClass(`blank`);
                                hasBlank++;
                            }

                            if (explain) {
                                // cell.addClass('explained');
                                cell.append($('<i>', {class:"fas fa-info-circle text-muted explain-icon"}));
                                cell.attr('title', 'Lab specific explanation');
                                cell.attr('data-content', explain);
                            }
                        });

                        let diffHelp = null;
                        const uniqueIndex = Object.keys(uniqueValues).length;
                        if (uniqueIndex === 0) {
                            diffHelp = 'No values for any cell';
                            row.addClass('empty');
                            // for crtieria blank vs just 1 other value can be considered significant
                            // for most other fields you can consider blank to just not be filled in or uploaded etc
                        } else if (uniqueIndex > 1 || (uniqueIndex === 1 && hasBlank && isCriteria)) {
                            row.addClass('differences');

                            const effectiveHasBlank = hasBlank;
                            if (isCriteria) {
                                const effectiveUniqueCount = uniqueIndex + Object.values(blankValues).length;
                                diffHelp = `<span style="color:#A44">Cells have ${effectiveUniqueCount} unique values</span>`;
                            } else {
                                diffHelp = `<span style="color:#A44">Cells have ${uniqueIndex} unique values` + (hasBlank ? ' and blanks ' : '') + '</span>';
                            }

                            th.prepend($('<i class="fas fa-not-equal" style="color:#b88"></i>'));
                        } else {
                            row.addClass('no-differences');
                            let icon;
                            if (hasBlank) {
                                row.addClass('blanks');
                                icon = $('<i class="fas fa-check-square" style="color:#bbb"></i>');
                                diffHelp = '<span style="color:#888">All non-empty cells have the same value</span>';
                            } else {
                                icon = $('<i class="fas fa-check-square" style="color:#6d6"></i>');
                                diffHelp = '<span style="color:#4A4">All cells have the same value</span>';
                            }
                            th.prepend(icon);
                        }

                        const diffBreakdown = [];
                        if (hasBlank) {
                            if (isCriteria) {
                                const blankCounts = Object.values(blankValues);
                                blankCounts.sort((vc1,vc2) => {
                                    if (vc1.count !== vc2.count) {
                                        return vc1.count - vc2.count;
                                    }
                                    return vc1.value.localeCompare(vc2.value);
                                });
                                for (const valueCount of blankCounts) {
                                    diffBreakdown.push(`${valueCount.count} <span style="color:#888;font-size:smaller">x</span> <span style="color:#888">${valueCount.value}</span>`);
                                }
                            } else {
                                diffBreakdown.push(`${hasBlank} <span style="color:#888;font-size:smaller">x</span> <span style="color:#888">blank</span>`);
                            }
                        }
                        const valueCounts = Object.values(uniqueValues);
                        valueCounts.sort((vc1,vc2) => {
                            if (vc1.count !== vc2.count) {
                                return vc1.count - vc2.count;
                            }
                            return vc1.value.localeCompare(vc2.value);
                        });
                        for (const valueCount of valueCounts) {
                            let displayValue = valueCount.value;
                            displayValue = displayValue.replace(/</g,'&lt;').replace(/\n/g,' ');
                            if (displayValue.length > 40) {
                                displayValue = displayValue.substring(0, 40) + '<span style="color:#888">...</span>';
                            }
                            diffBreakdown.push(`${valueCount.count} <span style="color:#888;font-size:smaller">x</span> ${displayValue}`);
                        }
                        const diffBreakdownString = diffBreakdown.join('<br/>');

                        let keyDescription;
                        if (eKey.description) {
                            keyDescription = EKeys.fixDescription(eKey.description).prop('outerHTML').replaceAll("\n", "<br/>");
                        } else {
                            keyDescription = "";
                        }
                        const content = $('<div>', {style:'white-space: pre-wrap;'});
                        th.attr('title', eKey.label);
                        th.attr('data-content', `${diffHelp}<br>${diffBreakdownString}<br>---<br>${keyDescription}`);

                        if (show === 'value' && includedVersions.length >= 2 && hasMultiValues) {
                            // break down multiselect or dbRefs into their own checkboxes
                            const multiselect = includedVersions.map(v => {
                                if (v[key] && Array.isArray(v[key].value)) {
                                    return v[key].value.map(k => ({id: k, label: eKey.prettyValue([k]).val}));
                                }
                                return null;
                            });
                            const dbRefs = includedVersions.map(v => {
                                if (v[key] && v[key].db_refs) {
                                    return CitationsManager.normalizeInPlace(v[key].db_refs);
                                }
                                return null;
                            });

                            if (multiselect.find(x => x !== null)) {
                                row.children('td').addClass('text-center');
                                renderDbRefCompare(this.dbRefCompare(multiselect, 'multiselect'));
                            }
                        }
                    }
                });
            });
            
            // convert from map to array and sort
            all_db_refs = Object.keys(all_db_refs).map(k => all_db_refs[k]);
            all_db_refs.sort((aw,bw) => {
                const [a,b] = [aw.dbRef, bw.dbRef];
                let diff = a.db.localeCompare(b.db);
                if (diff === 0) {
                    try {
                        diff = parseInt(a.idx) - parseInt(b.idx);
                    } catch (e) {}
                }
                if (diff === 0) {
                    diff = a.idx.localeCompare(b.idx);
                }
                return diff;
            });

            const references = all_db_refs.filter(ref => Diff.LITERATURE_DBS[ref.dbRef.db] !== true);
            const citations = all_db_refs.filter(ref => Diff.LITERATURE_DBS[ref.dbRef.db] === true);

            let row = $('<tr>', {class: 'group no-compare'}).appendTo(table);
            $('<th>', {class: 'heading', text:"Links"}).appendTo(row);

            renderDbRefCompare(references);

            row = $('<tr>', {class: 'group no-compare'}).appendTo(table);
            $('<th>', {class: 'heading', text:"Citations"}).appendTo(row);

            renderDbRefCompare(citations);

            e.empty().append(table);
            
            Flags.instance.init(this.flagParams);
            this.applyCriteriaHiding();
        }
    };

    return Diff;

})();