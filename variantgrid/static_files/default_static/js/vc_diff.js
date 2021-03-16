//
// DIFF
//
// Render a table of rows that are different by VCForm fields

const Diff = (function() {

    const Diff = function(records, diffDom, eKeys, flagParams) {
        
        this.flagParams = flagParams;
        this.versions = records.map(r => {
            let merged = Object.assign({}, r, r.data); 
            delete merged['data'];
            return merged;
        });
        let comparingDiffRecords = false;
        let firstId = null;
        for (let record of records) {
            if (!firstId) {
                firstId = record.id;
            } else if (record.id !== firstId) {
                comparingDiffRecords = true;
                break;
            }
        }

        let groupMap = {};
        let criteriaMap = {};
        for (let familyKey of Object.keys(EKey.families)) {
            groupMap[familyKey] = {"label": EKey.families[familyKey], eKeys: []};
            criteriaMap[familyKey] = {eKeys: []};
        }
        
        eKeys.forEach(k => {
            let key = k.key;
            let isCriteria = k.value_type === 'C';
            let category = k.evidence_category;

            let dict = (isCriteria ? criteriaMap : groupMap)[category];
            if (dict) {
                dict.eKeys.push(k);
            }
        });

        let groups = Object.keys(groupMap).map(g => {
            let combined = Object.assign({}, groupMap[g]);
            combined.eKeys = criteriaMap[g].eKeys.concat(combined.eKeys);
            return combined;
        }).filter(g => g.eKeys.length);

        this.diffDom = diffDom;
        this.eKeys = eKeys;
        this.groups = groups;

        this.showAllAcmg = comparingDiffRecords;

        this.excluded = {};
    };

    Diff.LITERATURE_DBS = {'PubMed': true, 'NCBIBookShelf': true};

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
            let sortDb = db1.db.localeCompare(db2.db);
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
            let modalContent = createModalShell('configureColumnsModal', 'Choose which columns to compare');
            let content = modalContent.find('.modal-body');

            let compareBox = $('<table>').appendTo(content);
            this.versions.forEach((v, index) => {
                let tr = $('<tr>').appendTo(compareBox);
                let ctd = $('<td>').appendTo(tr);
                let uniqueId = `${v.id}`;
                if (v.version) {
                    uniqueId += `.${v.version}`;
                }
                
                let checkbox = $('<input>', {type:'checkbox', id:`ch-${uniqueId}`, change:event => {
                    let checked = $(event.target).prop('checked');
                    this.excluded[uniqueId] = !checked;
                    this.render();
                    
                }}).appendTo(ctd);
                if (!this.excluded[uniqueId]) {
                    checkbox.prop('checked', true);
                }
                
                let ltd = $('<td>').appendTo(tr);
                let label = $('<label>', {for: `ch-${uniqueId}`, text: `${v.lab_name} / ${v.lab_record_id}`}).appendTo(ltd);
            });
            modalContent.on('');
            modalContent.modal({

            });
        },
               
        dbRefCompare(versionDbRefs, source) {
            // Move this out into its own function
            let allDbRefs = {};
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
            let orderedDbRefs = Object.keys(allDbRefs).map(k => allDbRefs[k]);
            return orderedDbRefs.sort((meta1, meta2) => Diff.compareDbRefs(meta1.dbRef, meta2.dbRef));
        },

        render: function() {
            let e = $(this.diffDom);
            let includedVersions = this.versions.filter(v => this.isIncluded(v));
           
            let table = $('<table>', {class: 'diff table'});
            let thead = $('<thead>').appendTo(table);
            let tbody = $('<tbody>').appendTo(table);
            let headerRow = $('<tr>', {class: 'versions no-compare'}).appendTo(thead);

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
                    let dbRef = dbRefC.dbRef;
                    let whoHas = dbRefC.whoHas;
                    let rowClass = dbRef.label ? 'option-row' : 'citation-row';
                    let row = $('<tr>', {class: `${rowClass}`}).appendTo(tbody);

                    let icon;
                    if (!dbRefC.everyoneHas) {
                        row.addClass('blanks');
                        icon = $('<i class="fas fa-not-equal" style="color:#b88"></i>');
                    } else {
                        icon = $('<i class="fas fa-check-square" style="color:#6d6"></i>');
                    }
                    let effectiveId = Citations.effectiveId(dbRef);
                    let isDetailed = !!effectiveId;
                    let th;
                    if (dbRef.label) {
                        th = $('<th>', {html: [icon, dbRef.label]}).appendTo(row);
                    } else {
                        th = $('<th>', {class: isDetailed ? `citation-row` : `simple-citation`, html: [icon, Citations.renderDbRef(dbRef)]}).appendTo(row);
                    }
                    let hasCount = 0;
                    let hasnotCount = 0;
                    includedVersions.filter(v => this.isIncluded(v)).forEach((v, index) => {
                        let uniqueId = `${v.id}`;
                        if (v.version) {
                            uniqueId += `.${v.version}`;
                        }
                        if (whoHas[index]) {
                            hasCount++;
                            let title = null;
                            let whoHasData = whoHas[index];
                            if (Array.isArray(whoHasData)) {
                                title = 'Referenced in ' + whoHasData.join(', ');
                            }
                            $('<td>', {
                                title: title,
                                html: $('<i class="far fa-check-square"></i>')
                            }).appendTo(row);
                           //$('<td>', {text: 'included', class: primeInclusion ? '' : 'addition'}).appendTo(row);
                       } else {
                            hasnotCount++;
                            $('<td>', {html: '<i class="far fa-square" style="opacity:0.2""></i>'}).appendTo(row);
                        }
                    });
                    let tooltip = [];
                    let referenceType = dbRefC.source === 'multiselect' ? 'value' : 'reference';
                    if (hasCount) {
                        let verb = hasCount === 1 ? 'record has this' : 'records have this';
                        tooltip.push(`${hasCount} <span style="color:#444;font-size:smaller">x</span> ${verb} ${referenceType}`);
                    }
                    if (hasnotCount) {
                        let verb = hasnotCount === 1 ? 'record does not have this' : 'records do not have this';
                        tooltip.push(`${hasnotCount} <span style="color:#444;font-size:smaller">x</span> ${verb} ${referenceType}`);
                    }
                    th.attr('title', dbRef.label || dbRef.id);
                    th.attr('data-content', tooltip.join('<br>'));
                });
            };
            renderDbRefCompare = renderDbRefCompare.bind(this);
            new Citations().refresh();

            //var percent = Math.floor(100 / Math.max(1, this.versions.length));
            includedVersions.forEach((v, index) => {
                let url = Urls.view_classification(v.id);
                if (v.version) {
                    url += `.${v.version}`;
                }

                let titleDom = $('<a>', {class:'hover-link text-center d-flex flex-column flex-align-center', href: url});
                
                let titlePart = $('<div>', {text: v.lab_name + ' / ' + v.lab_record_id}).appendTo(titleDom);
                if (v.first_seen) {
                    let first_seen_date = moment(v.first_seen * 1000).format('DD/MMM/YYYY HH:mm');
                    let last_seen_date = moment(v.version * 1000).format('DD/MMM/YYYY HH:mm');
                    if (first_seen_date !== last_seen_date) {
                        let title = `First uploaded at ${first_seen_date},<br/>Non-significant change at ${last_seen_date}`;
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
                let flagRow = $('<div>', {class:'d-flex mt-2 align-items-center', style:'min-height:20px'}).appendTo(titleDom);
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
                let clin_sig_key = this.eKeys.key(SpecialEKeys.CLINICAL_SIGNIFICANCE);
                let clin_sig = clin_sig_key.prettyValue((v.clinical_significance || {}).value);
                let clinSigRow = $('<div>', {class:'text-center my-1', text:clin_sig.val});
                /*
                let conditionRow = $('<div>', {class:'my-1'});
                if (v.resolved_condition) {
                    VCForm.format_condition(v.resolved_condition).appendTo(conditionRow);
                } else {
                    conditionRow.text((v.condition || {}).value);
                }
                */
                $('<div>', {'class': 'ml-2 text-center', 'data-flags': v.flag_collection, text: ''}).appendTo(flagRow);
                $('<div>', {class:'flex-grow'}).appendTo(flagRow);

                let th = $('<th>', {html: [
                        titleDom,
                        clinSigRow,
                        flagRow,
                        // conditionRow
                    ]}).appendTo(headerRow);
            });
            let all_db_refs = {};

            // loop through groups of keys
            this.groups.forEach(group => {
                let row = $('<tr>', {class: 'group no-compare'}).appendTo(table);
                $('<th>', {text:group["label"]}).appendTo(row);
                $('<td>', {class:'filler', colspan: includedVersions.length}).appendTo(row);

                group.eKeys.forEach(eKey => {
                    let key = eKey.key;

                    includedVersions.forEach((v, index) => {
                        let db_refs = (v[key] || {}).db_refs;
                        if (db_refs) {
                            for (let db_ref of db_refs) {
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

                    for (let show of ['value', 'note', 'processed']) {
                        let labelText = eKey.label;
                        if (show !== 'value') {
                            labelText += ` ${show}`;
                        }
                        let label = $('<span>', {text: labelText});

                        let row = $('<tr>', {class: `${show}`}).appendTo(table);

                        let th = $('<th>', {html: label}).appendTo(row);

                        let uniqueValues = {};
                        let blankValues = {};
                        let hasBlank = 0;
                        let hasMultiValues = false;
                        let isCriteria = show === 'value' && eKey.value_type === 'C' && eKey.namespace() === null;
                        if (eKey.value_type === 'T') {
                            row.addClass('textarea');
                        }
                        if (isCriteria) {
                            row.addClass('criteria');
                        }

                        includedVersions.forEach((v, index) => {

                            let blob = v[key] || {};
                            let {note, explain, processed, value, hidden, db_refs} = blob;

                            let cell = $('<td>').appendTo(row);
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
                            } else if (show === 'processed') {
                                val = Diff.emptyToNull(processed);
                                cell.text(val);
                            }

                            let valueKey = Diff.hash(val);
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
                        let uniqueIndex = Object.keys(uniqueValues).length;
                        if (uniqueIndex === 0) {
                            diffHelp = 'No values for any cell';
                            row.addClass('empty');
                            // for crtieria blank vs just 1 other value can be considered significant
                            // for most other fields you can consider blank to just not be filled in or uploaded etc
                        } else if (uniqueIndex > 1 || (uniqueIndex === 1 && hasBlank && isCriteria)) {
                            row.addClass('differences');

                            let effectiveHasBlank = hasBlank;
                            if (isCriteria) {
                                let effectiveUniqueCount = uniqueIndex + Object.values(blankValues).length;
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

                        let diffBreakdown = [];
                        if (hasBlank) {
                            if (isCriteria) {
                                let blankCounts = Object.values(blankValues);
                                blankCounts.sort((vc1,vc2) => {
                                    if (vc1.count !== vc2.count) {
                                        return vc1.count - vc2.count;
                                    }
                                    return vc1.value.localeCompare(vc2.value);
                                });
                                for (let valueCount of blankCounts) {
                                    diffBreakdown.push(`${valueCount.count} <span style="color:#888;font-size:smaller">x</span> <span style="color:#888">${valueCount.value}</span>`);
                                }
                            } else {
                                diffBreakdown.push(`${hasBlank} <span style="color:#888;font-size:smaller">x</span> <span style="color:#888">blank</span>`);
                            }
                        }
                        valueCounts = Object.values(uniqueValues);
                        valueCounts.sort((vc1,vc2) => {
                            if (vc1.count !== vc2.count) {
                                return vc1.count - vc2.count;
                            }
                            return vc1.value.localeCompare(vc2.value);
                        });
                        for (let valueCount of valueCounts) {
                            let displayValue = valueCount.value;
                            displayValue = displayValue.replace(/</g,'&lt;').replace(/\n/g,' ');
                            if (displayValue.length > 40) {
                                displayValue = displayValue.substring(0, 40) + '<span style="color:#888">...</span>';
                            }
                            diffBreakdown.push(`${valueCount.count} <span style="color:#888;font-size:smaller">x</span> ${displayValue}`);
                        }
                        let diffBreakdownString = diffBreakdown.join('<br/>');

                        let keyDescription = eKey.description || '';
                        th.attr('title', eKey.label);
                        th.attr('data-content', `${diffHelp}<br>${diffBreakdownString}<br>---<br>${keyDescription}`);

                        if (show === 'value' && includedVersions.length >= 2 && hasMultiValues) {
                            // break down multiselect or dbRefs into their own checkboxes
                            let multiselect = includedVersions.map(v => {
                                if (v[key] && Array.isArray(v[key].value)) {
                                    return v[key].value.map(k => ({id: k, label: eKey.prettyValue([k]).val}));
                                }
                                return null;
                            });
                            let dbRefs = includedVersions.map(v => {
                                if (v[key] && v[key].db_refs) {
                                    return v[key].db_refs;
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
                let [a,b] = [aw.dbRef, bw.dbRef];
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

            let references = all_db_refs.filter(ref => Diff.LITERATURE_DBS[ref.dbRef.db] !== true);
            let citations = all_db_refs.filter(ref => Diff.LITERATURE_DBS[ref.dbRef.db] === true);

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