/* jshint esversion: 6 */
//
// VCForm
//
// Render vcform form
const VCForm = (function() {

    // need a wrapper as disabled elements don't have tool tips
    const disableButton = button => {
        let wrapper = $('<div>', {class: 'btn-disabled-wrapper', style:'position:relative; cursor: not-allowed'});
        button.attr('click', null);
        button.attr("disabled", true);
        button.addClass('disabled');
        button.addClass('h-100');

        wrapper.append(button);
        let bonusWrapper = $('<div>', {style:'z-index:1; position:absolute; left:0; top:0; right:0; bottom:0', title: button.attr('title'), 'data-toggle':"tooltip"});
        bonusWrapper.tooltip({html:true, focus:'hover'});
        wrapper.append(bonusWrapper);

        return wrapper;
    };

    /* private methods */
    function emptyToNull(val) {
        if (typeof(val) === 'string') {
            val = val.trim();
            if (val.length === 0) {
                return null;
            }
            return val;
        } else if (typeof(val) === 'undefined') {
            return null;
        }
        return val;
    }
    
    /* private variables */
    let filtering = false;
    let filtered = {};
    
    let jContent = null;
    let jSyncStatus = null;
    let jHelp = null;
    let jCritTable = null;
    let jFilterBox = null;
    let jClearFilterButton = null;
    let jErrors = null;
    let jPublishHistory = null;
    let jLinks = null;
    let jShareButtons = null;

    let eKeys = null;
    let vcLinks = null;

    let VCForm = function() {};

    VCForm.prototype = {

        url: null,
        record: null,
        data : {},
        delta: {},
        publish_level: null,
        deleteRequest: false,
        delete_reason: null,
        undeleteRequest: false,
        messages: [],
        delayedPatch: {},

        renderReference(ref) {
            let text = ref.id;
            if (ref.db === "HTTP" || ref.db === "HTTPS") {
                text = `${ref.db.toLowerCase()}:${ref.idx}`;
            }

            return $('<div>', {
                class: 'ref',
                'data-preview-db': ref.db,
                'data-preview-id': ref.idx,
                html: [
                    $('<a>', {
                        class: 'hover-link external-link',
                        href: ref.url,
                        text: text,
                        target: '_blank',
                    }),
                    $('<div>', {class: 'ref-summary ml-1 my-1 d-inline', text: '...'})
                ]
            });
        },

        isEditMode() {
            return this.record.can_write && window.location.toString().indexOf("edit=true") !== -1;
        },

        editMode() {
            let url;
            if (document.location.toString().indexOf('?') !== -1) {
                url = document.location.toString() + "&edit=true";
            } else {
                url = document.location.toString() + "?edit=true";
            }
            window.open(url, '_self');
        },

        clear() {
            let res = confirm('Are you sure you wish to clear the form data?');
            if (res) {
                Object.keys(this.data).forEach(key => {
                    this.data[key] = null;
                });
                this.clearFilter();
                this.populateForm();
            }
            this.dataUpdated();
        },

        baseUrl() {
            let url = window.location.href;
            let questionIndex = url.indexOf("?");
            if (questionIndex !== -1) {
                return url.substring(0, questionIndex);
            }
            return url;
        },

        csv() {
            window.open(this.baseUrl() + '/classification.csv');
            return false;
        },
        
        report() {
            window.open(this.baseUrl() + '/report.html', '_blank');
            return false;
        },
        
        generateLink(href, html, new_tab) {
            let a = $('<a>', {href: href, html: html });
            if (new_tab) {
                a.attr('target', '_blank');
            }
            return a;
        },

        // init data with a record (or from localStorage if no record is there)
        initData(record) {
            let data = {};
            if (record === null) {
                // record = JSON.parse( localStorage.editingSubmission || '{}' );
                record = {};
            }
            this.record = record;
            let recordData = Object.assign({}, record['data'] || {});
            let recordMessages = record['messages'] || [];
            
            eKeys.forEach(eKey => {
                let key = eKey.key;
                data[key] = recordData[key] || null;
                delete recordData[key]; 
            });
            // only unknown keys are left
            Object.keys(recordData).forEach(unknownKey => {
               eKeys.key(unknownKey); // causes the key to get created
               data[unknownKey] = recordData[unknownKey]; 
            });
            
            this.messages = recordMessages;
            this.data = data;
            this.updateLinks();
            this.updateErrors();
            this.updateCitations();
            this.updateTitle();
            this.updatePublishHistory();
            
            if (!this.isEditMode()) {
                this.searchAsterisk();
            } else {
                /*
                window.setTimeout(() => {
                    $('#vc-form .collapse').first().collapse('show');
                });
                */
            }
        },
        
        labelForLevel(share_level) {
            return EKeys.shareLevelInfo(share_level, this.record, false);
        },
                
        generateExportButtons() {
            let wrapper = $('<div>', {html: $('<h5>', {text:'Export as', class: 'mt-4'})});
            let buttons = $('<div>', {class: 'btn-toolbar'}).appendTo(wrapper);
            let csvButton = $('<button>', {class:'btn btn-outline-primary btn-lg', id: 'export-csv', html: '<i class="fas fa-file-csv"></i> CSV', click: () => {this.csv()}});
            csvButton.appendTo(buttons);

            if (this.reportEnabled) {
                let reportButton = $('<button>', {class:'btn btn-outline-primary btn-lg', id: 'export-report', html: '<i class="far fa-file"></i> Report', click: () => {this.report()}});
                reportButton.appendTo(buttons);
            } else {
                buttons.append($('<div>'));
            }
            
            return wrapper;
        },
        createModal(headerText, bodyText, confirmButtonText, confirmCallback, cancelCallback, isWithdrawal = false) {
            let dialogContent = createModalShell('confirmation-modal', headerText, 'lg');
            let dialogBody = dialogContent.find('.modal-body');
            let dialogFooter = dialogContent.find('.modal-footer');
            let idPrefix = "modal-action";

            if (isWithdrawal) {
                let withdrawReasonDropdown = $('<select class="form-control" id="withdrawReasonDropdown" required></select>');
                    withdrawReasonDropdown.append('<option value="">Select Withdrawal Reason</option>');

                    for (const [value, display] of this.withdrawReasons) {
                        withdrawReasonDropdown.append(`<option value="${value}">${display}</option>`);
                    }

                dialogBody.append(bodyText, '<br>', withdrawReasonDropdown);
                    if (!this.deleteEnabled){
                        dialogBody.append('<p style="padding:16px 4px 0 4px"><i class="fas fa-exclamation-triangle text-warning"></i> Withdrawing a record is intended for a limited set of circumstances, ' +
                        'listed as available reasons for withdrawal above.<br/><br/>If this record was uploaded from your system, ' +
                        'and you believe the data to be incorrect or incomplete, please note this is considered a ' +
                        'normal circumstance, and it is preferred that you update the classification within your ' +
                        'curation system and then resubmit the classification in your next upload, rather than withdraw it.</p>');
                    }

            } else {
                dialogBody.html(bodyText);
            }

            dialogFooter.html([
                $('<button>', {
                    type: "button",
                    class: "btn btn-secondary",
                    id: `${idPrefix}-cancel`,
                    'data-dismiss': "modal",
                    text: 'Cancel',
                    click: cancelCallback
                }),
                $('<button>', {
                    type: "button",
                    class: "btn btn-primary",
                    text: confirmButtonText,
                    id: `${idPrefix}-confirm`,
                    click: function () {
                        if (isWithdrawal && !$('#withdrawReasonDropdown').val()) {
                            alert('Please select a withdrawal reason.');
                            return;
                        }
                        dialogBody.LoadingOverlay('show', {zIndex: 9999});
                        confirmCallback();
                    }
                })
            ]);

            dialogContent.modal();

            dialogContent.on('hidden.bs.modal', function () {
                dialogContent.modal('dispose');
                dialogContent.remove();
            });
        },

        trash() {
            let confirmAction, cancelAction;
            if (this.record.withdrawn) {
                confirmAction = () => {
                    this.delete_reason = null;
                    this.undeleteRequest = true;
                    this.dataUpdated();
                };
                this.createModal(
                    'Un-withdraw Record',
                    'Are you sure you wish to un-withdraw this record?',
                    'Un-withdraw',
                    confirmAction
                );
            } else {
                let is_withdraw = this.record.publish_level === 'logged_in_users' ||
                    this.record.publish_level === 'public' || !this.deleteEnabled;
                let mode = is_withdraw ? "withdraw" : "delete";
                confirmAction = () => {
                    if (is_withdraw) {
                        this.delete_reason = $('#withdrawReasonDropdown').val();
                        this.deleteRequest = "withdraw";
                    } else {
                        this.deleteRequest = true;
                    }
                    this.dataUpdated();
                };

                cancelAction = () => {
                };

                this.createModal(
                    `${mode.charAt(0).toUpperCase() + mode.slice(1)} Record`,
                    `Are you sure you wish to ${mode} this record?`,
                    mode.charAt(0).toUpperCase() + mode.slice(1),
                    confirmAction,
                    cancelAction,
                    is_withdraw
                );
            }
        },
        
        submit() {
            this.publish(this.record.publish_level);
        },
        
        share() {
            let dialogContent = createModalShell('share', 'Submit and Share This With');
            let dialogBody = dialogContent.find('.modal-body');
            let dialogFooter = dialogContent.find('.modal-footer');
            dialogContent.on('hidden.bs.modal', function() {
                dialogContent.modal('dispose');
                dialogContent.remove();
            });

            let shareLevelDivs = [];
            let newShareLevel = this.record.publish_level;
            for (let share_level of VcSettings.SHARE_LEVELS) {
                let label = this.labelForLevel(share_level);
                
                let selectionDiv = $('<li>', {class: 'list-group-item share-list', 'data-current': label.included || label.current}).appendTo(dialogBody);
                shareLevelDivs.push(selectionDiv);
                $('<i class="far fa-circle d-inline-block"></i>').appendTo(selectionDiv);
                $('<img>', {class: 'ml-2 d-inline-block share-icon', width: '16px', height: '16px', src: label.icon }).appendTo(selectionDiv);
                $('<div>', {class: 'ml-2 d-inline', text: label.title}).appendTo(selectionDiv);
                
                if (label.included) {
                    selectionDiv.addClass('included list-group-item-primary');
                    selectionDiv.find('i').removeClass('fa-circle').addClass('fa-check-circle');
                    selectionDiv.attr('title', 'This classification record has already been shared at a higher level');
                } else {
                    selectionDiv.addClass('list-group-item-action');
                    if (label.current) {
                        selectionDiv.find('i').removeClass('fa-circle').addClass('fa-check-circle');
                        selectionDiv.addClass('list-group-item-primary');
                    }
                    selectionDiv.click(() => {
                        newShareLevel = share_level;
                        let hitMeYet = false;
                        for (let divy of shareLevelDivs) {
                            if (hitMeYet) {
                                divy.removeClass('list-group-item-secondary');
                                divy.find('i').removeClass('fa-check-circle').addClass('fa-circle');
                            } else if (divy.attr('data-current') !== 'true') {
                                divy.addClass('list-group-item-secondary');
                                divy.find('i').removeClass('fa-circle').addClass('fa-check-circle');
                            }

                            if (divy === selectionDiv) {
                                hitMeYet = true;
                            }
                        }
                    });
                }
            }
            $('<p>', {class: 'mt-2', text: 'Note that once shared at a certain level, this classification record can only be shared at the same or higher level'}).appendTo(dialogBody);

            dialogFooter.html([
                $('<button>', {type:"button", class:"btn btn-secondary", 'data-dismiss':"modal", text:'Cancel', id:'share-cancel'}),
                $('<button>', {type:"button", class:"btn btn-primary", text:'Submit', id:'share-confirm', click: () => {
                    this.publish(newShareLevel);
                    dialogContent.modal('hide');
                }})
            ]);
            dialogContent.modal();
        },

        updateSubmitButton() {
            // The submit button appears next to the quick summary now so it's always at hand
            let quickSubmitWrapper = $('#vc-quick-submit');
            quickSubmitWrapper.empty();
            if (quickSubmitWrapper.length) {
                let message = null;
                let provideSubmitButton = true;

                if (this.record.can_write) {
                    if (this.record.withdrawn) {
                        message = "You must unwithdraw to perform any changes.";
                    //} else if (this.hasSendingStatus()) {
                        // message = "Changes uploading.";
                        // shouldn't happen as when we have a sending status it just says Sending
                    } else if (this.hasErrors()) {
                        message = "Messages must be resolved before submitting.";
                        provideSubmitButton = false;
                    } else {
                        if (this.record.has_changes) {
                            message = "This record has unsubmitted changes.<br/>Changes won't be visible to other users.";
                        } else if (this.record.publish_level === 'lab') {
                            message = `This record is only visible to users in your lab group.`;
                        } else if (this.record.publish_level === 'institution') {
                            message = `This record is only visible to users in your organisation.`;
                        } else {
                            let message_suffix = VcSettings.LOGGED_IN_USERS_MESSAGE || "This record is shared to Shariant users.";
                            $('<div>', {class: 'text-center mt-3 mb-2', style:'font-size:14px', html: `<i class="fas fa-check-circle text-success"></i> ${message_suffix}`}).appendTo(quickSubmitWrapper);
                        }
                    }
                    if (message) {
                        $('<div>', {class: 'text-center mt-3 mb-2 font-weight-bold text-danger', style:'font-size:14px', html: '<i class="fas fa-exclamation-triangle text-warning"></i>' + message}).appendTo(quickSubmitWrapper);
                    }

                    if (!this.isEditMode()) {
                        $('<button>', {class:"mt-1 btn btn-warning w-100", html:`<i class=\"fas fa-unlock-alt\"></i> EDIT`, click:this.editMode}).appendTo(quickSubmitWrapper);
                    } else {
                        if (this.sendStatus.error) {
                            $('<div>', {
                                class: 'mt-1 btn btn-danger w-100',
                                disabled: true,
                                html: '<i class="fas fa-bomb"></i> Error. Please Reload Page.'
                            }).appendTo(quickSubmitWrapper);

                        } else if (this.hasSendingStatus()) {
                            $('<div>', {
                                class: 'mt-1 btn btn-secondary w-100',
                                disabled: true,
                                html: '<i class="fas fa-clock"></i> Uploading Changes'
                            }).appendTo(quickSubmitWrapper);
                        } else {
                            if (provideSubmitButton) {
                                $('<button>', {
                                    class: 'mt-1 btn btn-primary w-100',
                                    html: '<i class="fas fa-upload"></i> Submit',
                                    title: 'Submit/Share',
                                    id :'action-submit',
                                    click: this.share.bind(this)
                                }).appendTo(quickSubmitWrapper);
                            } else if (this.hasErrors()) {
                                $('<div>', {
                                    class: 'mt-2 border border-information rounded w-100 text-center p-2 text-secondary bg-light hover-detail',
                                    title: `Errors in the messages section must be fixed before you can submit.`,
                                    'data-toggle': 'hover',
                                    html: '<i class="fas fa-exclamation-circle"></i> Resolve Messages to Submit'
                                }).appendTo(quickSubmitWrapper);
                            }
                        }
                    }
                } else if (!this.record.is_last_published || this.record.can_write_latest) {
                    let goToLatest = () => {
                        window.open(`/classification/classification/${this.record.id}`, '_self');
                    };

                    $('<button>', {
                        class: 'mt-1 btn btn-primary w-100',
                        html: '<i class="fas fa-stopwatch"></i> Go to Latest Version',
                        title: 'Submit/Share',
                        id: 'action-latest',
                        click: goToLatest
                    }).appendTo(quickSubmitWrapper);
                } else {
                    $('<div>', {
                        class: 'mt-2 border border-information rounded w-100 text-center p-2 text-secondary bg-light hover-detail',
                        title: `Only users belonging to "${this.record.lab_name}" can edit this record.`,
                        'data-toggle': 'hover',
                        html: '<i class="fas fa-eye"></i> Read Only'
                    }).appendTo(quickSubmitWrapper);
                }
            }
        },

        generateActionButtons() {
            this.updateSubmitButton();

            let wrapper = $('<div>', {class:'mt-4'});
            wrapper.append($('<h5>', {text: 'Actions'}));
            let buttons = $('<div>', {class: 'btn-toolbar'}).appendTo(wrapper);

            let butt = $('<button>', {class: 'btn btn-danger btn-lg', html: '<i class="fas fa-trash-alt"></i> Delete', title: 'Delete', click: this.trash.bind(this)});
            if (this.record.withdrawn) {
                butt.attr('id', 'action-un-withdraw');
                butt.attr('title', 'Un-withdraw Classification');
                butt.html('<i class="fas fa-trash-restore-alt"></i> Unwithdraw');
                butt.removeClass('btn-danger');
                butt.addClass('btn-secondary');
            } else if (this.record.publish_level === 'logged_in_users' ||
                this.record.publish_level === 'public' ||
                !this.deleteEnabled) {
                //butt.attr('title', 'Cannot delete classification after it has been shared with logged in users or 3rd parties');
                //butt = disableButton(butt);
                butt.attr('id', 'action-withdraw');
                butt.attr('title', 'Withdraw Classification');
                butt.html('<i class="fas fa-trash-alt"></i> Withdraw');
            }
            if (this.hasSendingStatus()) {
                butt.attr('title', 'Saving... Please wait');
                butt = disableButton(butt);
            }

            butt.appendTo(buttons);

            return wrapper;
        },
        
        updateLinks() {
            jLinks.empty();
            let linkData = {};
            for (let key of VCLinks.ALL_KEYS) {
                linkData[key] = this.value(key);
            }

            linkData[SpecialEKeys.VARIANT_COORDINATE] = this.variantCoordinate();
            linkData[SpecialEKeys.C_HGVS] = this.cHGVS();

            let allLinks = vcLinks.generateLinks(linkData).map(vcLink => {return vcLink.asAnchor("bootstrap").addClass('list-group-item').addClass('list-group-item-action')});
            for (let link of allLinks) {
                link.attr('data-placement', 'left');
            }
            let length = allLinks.length;
            let colA = allLinks;
            let colB = colA.splice(0, Math.ceil(length/2));

            let columnLinks = [colB, colA];
            let row = $('<div>', {class:'row no-gutters'}).appendTo(jLinks);
            for (let i=0; i < 2; i++) {
                let col = $('<div>', {class: 'col-6', html:
                    $('<ul>', {class: 'list-group', html: columnLinks[i]})
                }).appendTo(row);
                if (i === 1) {
                    col.css({'border-left': '1px solid lightgray'});
                }
            }
        },
        
        updatePublishHistory() {
            jPublishHistory.show();
            jPublishHistory.empty();
            
            let publishUL = $('<ul>', {class: 'list-group'}).appendTo(jPublishHistory);
            
            let versions = [];
            if (this.record.version !== this.record.last_edited && this.record.version !== this.record.published_version) {
                versions.push({
                    label: 'this version',
                    icon: 'current',
                    timestamp: this.record.version
                });
            }
            
            if (this.record.last_edited) {
                versions.push({
                    label: 'last edited',
                    icon: 'draft',
                    timestamp: this.record.last_edited,
                    editable: true,
                    title: 'The time this record was last edited'
                });
            }
            
            versions.push({
                label: 'last shared ver.',
                icon: this.record.publish_level,
                timestamp: this.record.published_version,
                title: 'Version that was last shared with ' + this.labelForLevel(this.record.publish_level).title
            });
            
            for (let version of versions) {
                let content = $('<span>', {html:[
                    $('<img>', {src: '/static/icons/share_level/' + version.icon + '.png', class:'share-icon'}),
                    version.label
                ]});
                let isCurrentVersion = version.timestamp === this.record.version;

                if (!version.timestamp) {
                    $('<span>', {class: 'timestamp', text: '---' }).appendTo(content);
                    $('<li>', {class: 'list-group-item', title: 'this record has not been shared at any level yet', 'data-toggle': 'tooltip'}).appendTo(publishUL);
                } else {
                    let href = `/classification/classification/${this.record.id}`;
                    if (!version.editable) {
                        href += `.${version.timestamp}`;
                    }
                    $('<span>', {class: 'timestamp', text: ' ' + moment(version.timestamp * 1000).format('DD/MMM/YYYY HH:mm') }).appendTo(content);
                    if (this.baseUrl().endsWith(href)) {
                        $('<li>', {class: 'list-group-item list-group-item-success', html: content, title: version.title, 'data-toggle': 'tooltip'}).appendTo(publishUL);
                    } else {
                        $('<a>', {class: `list-group-item list-group-item-action ${isCurrentVersion ? '' : ''}`, title: version.title, html: content, 'data-toggle': 'tooltip', href: href }).appendTo(publishUL);
                    }
                }
            }

            if (this.record.can_write_latest) {
                this.generateLink(
                    `/classification/classification/${this.record.id}/history`,
                    '<i class="fas fa-history"></i> View change log', true
                ).addClass('list-group-item').addClass('list-group-item-action').appendTo(publishUL);
            }
            if (this.record.published_version) {
                this.generateLink(
                    `/classification/diff/?history=${this.record.id}&classification_id=${this.record.id}`,
                    '<i class="fas fa-history"></i> Compare with historical versions of this record', true
                ).addClass('list-group-item').addClass('list-group-item-action').appendTo(publishUL);
            }
            if (this.otherClassificationsSummary) {
                let viewSummary = "<i class=\"fas fa-columns\"></i> Compare with other classifications for this variant (" + this.otherClassificationsSummary + ")";
                this.generateLink(
                    `/classification/diff/?variant_compare=${this.record.id}&allele_id=${(this.record.allele || {}).id}`,
                    viewSummary, true
                ).addClass('list-group-item').addClass('list-group-item-action').appendTo(publishUL);
            }

            if (this.isEditMode()) {
                jShareButtons.find('[title]').tooltip('hide').tooltip('dispose');
                jShareButtons.empty();
                
                this.generateActionButtons().appendTo(jShareButtons);
            } else {
                this.updateSubmitButton();
            }
            this.generateExportButtons().appendTo(jShareButtons);
        },
        
        updateCitations() {
            let seenIds = {};
            let elements = [];
            Object.keys(this.data).forEach(k => {
               let blob = this.data[k];
               if (blob && blob['db_refs']) {
                   for (let ref of blob['db_refs']) {
                       let id = (ref.id || '').replace('\s*', '');
                       if (!seenIds[id]) {
                           let dom = CitationsManager.defaultManager.citationDomFor(ref);
                           if (dom) {
                               elements.push(dom);
                           }
                           seenIds[id] = true;
                       }
                   }
               } 
            });
            this.citations.empty().html(elements);
        },

        iconForSeverity(severity) {
            switch (severity) {
                case 'error':
                    return $('<i class="fas fa-exclamation-circle text-danger"></i>');
                case 'warning':
                    return $('<i class="fas fa-exclamation-triangle text-warning"></i>');
                case 'info':
                    return $('<i class="fas fa-info-circle text-primary"></i>');
            }
        },

        updateErrors() {
            $('.entry .inline-error').remove();
            if (this.messages.length === 0) {
                jErrors.closest('.card').hide();
                return;
            }
            this.messages.sort((a, b) => eKeys.key(a.key).label.localeCompare(eKeys.key(b.key).label));
            
            jErrors.closest('.card').show();
        
            let ul = $('<ul>', {class:'list-group'});
            this.messages.forEach(error => {
                let eKey = eKeys.key(error.key);
                
                // add the warning next to the field as well
                let inlineError = this.iconForSeverity(error.severity);
                inlineError.addClass('inline-error');
                inlineError.attr('title', error.message);

                $(`#label-${error.key}`).prepend(inlineError);

                switch (error.severity) {
                    case 'error': icon = '<i class="fas fa-exclamation-circle text-danger"></i>'; break;
                    case 'warning': icon = '<i class="fas fa-exclamation-triangle text-warning"></i>'; break;
                    case 'info': icon = '<i class="fas fa-info-circle text-primary"></i>'; break;
                }

                let listItem = $('<a>', {class: 'list-group-item list-group-item-action',  target: '_blank', click: () => {
                    jFilterBox.val('#' + eKey.key);
                    jFilterBox.keyup();
                    $(`[entry="${eKey.key}"]`).trigger('click');
                    // return false;
                }, html: [
                    this.iconForSeverity(error.severity),
                    $('<span>', {class: 'font-weight-bold', text: eKey.label}),
                    $('<span>', {class: 'message', text: ` - ${error.message}`})
                ]});
                ul.append(listItem);

                if (error.link) {
                    listItem.append($('<a>', {class: 'hover-link', href: error.link, target: '_blank', text: 'Click here for more details'}));
                }
            });
            jErrors.empty().append(ul);
        },
        
        hasErrors() {
            return this.messages.find(m => m.severity === 'error');
        },
        
        severity(key) {
            if (typeof(key) !== 'string') {
                key = key.key;
            }
            let messages = this.messages.filter(message => message.key === key);
            for (let severity of ['error', 'warning', 'info']) {
                if (messages.find(m => m.severity === severity)) { return severity; }
            }
            return null;
        },
        
        sendStatus: {
            modifications: false,
            sending: false,
            error: null
        },

        hasSendingStatus() {
            return this.sendStatus.modifications || this.sendStatus.sending || this.sendStatus.error;
        },
        
        setupUnloadWarning() {
            window.addEventListener("beforeunload", (e) => {
                if (!this.sendStatus.modifications && !this.sendStatus.sending) {
                   return undefined;
                }
                // Note custom message has almost no browser support
                let confirmationMessage = `
                You have unsaved changes (they will automatically be uploaded shortly).\n
                If you leave while changes are being uploaded, they will be lost.`;
                
                if (this.deleted || this.withdrawn) {
                    return null;
                } else {
                    (e || window.event).returnValue = confirmationMessage; //Gecko + IE
                    return confirmationMessage; //Gecko + Webkit, Safari, Chrome etc.
                }
            });
        },
        
        updateSendingStatus(delta) {
            Object.assign(this.sendStatus, delta);
            this.updateTitle();
            this.updateSubmitButton();
            this.updatePublishHistory();
        },

        publish(publish_level) {
            this.publish_level = publish_level;
            this.dataUpdated();
        },

        sendChanges: debounce(function() {
            if (this.sendStatus.error) {
                return;
            }
            let sendingDelta = Object.assign({}, this.delta);
            this.updateSendingStatus({sending: true});
            this.delta = {};

            // have to request config if we change a key responsible for which namespaces get used
            let requestConfig = !!sendingDelta[SpecialEKeys.ASSERTION_METHOD] || !!sendingDelta[SpecialEKeys.ALLELE_ORIGIN];
            let envelope = {
                'id': {'record_id': this.record.id},
                source: 'form',
                patch: sendingDelta,
                return_data: 'changes',
                config: requestConfig
            };
            if (this.publish_level) {
                envelope['publish'] = this.publish_level;
            }
            if (this.deleteRequest) {
                envelope['delete'] = this.deleteRequest;
                envelope['delete_reason'] = this.delete_reason;
            } else if (this.undeleteRequest) {
                envelope['delete'] = false;
                envelope['delete_reason'] = null;
            }
            
            $.ajax({
                headers: {
                    'Accept' : 'application/json',
                    'Content-Type' : 'application/json'
                },
                url: this.url,
                type: 'POST',
                data: JSON.stringify(envelope),
                error: (call, status, text) => {
                    this.updateSendingStatus({
                        sending: false,
                        error: text || status || 'Connection error'
                    });
                    // send failed, so re-populate the delta
                    for (let key of Object.keys(sendingDelta)) {
                        if (!this.delta.hasOwnProperty(key)) {
                            this.delta[key] = sendingDelta[key];
                        }
                    }
                    // try again
                    this.dataUpdated();
                },
                success: (record) => {
                    if (record.deleted) {
                        this.deleted = true;
                        document.location.href = '../../../../classification/classifications';
                        return;
                    } else if (record.withdrawn !== this.record.withdrawn) {
                        // changing of withraw status required complete reload
                        this.deleted = true; // not necessarily true but stops the reload warning
                        document.location.reload(true);
                        return;
                    }

                    let updateSet = new Set();

                    if (record.config) {
                        eKeys = eKeysBase.configCopy(record.config);
                        eKeys.forEach(eKey => {
                            if (eKey.config_updates) {
                                let original = jContent.find(`[entry="${eKey.key}"]`);
                                let replacement = this.createEntry(eKey);
                                replacement.css('display', original.css('display'));
                                original.replaceWith(replacement);
                                updateSet.add(eKey.key);
                            }
                        });
                        this.updateSummaryTable();
                    }
                    
                    this.publish_level = null;
                    record.evidence = this.record.evidence;
                    this.record = record;
                    this.messages = record.messages || [];
                    
                    Object.assign(this.delayedPatch, record.data);
                    let delayedPatch = {};
                    for (let key of Object.keys(this.delayedPatch)) {
                        if (typeof(this.delta[key]) !== 'undefined') {
                            // only override values that the user hasn't changed since
                            // presumable values the user changed will be updated
                            delayedPatch[key] = this.delayedPatch[key];
                        } else {
                            updateSet.add(key);
                            this.data[key] = this.delayedPatch[key];
                            if (key === this.currentContextKey) {
                                this.updateContext(key);
                            }
                        }
                    }
                    this.delayedPatch = delayedPatch;
                    if (updateSet.size) {
                        this.populateForm(updateSet);
                    }

                    this.updateSendingStatus({
                        sending: false, modifications: Object.keys(this.delta).length > 0,
                        error: null
                    });

                    this.updatePublishHistory();
                    this.updateLinks();
                    this.updateErrors();
                    this.updateCitations();
                    
                    try {
                        Flags.instance.refresh();
                    } catch (e) {
                        console.log(e);
                    }
                }
            });
            // only send updates a maximum of once every 2 seconds
        }, 2000),

        dataUpdated(force) {
            // localStorage.editingSubmission = JSON.stringify(this.data);
            // FIXME only send patch of data
            if (!force && (Object.keys(this.delta).length === 0 && !this.publish_level && !this.deleteRequest && !this.undeleteRequest)) {
                return;
            }
            this.updateSendingStatus({modifications: true});
            this.sendChanges();
        },

        updateTitle() {
            let appendLabelHeading = (label, valueElement) => {
                return $('<div>', {class: 'row no-gutters', html:[
                    $('<label>', {class:'col-3 text-right align-self-center', text: label}),
                    $('<div>', {class:'col-9 pl-3 align-self-center', html:valueElement})
                ]}).appendTo(jSyncStatus);
            };
            let appendLabelHeadingForKey = (key, showBlank, labelOverride, extra) => {
                let value = this.value(key);
                if (value !== null || showBlank) {
                    let eKey = eKeys.key(key);
                    let label = labelOverride || eKey.label;
                    let valueElement;
                    if (value === null) {
                        valueElement = $('<span>', {class: 'no-value', text: '-'});
                    } else {
                        value = eKey.prettyValue(value).val;
                        if (extra) {
                            value += extra;
                        }
                        valueElement = $('<span>', {text: value});
                    }
                    appendLabelHeading(label, valueElement);
                }
            };

            jSyncStatus.empty();

            if (this.lab_record_id) {
                appendLabelHeading('Lab ID', $('<span>', {text: this.lab_record_id}));
            }
            appendLabelHeadingForKey(SpecialEKeys.GENOME_BUILD, true, 'Build');

            if (this.record.allele && this.record.allele.resolved) {
                let resolved = this.record.allele.resolved;
                let hgvsDom = VCTable.format_hgvs(resolved);
                appendLabelHeading("Variant", hgvsDom);
                if (resolved.allele_info_id) {
                    appendLabelHeading('', $('<a>', {'data-toggle':'ajax-modal', href: Urls.view_imported_allele_info_detail(resolved.allele_info_id), text:'Resolution Details'}));
                }
            }

            let p_hgvs = this.value(SpecialEKeys.P_HGVS);
            if (p_hgvs) {
                let p_dot = p_hgvs.indexOf('p.');
                if (p_dot !== -1) {
                    p_hgvs = p_hgvs.substring(p_dot);
                }
                appendLabelHeading('p.HGVS', $('<span>', {text: p_hgvs}));
            }

            let alleleOriginBucket = this.record.allele_origin_bucket;
            let alleleOriginDisplay = {
                "U": "not-set",
                "G": "GERMLINE (default)",
                "S": "SOMATIC (default)"
            }[alleleOriginBucket];
            let alleleOriginValue = this.value(SpecialEKeys.ALLELE_ORIGIN);
            if (alleleOriginValue) {
                let eKey = eKeys.key(SpecialEKeys.ALLELE_ORIGIN);
                alleleOriginDisplay = eKey.prettyValue(alleleOriginValue).val;
                // TODO hopefully replace this with an attribute of the allele origin drop down
                // also see classification.py def calc_allele_origin_bucket(self) -> AlleleOriginBucket:
                if (alleleOriginValue.indexOf('germline') !== -1) {
                    alleleOriginBucket = "G";
                } else if (alleleOriginValue.indexOf('somatic') !== -1) {
                    alleleOriginBucket = "S";
                } else {
                    alleleOriginBucket = "U";
                }
            }

            let bucketHtml = $(`<div class="allele-origin-box horizontal allele-origin-${alleleOriginBucket}">
                <div class="allele-origin-text">
                    ${alleleOriginDisplay}
                </div>
            </div>`);

            $('<hr/>').appendTo(jSyncStatus);
            appendLabelHeading("Origin", bucketHtml);

            appendLabelHeadingForKey(SpecialEKeys.ASSERTION_METHOD, true, "Method")

            if (this.record.resolved_condition) {
                let label = "Condition";
                if (alleleOriginBucket === "S") {
                    label = "Cancer";
                }
                appendLabelHeading(label, VCForm.format_condition(this.record.resolved_condition));
            }

            appendLabelHeadingForKey(SpecialEKeys.CLINICAL_SIGNIFICANCE, true, 'Class.');

            if (this.value(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)) {
                let extra = null;
                for (let level of ['a', 'b', 'c', 'd']) {
                    let value = this.value(`amp:level_${level}`);
                    if (value) {
                        extra = ` ${level.toUpperCase()}`;
                        break;
                    }
                }
                appendLabelHeadingForKey(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE, true, 'Somatic Sig', extra);
            }

            if (this.record.sample_id) {
                let href = Urls.view_sample(this.record.sample_id);
                let sampleElement = $('<a>', {class:'hover-link', text: this.record.sample_name, href:href});
                sampleElement = $('<span>', {html: [
                    sampleElement,
                ]})
                appendLabelHeading('Sample', sampleElement);
            }

            if (this.userAdmin) {
                let adminLinkHref = Urls["admin:classification_classification_change"](this.record.id);
                let adminLink = $('<a>', {class:'admin-link', target: '_blank', html: 'Manage', href: adminLinkHref});
                appendLabelHeading('Admin', adminLink);
            }

            $('<div>', {id:'vc-quick-submit'}).appendTo(jSyncStatus);
        },

        variantCoordinate() {
            let vc_value = this.value(SpecialEKeys.VARIANT_COORDINATE);
            if (vc_value) {
                return vc_value;
            }
            return this.record.allele.variant_coordinate;
        },

        cHGVS() {
            vc_value = this.value(SpecialEKeys.C_HGVS);
            if (vc_value) {
                return vc_value;
            }
            if (this.allele && this.allele.resolved) {
                return this.allele.resolved.full;
            }
            return this.value(SpecialEKeys.G_HGVS);
        },

        populateForm(restrictToKeysSet) {
            let vcform = this;
            jContent.find('[name]').each(function(index, elem) {
                elem = $(elem);
                let name = elem.attr('name');
                if (restrictToKeysSet && !restrictToKeysSet.has(name)) {
                    return;
                }
                
                if (elem.closest('#fileupload').length > 0) {
                    return;
                }

                let visualElement = elem;
                if (elem.prop('tagName') === 'SELECT') {
                    let container = elem.siblings('.chosen-container');
                    if (container.length) {
                        visualElement = container;
                    }
                }

                let val = vcform.value( name );
                let immutable = vcform.immutable( name );
                
                let immutableDiv = visualElement.siblings('.immutable');
                if (immutable) {
                    if (!vcform.isEditMode()) {
                        visualElement.parent().addClass('entry-readonly');
                    }
                    visualElement.hide();
                    if (!immutableDiv.length) {
                        immutableDiv = $('<div>', {class: 'immutable'});
                        let isTextArea = elem.prop('tagName') === 'TEXTAREA';
                        if (isTextArea) {
                            immutableDiv.addClass('immutable-textarea');
                        }
                        immutableDiv.insertAfter( visualElement );
                    }
                    if (vcform.hidden(name)) {
                        immutableDiv.addClass('hidden-value');
                        val = 'hidden';
                    }
                    eKeys.key(name).formatValue(val, immutableDiv);
                    
                } else {
                    immutableDiv.hide();
                    visualElement.show();
                    visualElement.parent().removeClass('entry-readonly');
                }

                // populate inputs even if they're hidden, as they might have
                // widget functionality (such as scientific notation)
                if (elem.attr('customPopulate')) {
                    elem.trigger('onpopulate', val);
                } else {
                    if (elem[0].nodeName === 'SELECT') {
                        vcform.populateSelect(elem[0], val);
                    } else if (elem[0].nodeName === 'SPAN') {
                        elem.age('value', val);
                    } else {
                        if (elem.val().trim() === val) {
                            // don't remove whitespace that the user has entered
                            // e.g. "One day I " would be changed to "One day I" before the user had a chance
                            // to finish (if the auto-update kicked in).
                        } else {
                            elem.val(val);
                            if (elem.refresh) {
                                elem.refesh();
                            }
                        }
                        // elem.keyup(); //keyup to trigger scientific notation
                    }
                }

            });
            if (!restrictToKeysSet) {
                this.filter('');
                this.dataUpdated();
            }
        },

        populateSelect(select, val) {
            select = $(select);
            if (val != null) {
                let options = Array.apply(null, select[0].options);
                
                if (!Array.isArray(val)) {
                    val = [val];
                }
                
                // find all the values that have no matching options
                let customValues = val.filter(v => !options.some(o => o.value === `${v}`));

                if (customValues.length) {
                    let appendTo = select;
                    let optGroups = select.children('optgroup');

                    if (optGroups && optGroups.length > 0) {
                        appendTo = $(optGroups[ optGroups.length - 1]);
                    } else {
                        appendTo = $('<optgroup>', {label: 'Illegal Value'});
                        select.append(appendTo);
                    }
                    for (const custom of customValues) {
                        appendTo.prepend($('<option>', {value: custom, text: `${custom}`}));
                    }
                }
            }
            select.val( val );
            select.trigger("chosen:updated");
        },

        immutable(key) {
            if (!this.isEditMode() || this.record.withdrawn) {
                return true;
            }
            let blob = this.data[key] || {};
            return !!blob.immutable || !!blob.hidden;
        },
        
        hidden(key) {
            let blob = this.data[key] || {};
            return !!blob.hidden;
        },

        value(key, value) {
            if (value === undefined) {
                return emptyToNull( (this.data[key] || {}).value );
            }
            if (this.immutable(key)) {
                // double check that we're not setting values for immutable values
                return;
            }

            value = emptyToNull(value);
            let eKey = eKeys.key(key);
            if (eKey.value_type === 'B' && value != null) {
                value = value == 'true';
            }
            let blob = this.data[key] || {};
            if (blob.value === value) {
                return;
            }
            blob.value = value;
            
            this.data[key] = blob;
            this.delta[key] = blob;
            
            this.dataUpdated();
            this.updateContext(key);
        },
        
        explain(key) {
            return emptyToNull( (this.data[key] || {}).explain );
        },
        
        valueLabel(key) {
            let value = this.value(key);
            let ekey = eKeys.key(key);
            
            return ekey.prettyValue(value).val;
        },
        
        note(key, note) {
            if (note === undefined) {
                return emptyToNull( (this.data[key] || {}).note );
            }
            note = emptyToNull(note);
            let noteDiv = $(`#note-${key}`);

            this.fixNotePopover(noteDiv, note);
            
            let blob = this.data[key] || {};
            blob.note = note;
            
            this.data[key] = blob;
            this.delta[key] = blob;
            this.updateContext(key);
            
            this.dataUpdated();
        },

        fixNotePopover(noteDom, note) {
            noteDom.attr('data-content', note);
            if (note) {
                noteDom.attr('data-content', note);
                noteDom.attr('data-toggle', 'popover');
                noteDom.attr('title', 'Note');
                noteDom.popover({trigger:'hover'});
                noteDom.removeClass('far');
                noteDom.addClass('fas');
                noteDom.removeClass('d-none');
            } else {
                noteDom.attr('data-toggle', null);
                noteDom.attr('title', null);
                noteDom.popover('dispose');
                noteDom.removeClass('fas');
                noteDom.addClass('far');
                if (!this.isEditMode()) {
                    noteDom.addClass('d-none');
                }
            }
        },
        
        refs(key) {
            let blob = this.data[key] || {};
            return blob.db_refs || [];
        },

        init: function(params, _eKeys, record) {

            let vcform = this;

            jContent = $(params.content);
            jSyncStatus = $(params.summary);
            jErrors = $(params.errors);
            jPublishHistory = $(params.publishHistory);
            jLinks = $(params.links);
            jShareButtons = $(params.shareButtons);
            this.userAdmin = params.userAdmin;
            this.lab_record_id = params.lab_record_id;
            this.citations = params.citations;
            this.attachmentsEnabled = params.attachmentsEnabled || false;
            
            jHelp = $(params.help);
            jCritTable = $(params.critTable);
            jFilterBox = $(params.filterBox);
            jClearFilterButton = $(params.clearFilterButton);
            jClearFilterButton.hide();

            eKeysBase = _eKeys;  // so we can copy it again when config changes
            eKeys = _eKeys.configCopy(record.config);
            vcLinks = new VCLinks(eKeys);
            vcform.url = Urls.classification_api();
            vcform.genomeBuild = params.genomeBuild;

            this.otherClassificationsSummary = params.otherClassificationsSummary;
            this.reportEnabled = params.reportEnabled;
            this.deleteEnabled = params.deleteEnabled;
            this.withdrawReasons = params.withdrawReasons;

            this.initData(record);
            // Create the sections
            let firstFamily = true;
            for (let familyKey of Object.keys(EKey.families)) {
            
                if (familyKey === 'U' && !eKeys.hasUnknown) {
                    continue;
                }

                let firstCss = firstFamily ? " in" : "";
                let label = EKey.families[familyKey];
                jContent.append(
                    $('<div>', {class: 'card', family: familyKey, html:[
                        $('<div>', {class:'card-header collapsed' + firstCss, "data-toggle":"collapse", "data-target": `#section-${familyKey}`, html:[
                            $('<a>', {class:'card-title', text: label})
                        ], click: (e) => {this.toggleCheck(e)}}),
                        $('<div>', {id:`section-${familyKey}`, class: "panel-collapse collapse" , 'data-parent': `#${jContent.attr('id')}`, html: [
                            $('<div>', {class:`card-body sub-content sub-content-${familyKey}`, family: familyKey})
                        ]})
                    ]})
                );
                firstFamily = false;
            }
            // fill the data in the sections
            jContent.find(".sub-content[family]").each((index, elem) => {
                this.populateFamily(elem);
            });

            if (this.attachmentsEnabled) {
                let familyKey = 'uploads';
                jContent.append(
                    $('<div>', {family: familyKey, class: 'card', html:[
                        $('<div>', {class:'card-header collapsed', "data-toggle":"collapse", "data-target": `#section-${familyKey}`, html:[
                            $('<a>', {class:'card-title', text: 'Uploads'})
                        ], click: (e) => {this.toggleCheck(e)}}),
                        $('<div>', {id:`section-${familyKey}`, class: "panel-collapse collapse", 'data-parent': `#${jContent.attr('id')}`, html: [
                            $('<div>', {class:`card-body sub-content sub-content-${familyKey}`, family: familyKey})
                        ]})
                    ]})
                );
                $('#upload-section').detach().appendTo($(`#section-${familyKey}`));
            }

            jContent.find('[data-parent]').on('shown.bs.collapse', function() {
                if (!filtering) {
                    let scrollMe = $('.main-content');
                    let scrollTo = $(this).closest('.card');
                    let scrollToOffset = scrollTo.position().top;
                    // scrollMe.scrollTop(scrollToOffset + 15);
                    scrollMe.animate({scrollTop: scrollToOffset + 15}, 200);
                } else {
                    // console.log('STILL FILTERING, NO JUMPINT');
                }
            });

            jContent.find("h3").attr('tabindex', 0);
            let entries = jContent.find("[entry]");
            entries.focusin(function() {
                vcform.updateContext($(this).attr('entry'));
            });
            entries.click(function() {
                vcform.updateContext($(this).attr('entry'));
            });

            jContent.find("label").click(function() {
                let entry = $(this).closest('[entry]');
                entry.find('input').focus();
                entry.find('[tabindex]').focus();
            });

            jFilterBox.keyup(function() {
                let filterValue = $(this).val().trim();
                vcform.filter( filterValue );
                if (filterValue.length) {
                    jClearFilterButton.show();
                } else {
                    jClearFilterButton.hide();
                }
            });
            jFilterBox.focusin(function() {
                window.setTimeout(() => {
                    jFilterBox.select();
                }, 1);
                jHelp.closest('.card').find('.card-title').text('Filter');
                jHelp.empty()
                    .append($('<div>', {class: 'description', html: 'Type in the filter to search on field names<br/>Type * to see the entire populated form at once.<br/>Type ** to see all possible fields at once.'}));
            });

            this.populateForm();
            this.updateErrors();
            this.updateSummaryTable();
            this.setupUnloadWarning();

            $('#vc-extras').on('mousewheel', function(event) {
                let scrollMe = $('.main-content');
                let scroll = scrollMe.scrollTop();
                scrollMe.scrollTop(scroll - event.originalEvent.wheelDeltaY);
                return false;
            });
        },

        toggleCheck: function(e) {
            if (filtering) {
                e.stopPropagation();
            }
        },

        // The actual show/hide of elements on the form
        filter: function(val) {
            let oldFiltering = filtering;
            let showUploads = null;
            if (val.length === 0) {
                if (filtering) {
                    showUploads = true;
                }
                filtering = false;
            } else {
                showUploads = 'uploads'.startsWith(val.toLowerCase());
            }
            if (!this.attachmentsEnabled) {
                showUploads = false;
            }

            val = val.toLowerCase();
            //let inst = jContent.accordion('instance');
            eKeys.forEach(eKey => {
                let matches = false;
                if (eKey.exclude_namespace && this.value(eKey.key) == null) {
                    matches = false;
                } else {
                    if (val.length === 0 || val === '**') {
                        matches = true;
                    } else if (val === '*') {
                        matches = this.value(eKey.key) !== null || this.note(eKey.key) !== null || this.severity(eKey.key) !== null;
                    } else if (eKey.matchesFilter(val)) {
                        matches = true;
                    }
                }
                let elem = $(`[entry="${eKey.key}"]`);

                filtered[eKey.key] = !matches;
                let shouldHide = filtered[eKey.key];

                if (!shouldHide &&
                    !filtering &&
                    eKey.hide === true &&
                    this.value(eKey.key) === null &&
                    this.note(eKey.key) === null) {
                    shouldHide = true;
                }
                if (!shouldHide) {
                    elem.show();
                } else {
                    elem.hide();
                }
            });

            if (showUploads === true) {
                $(`[family=uploads]`).show();
            }

            if (val.length) {
                filtering = true;
            } else {
                filtering = false;
            }
            this.showHideFamilies();
            if (oldFiltering !== filtering) {
                let allPanels = $('#vc-form [data-parent].collapse');
                let scrollMe = $('.main-content');
                if (filtering) {
                    this.previously_expanded = allPanels.filter('.show');
                    $('#vc-form .accordion').removeClass('indicator-plus-before');
                    allPanels.addClass('show');
                    scrollMe.scrollTop(0);
                } else {
                    $('#vc-form .accordion').addClass('indicator-plus-before');
                    allPanels.not(this.previously_expanded).removeClass('show');
                    if (this.previously_expanded && this.previously_expanded.length) {
                        //this.previously_expanded.get(0).scrollIntoView();
                        window.setTimeout(() => {
                            scrollMe.scrollTop(this.previously_expanded.closest('.accordion').position().top + 15);
                        },1);
                    }
                }
            }
            if (showUploads === false) {
                $(`[family=uploads]`).hide();
            }
        },

        showHideFamilies: function() {
            let showing = {};
            eKeys.forEach(eKey => {
                if (eKey.exclude_namespace && this.value(eKey.key) == null) {
                    return;
                }
                let shouldHide = filtered[eKey.key];
                if (!shouldHide) {
                    showing[eKey.evidence_category] = true;
                }
            });
            Object.keys(EKey.families).forEach(family => {
                let members = jContent.find(`.card[family='${family}']`);
                if (showing[family]) {
                    members.removeClass('no-results');
                } else {
                    members.addClass('no-results');
                }
            });
        },

        refreshContext: function() {
            this.updateContext(this.currentContextKey);
        },
        /**
         * Updates the help context box based on an entry element
         */
        updateContext: function(key) {
            let thisForm = this;
            this.currentContextKey = key;
            let value = this.value(key);
            let note = this.note(key);
            let explain = this.explain(key);
            let refs = this.refs(key);
            let eKey = eKeys.key(key);
            let immutable = this.immutable(key);

            jHelp.closest('.card').find('.card-title').text(eKey.label);
            jHelp.empty();
            let content = jHelp;
            
            let valueElement = $('<span>');
            eKey.formatValue(value, valueElement);
            if (valueElement.text() === '') {
                valueElement = $('<span>', {class:'no-value', text:'blank'});
            }

            $('<div>', {class: 'my-1', style:'white-space:pre-wrap', html:[
                $('<label>', {class:'mr-2', text: 'Value: '}),
                valueElement
            ]}).appendTo(content);

            if (note !== null) {
                $('<div>', {class: 'note text-secondary my-1', html:[
                    $('<label>', {class:'mr-2', text: 'Note: '}),
                    $('<span>', {text:note})
                ]}).appendTo(content);
            }
            if (refs !== null && refs.length) {
                let refsDom = $('<div>', {class: 'refs'}).appendTo(content);
                for (let ref of refs) {
                    this.renderReference(ref).appendTo(refsDom);
                }
            }
            
            $('<div>', {class: 'value-end'}).appendTo(content);
    
            if (eKey.version) {
                content.append(
                    $('<div>', {class: "my-1"}).append($('<label>', {text: 'Version: '}), $('<span>').html(eKey.version))
                );
            }
            if (eKey.examples && !immutable) {
                eKey.examples.forEach(example => {
                    content.append(
                        $('<div>', {class: "my-1"}).append($('<label>', {class:'mr-2', text: 'Example: '}), $('<span>').html(example))
                    );
                });
            }
            
            if (explain) {
                content.append(
                    $('<div>', {class: "description mt-2"}).html([
                        $('<div>', {class: 'explain full'}),
                        $('<div>', {text: 'Lab specific description'}).css('font-style', 'italic').css('display', 'inline-block'),
                        $('<div>', {class: 'explain-content', text: explain})
                    ])
                );
            } else {
                descriptionSpan = eKey.description ? EKeys.fixDescription(eKey.description) : $('<i>', {text:'No help is provided for this field'});
                if (this.isEditMode()) {
                    if (eKey.hide === true) {
                        descriptionSpan.append($('<br/><br/><i>This field is not shown by default for your lab.</i>'));
                    }
                }
                let description = $('<div>', {class: "description mt-2"}).appendTo(content);
                descriptionSpan.appendTo(description);

                if (eKey.see) {
                    content.append([
                        $('<br/>'),
                        $('<a>', {class: 'hover-link external-link', href: eKey.see, text: eKey.see, target: eKey.key})
                    ]);
                }
            }
            jHelp.css({'max-height': 660});

            if (this.helpOverflowAmount() > 0) {
                $('<div>', {
                    style: 'position:absolute; bottom: 0px; right: 0px; padding: 4px; background-color: white; border-radius: 4px',
                    html: $('<a>', {text: 'Click for full content', class: 'hover-link'}).click(function(e) {

                        function titledValue(title, value) {
                            return $("<div>", {class: 'mb-4', html:[
                                $("<label>", {"text": title, "style": "font-weight:600"}),
                                $("<hr>", {"style": "margin-top:0.5rem;margin-bottom:0.5rem"}),
                                $("<div>", {class:'text-body', html: value})
                            ]});
                        }
                        let helpHtml = eKey.description ? EKeys.fixDescription(eKey.description) : $('<i>', {text:'No help is provided for this field'});
                        let popupContent = $('<div>');
                        popupContent.append(titledValue("Description", helpHtml));

                        if (eKey.see) {
                            popupContent.append(titledValue("Relevant Link", $('<a>', {class: 'hover-link external-link', href: eKey.see, text: eKey.see, target: eKey.key})));
                        }
                        if (explain) {
                            popupContent.append(titledValue("Lab Specific Details", explain));
                        }
                        if (note) {
                            popupContent.append(titledValue("Note", explain));
                        }
                        let valueHtml;
                        if (value) {
                            valueHtml = $("<div>", {style:'white-space:pre-wrap;word-break: break-word;', html: eKey.formatValue(value)});
                        } else {
                            valueHtml = $("<i>", {class:'no-value', text:"No Value"});
                        }
                        popupContent.append(titledValue("Value", valueHtml));

                        if (refs !== null && refs.length) {
                            let refsDom = $('<ul>', {class: 'refs'});
                            for (let ref of refs) {
                                $('<li>', {html: thisForm.renderReference(ref) }).appendTo(refsDom);
                            }
                            popupContent.append(titledValue("References", refsDom));
                        }

                        createModal('info', eKey.label, popupContent);
                    })
                }).appendTo(jHelp);
            }
            // if (this.helpOverflowAmount() > 0) {
            //    jHelp.css({'max-height': jHelp[0].scrollHeight});
            //}

        },

        helpOverflowAmount: function() {
            let box = jHelp[0];
            return box.scrollHeight - box.clientHeight;
        },
        
        promptNote: function(key) {
            if (!this.isEditMode() || this.record.withdrawn) {
                return;
            }
        
            let vcform = this;
            let eKey = eKeys.key(key);

            let dialogContent = createModalShell('note',eKey.label, 'lg');
            let body = dialogContent.find('.modal-body');
            $('<p>', {text: 'Edit the note for this field below'}).appendTo(body);
            let textArea = $('<textarea>', {text: this.note(key), class:'form-control', rows:5}).appendTo(body);
            let footer = dialogContent.find('.modal-footer');
            footer.html([
                $('<button>', {type:"button", class:"btn btn-secondary", 'data-dismiss':"modal", text:'Cancel'}),
                $('<button>', {type:"button", class:"btn btn-primary", text:'OK', click: () => {
                    vcform.note(key, textArea.val());
                    dialogContent.modal('hide');
                }}),
            ]);

            dialogContent.on('hidden.bs.modal', function() {
                dialogContent.modal('dispose');
                dialogContent.remove();
            });
            dialogContent.on('shown.bs.modal', function() {
                textArea.focus();
            });
            dialogContent.modal({focus:true});
        },

        createEntry: function(eKey) {
            let key = eKey.key;
            let vcform = this;
            
            let type = eKey.value_type;
            switch (type) {
                case 'B': type = 'bool'; break;
                case 'C': type = 'crit'; break;
                case 'D': type = 'date'; break;
                case 'S': type = 'select'; break;
                case 'M': type = 'multiple'; break;
                case 'A': type = 'age'; break;
                case 'T': type = 'textarea'; break;
                case 'U': type = 'user'; break;
                case 'N': type = 'unit'; break;
                // case N Unit (0 to 100)
                default: type = 'input';
            }


            let label = eKey.label;
            if (eKey.mandatory) {
                label = '*' + label;
            }
            let labelDom =  $('<div>', {class: 'text-right align-self-center', id: `label-${key}`, html:[
                $('<label>', {text: label})
            ]});
            if (type === 'textarea') {
                labelDom.addClass('col-12 mb-2');
            } else {
                labelDom.addClass('col-4');
            }

            if (eKey.sub_label) {
                labelDom.append($('<div>', {class: 'text-info', text: eKey.sub_label}));
            }

            let noteText = this.note(key);
            let explain = this.explain(key);
            let hasNote = !!noteText;
            let hasExplain = !!explain;

            let widgetDiv = $('<div>', {class: "col-7"});

            let noteDiv = $('<i/>', {
                class: 'text-muted fa-comment-alt note',
                id: `note-${key}`,
                click: () => {this.promptNote(key);}
            });
            this.fixNotePopover(noteDiv, noteText);

            let infoDiv = $('<i/>', {
                id: `explain-${key}`,
                title: 'This field has a custom description. Click on it and view the blue help box for details.',
                class: `fas fa-info-circle text-muted explain`,
            });
            if (!hasExplain) {
                infoDiv.addClass('d-none');
            }

            let entry = $('<div>', {class:`row entry entry-${type} form-group`, entry:key, html:[
                labelDom,
                widgetDiv,
                $('<div>', {class: 'col-1 text-nowrap p-0 d-flex', html: [noteDiv, infoDiv]})
            ]});
            entry.css('align-items', 'center');

            switch (type) {
            
                case 'bool':
                case 'crit':
                case 'select':
                case 'multiple': {
                    let value = this.value(key);
                    let optionSources = eKey.options || [];

                    // move all disabled options to the bottom
                    optionSources = _.flatten(_.partition(optionSources, o => !o.exclude_namespace));
                    if (type === 'bool') {
                        optionSources = [
                            {key: 'true', label: 'True'},
                            {key: 'false', label: 'False'}
                        ];
                    }

                    let optGroups = [];

                    let standardPrefix = ""
                    let overridePrefix = ""
                    let emptyValuePrefix = "";
                    if (type === "crit") {
                        overridePrefix = "❗";
                        standardPrefix = "⭐ ";
                        emptyValuePrefix = "⚪ ";
                    }


                    let blankOption = optionSources.find(o => !('key' in o));

                    let makeOption = (option, overridePrefix) => {
                        let prefix = "";
                        if (option.key === "NM" || option.key === "NA") {
                            prefix = emptyValuePrefix;
                        } else if (option.key) {
                            prefix = standardPrefix;
                        }
                        if (overridePrefix) {
                            prefix = overridePrefix;
                        }

                        let optDom = $('<option>', {
                            value: option.key || '',
                            text: prefix + (option.label || EKey.prettyKey(option.key))
                        });
                        optDom.prop("disabled", !!option.exclude_namespace);

                        return optDom;
                    }

                    if (optionSources.find(o => o.override)) {
                        let optGroupNormal = optionSources.filter(o => !o.override).map(option => {
                            return makeOption(option);
                        });
                        let optGroupOverride = optionSources.filter(o => o.override).map(option => {
                            return makeOption(option, overridePrefix);
                        });
                        optGroups.push($('<optgroup>', {label: 'standard values', html: optGroupNormal}));
                        optGroups.push($('<optgroup>', {label: 'override values', html: optGroupOverride}));
                    } else {
                        options = optionSources.map(option => {
                            return makeOption(option);
                        });
                        optGroups = optGroups.concat(options);
                    }

                    if (eKey.allow_custom_values) {
                        optGroups.push($('<optgroup>', {
                            label: 'custom values',
                            html: $('<option>', {value: '!CUSTOM!', text: 'Enter Custom Value'})
                        }));
                    }

                    if (!blankOption) {
                        options.unshift($('<option>', {value: '', text: ''}));
                    }

                    let select = $('<select>', {name: key, html: options.concat(optGroups)});
                    if (type === 'multiple') {
                        select.attr('multiple', true);
                    }
                    widgetDiv.append(select);

                    select.chosen({
                        allow_single_deselect: true,
                        width: '288px',
                        placeholder_text_single: ' '
                    }).change(function () {

                        let custom = false;
                        let value = emptyToNull($(this).val());

                        let valueArray = value;
                        if (!Array.isArray(valueArray)) {
                            valueArray = [valueArray];
                        }
                        if (valueArray.indexOf('!CUSTOM!') !== -1) {
                            valueArray = valueArray.filter(v => v !== '!CUSTOM!');
                            custom = true;
                            let customVal = emptyToNull(prompt('Please enter the custom value'));
                            valueArray.push(customVal);
                        }

                        if (type !== 'multiple') {
                            value = valueArray[0];
                        } else {
                            value = valueArray;
                        }
                        let key = $(this).attr('name');

                        vcform.value(key, value);

                        if (custom) {
                            this.blur();
                            vcform.populateSelect(this, value);
                        }
                        if (type === 'crit' || $(this).prop('name') === SpecialEKeys.ASSERTION_METHOD) {
                            vcform.updateSummaryTable();
                        }
                    });

                    break;
                }
                case 'textarea': {
                    labelDom.removeClass('text-right');
                    let input =
                        $('<textarea>', {
                            class: 'form-control',
                            name: key,
                            rows: 4,
                            placeholder: eKey.placeholder || '',
                            keyup: function () {
                                vcform.value(key, $(this).val());
                            },
                            change: function (val) {
                                vcform.value(key, $(this).val());
                            },
                        });
                    widgetDiv.append(input);
                    widgetDiv.removeClass('col-6').addClass('col-11');
                    break;
                }
                case 'age':
                {
                    let ageSpan = $('<span>', {name: key});
                    widgetDiv.append(ageSpan);
                    ageSpan.age({
                        updated: (event, data) => {
                            vcform.value(key, data.value);
                        }
                    });
                }
                break;
                default: {
                    let input =
                        $('<input>', {
                            class: 'form-control',
                            name: key,
                            placeholder: eKey.placeholder || '',
                            keyup: function () {
                                vcform.value(key, $(this).val());
                            },
                            change: function (val) {
                                vcform.value(key, $(this).val());
                            }
                        });
                    if (type === 'date') {
                        input.datepicker({dateFormat: "yy-mm-dd"});
                    }

                    widgetDiv.append(input);

                    if (type === 'unit') {
                        //$('<div>', {text: '(0-1)', class: 'unit', title: 'This value should be between 0 and 1'}).appendTo(entry);
                        input.scientific({
                            placeholder_short: '(0-1)',
                            placeholder: 'Enter a frequency between 0 and 1. You will be shown a % representation to the right.',
                            tooltip: 'This is the % representation of the inputted frequency'
                        });
                    }
                }
            }

            return entry;
        },

        populateFamily: function(divy) {
            divy = $(divy);
            let family = divy.attr('family');
            eKeys.forEach(eKey => {
                if (eKey.evidence_category === family) {
                    const entry = this.createEntry(eKey);
                    divy.append(entry);
                }
            });
        },

        clearFilter: function() {
            jFilterBox.val('');
            jFilterBox.keyup();
        },
        
        searchAsterisk: function() {
            window.setTimeout(() => {
                jFilterBox.val('*');
                jFilterBox.keyup();
            },0);
        },

        setAccordionIndex: function(family) {
            let segment = $(`.card[family=${family}] .collapse`);
            if (filtering) {
                this.clearFilter();
                window.setTimeout(() => {
                    segment.collapse('show');
                }, 1);
            } else {
                segment.collapse('show');
            }
        },

        hasStrength: function(val) {
            return val !== null && val !== 'NM' && val !== 'NA' && val !== 'N';
        },

        evidenceWeights: function(strs) {
            let weights = [];
            let critValuesExtra = Object.assign({"BM": "Benign Moderate"}, EKey.critValues);
            for (let strength of Object.keys(critValuesExtra)) {
                if (!this.hasStrength(strength)) {
                    continue;
                }
                let count = strs[strength] || 0;
                if (count > 0) {
                    weights.push(`${count}x${strength}`);
                }
            }
            return weights.join(', ');
        },

        calculateOverall: function(strs) {
            const getStr = function(key) {
                return strs[key] || 0;
            };
            const PVS = getStr('PVS');
            const PS = getStr('PS');
            const PM = getStr('PM');
            const PP = getStr('PP');
            const BA = getStr('BA');
            const BS = getStr('BS');
            const BP = getStr('BP');
            const BM = getStr('BM'); // note this aren't standard, so can't calculate anything with it
            const BX = getStr('BX');
            const PX = getStr('PX');

            let result = {
                P: null,
                LP: null,
                B: null,
                LB: null,
                US: null,
                X: null
            };

            // can't do ACMG calculation with Benign Moderate, non standard strength
            if (BM || BX || PX) {
                result.X = 1;
                return result;
            }

            if ((PVS || PS || PM || PP) && (BA || BS || BP)) {
                result.US = 2;
            } else {
                // Check for pathogenic
                if (
                    PVS >= 1 && (
                        PS >= 1 ||
                        PM >= 2 ||
                        (PM === 1 && PP === 1) ||
                        PP >= 2)) {
                    result.P = 1;
                } else if (
                    PS >= 2) {
                    result.P = 2;
                
                } else if (
                    PS === 1 && (
                        PM >= 3 ||
                        (PM === 2 && PP >= 2) ||
                        (PM === 1 && PP >= 4))
                    ) {
                    result.P = 3;
                }
                // Check likely pathogenic
                else if (PVS >= 1 && PM === 1) {
                    result.LP = 1;
    
                } else if (PS === 1 && (PM === 1 || PM === 2)) {
                    result.LP = 2;
    
                } else if (PS === 1 && PP >= 2) {
                    result.LP = 3;
    
                } else if (PM >= 3) {
                    result.LP = 4;
    
                } else if (PM === 2 && PP >= 2) {
                    result.LP = 5;
    
                } else if (PM === 1 && PP >= 4) {
                    result.LP = 6;
                }
    
                // check benign
                if (BA >= 1) {
                    result.B = 1;
                } else if (BS >= 2) {
                    result.B = 2;
                }
    
                // check likely benign
                else if (BS === 1 && BP === 1) {
                    result.LB = 1;
                } else if (BP >= 2) {
                    result.LB = 2;
                }
                
                if (Object.keys(result).filter(key => !!result[key]).length === 0) {
                    result.US = 1;
                }
            }

            return result;
        },

        classificationLabels: {
            P: 'Pathogenic',
            LP: 'Likely pathogenic',
            B: 'Benign',
            LB: 'Likely benign',
            US: 'Uncertain significance',
            X: 'Custom assertion method used (includes non-standard strength)'
        },
        roman: ['-','i','ii','iii','iv','v','vi'],

        formatOverall: function(result) {
            let footnote = (classification, note) => {
                let assertionMethod = this.value(SpecialEKeys.ASSERTION_METHOD);
                // FIXME maybe just check for "ACMG" in assertionMethod
                let row = $('<span>');
                if (assertionMethod && assertionMethod.indexOf('VCGS') !== -1) {
                    row.text(`Custom assertion method used.`);
                } else {
                    row.text('Calculation: ' + this.classificationLabels[classification]);
                    let noteLabel = '';
                    if (classification !== 'X') {
                        noteLabel = '(' + this.roman[note] + ')';
                    }
                    let calculation = $('<span>', {text: ' ' + noteLabel, title: `See ${row.text()} ${noteLabel} in ACMG Standards and Guidelines Table 5. Rules for combining criteria to classify sequence variants`}).appendTo(row);
                }
                return row;
            };

            let dom = $('<div>', {class:'calculated-result m-2'});
            for (let key of Object.keys(result)) {
                let value = result[key];
                if (value) {
                    dom.append(footnote(key, value));
                }
            }
            return dom;
        },

        updateSummaryTable: function() {
            let cellCount = {};
            let strengthCount = {};
            let score = 0;

            let familiesWithCrits = {};
            for (let crit of eKeys.criteria()) {
                familiesWithCrits[crit.evidence_category] = true;
                let critId = crit.key;
                let val = this.value(critId);
                let str = crit.defaultValue();

                let cell = crit.evidence_category + '-';
                if (this.hasStrength(val)) {
                    cell += val;
                    strengthCount[val] = (strengthCount[val] || 0) + 1;
                    let critPoints = EKey.strengthToPoints[val] || 0;
                    score += critPoints;
                } else {
                    cell += str;
                }
                let cellEntry = cellCount[cell] || {possible:[], neutral:[], notMet:[], actual:[]};
                cellCount[cell] = cellEntry;
                let label = crit.label;

                if (val == null && crit.hide) {
                    continue;
                }
                if (val == null) {
                    cellEntry.possible.push(label);
                } else if (val === 'N') {
                    cellEntry.neutral.push(label);
                } else if (!this.hasStrength(val)) {
                    cellEntry.notMet.push(label);
                } else {
                    cellEntry.actual.push(label);
                }
            }

            let assertionMethod = this.value(SpecialEKeys.ASSERTION_METHOD);
            let pointBased = false;
            // maybe look at the criteria to see if they're point based?
            if (assertionMethod && assertionMethod.indexOf('horak') !== -1) {
                pointBased = true;
            }

            let critTable = $('<table>');
            let headerRow = $('<tr>');
            headerRow.append($('<th></th>'));

            for (let key of Object.keys(EKey.critValues)) {
                if (!this.hasStrength(key)) {
                    continue;
                }
                let label = key;
                let th;
                if (pointBased) {
                    let points = EKey.strengthToPoints[key];
                    let title = "";
                    if (points != 0) {
                        let plural = Math.abs(points) > 1 ? "s" : "";
                        let direction = points < 0 ? "Benign" : "Oncogenic"
                        title = `(${Math.abs(points)} point${plural}) Towards ${direction}`;
                    } else {
                        title = "Neutral";
                    }
                    th = $(`<th>`, {text: `${(EKey.strengthToPoints[key])}`, title: title, class: 'col-header'});
                } else {
                    th = $(`<th>`, {text: key, class: 'col-header', title: EKey.critValues[key], 'data-toggle':"tooltip", 'data-placement': 'left'});
                }
                
                headerRow.append(th);
            }
            critTable.append(headerRow);

            let accordionWid = $(".content");
            Object.keys(EKey.families).forEach(fam => {
                if (!familiesWithCrits[fam]) {
                    return;
                }
                let familyLabel = EKey.families[fam];
                    
                let clickMethod = this.setAccordionIndex.bind(this, fam);
                let row = $('<tr>', {click: clickMethod});

                row.append($('<th>', {text: fam, class: 'row-header', title: familyLabel}));
                for (let str of Object.keys(EKey.critValues)) {
                    if (!this.hasStrength(str)) {
                        continue;
                    }
                    let key = fam + '-' + str;
                    let values = cellCount[fam + '-' + str] || {possible:[], notMet:[], neutral: [], actual:[]};
                    
                    let tooltip = [];
                    if (values.actual.length) {
                        tooltip.push( 'Met : ' + values.actual.join(', ') );
                    }
                    if (values.neutral.length) {
                        tooltip.push( 'Neutral : ' + values.neutral.join(', '));
                    }
                    if (values.notMet.length) {
                        tooltip.push( 'Not Met/Applicable : ' + values.notMet.join(', ') );
                    }
                    if (values.possible.length) {
                        tooltip.push( 'Not Evaluated : ' + values.possible.join(', ') );
                    }

                    let td = $('<td>', {title: tooltip.join('\n'), id:`table-${fam}-${str}`});
                    if (values.actual.length > 0) {
                        let direction = str[0];
                        td.addClass(`${direction}-cell`)
                        td.append($('<span>', {class: 'met', text: values.actual.length}));
                    }
                    if (values.notMet.length || values.possible.length) {
                        td.append($('<span>', {class: 'possible',
                            text: '/'.repeat(values.notMet.length) + '/'.repeat(values.neutral.length) + '?'.repeat(values.possible.length)
                        }));
                    }
                    if (values.actual.length === 0 && values.possible.length === 0) {
                        if (values.notMet.length === 0) {
                            td.addClass(`blank-cell`);
                        } else {
                            td.addClass('rejected-cell');
                        }
                    }
                    row.append(td);
                }
                critTable.append(row);
            });

            jCritTable.empty().append(critTable)

            if (!pointBased) {
                let result = this.calculateOverall(strengthCount);
                let evidenceWeightsDom = $('<div>', {class: 'm-2', text: this.evidenceWeights(strengthCount)});
                let resultDom = this.formatOverall(result);
                jCritTable.append(evidenceWeightsDom).append(resultDom);
            } else {
                let overallValue;
                if (score <= -7) {
                    overallValue = "Benign";
                } else if (score <= -1) {
                    overallValue = "Likely Benign";
                } else if (score <= 5) {
                    overallValue = "VUS";
                } else if (score <= 9) {
                    overallValue = "Likely Oncogenic";
                } else {
                    overallValue = "Oncogenic";
                }

                jCritTable.append(
                    $('<div>', {class: 'text-center my-2', html:`Calculated Score: ${score} ${overallValue}`})
                )
            }
        }
    };

    return VCForm;
})();

VCForm.format_condition = function(condition_json) {
    let dom = $('<span>');
    if (!condition_json) {
        return dom;
    }
    if (!condition_json.resolved_terms) {
        dom.append(condition_json.display_text);
        return dom;
    }

    let first = true;
    for (let term of condition_json.resolved_terms) {
        if (!first) {
            $('<br>').appendTo(dom);
        }
        first = false;
        $('<span>', {html: [
            $('<a>', {
                text: term.term_id,
                href:Urls.ontology_term(term.term_id.replace(':','_')),
                class: 'hover-link'
            }),
            " ",
            term.name
        ]}).appendTo(dom);
    }
    if (condition_json.resolved_terms.length > 1 && condition_json.resolved_join) {
        $('<span>', {class: 'font-italic', text:condition_json.resolved_join === 'C' ? ' Co-occurring' : ' Uncertain'}).appendTo(dom);
    }
    return dom;
};

let VCTable = (function() {
    let VCTable = function() {};
    VCTable.prototype = {};
    return VCTable;
})();

VCTable.format_hgvs = (parts) => {
    if (typeof(parts) == 'string') {
        parts = parts.trim();
        if (!parts.length) {
            return $('<span>', {text:'-', class:'no-value'});
        }
        return $('<span>', {text:limitLength(parts, 100)});
    } else if (!parts) {
        return $('<span>', {text:'-', class:'no-value'});
    }
    let genomeBuild = parts.genome_build;
    let transcript = parts.transcript;
    let geneSymbol = parts.gene_symbol;
    let cNomen = parts.c_nomen || parts.variant; // older code called cNomen 'variant'
    let allele = parts.allele;
    let variantId = parts.variant_id;
    let alleleId = parts.allele_id;
    let alleleInfoId = parts.allele_info_id;
    // let validationInclude = parts.validation_include;
    // let alleleInfoStatus = parts.allele_info_status;
    let icon = parts.icon;
    let tooltip = parts.tooltip;
    let pHgvs = parts.p_hgvs;
    let url = null;
    let error = parts.error;

    if (error) {
        return $(`<span><i class="fa-solid fa-circle-exclamation text-danger"></i> ${error}</span>`);
    }

    if (variantId) {
        url = Urls.view_allele_from_variant(variantId)
    } else if (alleleId) {
        url = Urls.view_allele(alleleId)
    }
    // also turn into a link

    let outterDom = $('<div>');
    let dom = $('<div>').appendTo(outterDom);

    if (allele) {
        outterDom.prepend($('<dom>', {class: 'font-weight-bold', text: allele}));
    }

    if (genomeBuild && (parts.desired === false || parts.normalized === false || parts.always_show_genome_build)) {
        let genomeBuildWrapper = $('<div>');
        let addNewLine = false;
        if (parts.normalized === false) {
            $('<span>', {html: 'not resolved, showing imported ', style:'color:#888'}).appendTo(genomeBuildWrapper);
        } else if (parts.desired === false) {
            $('<span>', {html: 'not lifted-over ', style:'color:#888'}).appendTo(genomeBuildWrapper);
        } else {
            addNewLine = true;
        }

        dom.append(genomeBuildWrapper);
        dom.append($('<span>', {text: genomeBuild, style:'font-weight:500; color:#888'}));
        if (addNewLine) {
            dom.append("<br/>");
        } else {
            dom.append(' ');
        }
        // <span style="white-space: nowrap"><span>{{ c_hgvs.transcript }}</span>{% if c_hgvs.gene_symbol %}(<span class="text-secondary" style="letter-spacing: 0.5px">{{ c_hgvs.gene_symbol }}</span>){% endif %}:</span><span style="display:inline-block;word-break: break-all">{{ c_hgvs.raw_c }}</span>
    }

    let cDom = $('<span>', {class:'c-hgvs-body'});
    if (transcript && cNomen) {
        cDom.append($('<span>', {text: transcript, class:'c-hgvs-transcript'}));
        if (geneSymbol) {
            cDom.append($('<span>', {class: 'c-hgvs-gene-symbol-b', html:[
                $('<span>', {text: "(", class: 'bracket'}),
                $('<span>', {class: 'c-hgvs-gene-symbol', text: geneSymbol}),
                $('<span>', {text: ")", class: 'bracket bracket-close'})
            ]}));
        }
        cDom.append($('<span>', {text: ":", class: 'colon'}));
        // used to be display:inline-block; but that doesn't underline
        cDom.append($('<span>', {class:'c-hgvs-nomen', text: limitLength(cNomen, 100)}));
    } else {
        cDom.append(limitLength(parts.full, 100));
    }
    if (url) {
        cDom = $('<a>', {href: url, html: cDom, class: 'hover-link'});
    }
    dom.append(cDom);
    if (icon) {
        let iconData = {"class": icon};
        if (tooltip) {
            iconData.title = tooltip;
        }
        let iconDom = $('<i>', iconData);
        iconDom.appendTo(dom);
    }

    if (pHgvs) {
        $('<span>', {class: 'd-block mt-1 text-secondary', text: limitLength(pHgvs)}).appendTo(dom);
    }
    return outterDom;
};

VCTable.hgvs = (data, type, row) => {
    return VCTable.format_hgvs(data).prop('outerHTML');
}

VCTable.condition = (data, type, row) => {
    if (data.resolved_terms) {
        return VCForm.format_condition(data).prop('outerHTML');
    } else {
        return data.display_text;
    }
};

VCTable.somatic_clinical_significance = (data, type, row) => {
    if (SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE in data) {
        let scs = data[SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE];
        let dom;
        if (scs) {
            let scsKey = EKeys.cachedKeys.key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE);
            let scsLabel = scsKey.prettyValue(scs);
            dom = $('<span>', {text: scsLabel.val, 'class': `c-pill scs scs-${scs}`});

            let highest_level = data["highest_level"];
            if (highest_level) {
                dom.append(` <span class="amp-level">${highest_level}</span>`);
            }

        } else {
            dom = $('<span>', {class: 'c-pill scs-none no-value', text: 'No Data'});
        }
        return dom.prop('outerHTML');
    } else {
        return "";
    }
}

VCTable.classification = (data, type, row) => {
    let dom;
    let cs = data[SpecialEKeys.CLINICAL_SIGNIFICANCE];
    let csKey = EKeys.cachedKeys.key(SpecialEKeys.CLINICAL_SIGNIFICANCE);
    let label = csKey.prettyValue(cs);
    let csClass =  `cs-` + (cs || '').toLowerCase()
    if (cs && cs.length) {
        dom = $('<span>', {class: `c-pill cs ${csClass}`, text: label.val});
    } else {
        dom = $('<span>', {class: 'c-pill cs-none no-value', text: 'No Data'});
    }
    return dom.prop('outerHTML');
};

VCTable.evidence_key = (key_name, data, type, row) => {
    let csKey = EKeys.cachedKeys.key(key_name);
    let span;
    if (data == null) {
        span = $('<span>', {class: 'no-value', text: '-'});
    } else {
        span = $('<span>');
        csKey.formatValue(data, span);
    }
    return span.prop('outerHTML');
};

VCTable.allele_origin_bucket_label = (allele_origin_bucket, override_text = "", alignment="vertical") => {
    let allele_origin_label = override_text;
    if (!allele_origin_label) {
        if (allele_origin_bucket == "S") {
            allele_origin_label = "SOMATIC"
        } else if (allele_origin_bucket == "G") {
            allele_origin_label = "GERMLINE"
        } else if (allele_origin_bucket == "U") {
            allele_origin_label = "UNKNOWN"
        } else {
            allele_origin_label = "???"
        }
    }
    return $('<div>', {
        class: `allele-origin-box ${alignment} allele-origin-${allele_origin_bucket}`,
        html: [
            $('<div>', {
                class: 'allele-origin-text',
                text: allele_origin_label
            })
        ]
    });
}

VCTable.identifier = (data, type, row) => {
    let id = data.id;
    let org_name = data.org_name;
    let lab_name = data.lab_name;
    let lab_record_id = data.lab_record_id;
    let shareLevel = data.share_level;
    let allele_origin_bucket = data.allele_origin_bucket;

    let alleleOriginDiv = VCTable.allele_origin_bucket_label(allele_origin_bucket);

    let shareInfo = EKeys.shareLevelInfo(shareLevel);
    let icon = $('<img>', {src: shareInfo.icon, class:'share-icon'});

    let link = $('<a>', {
        href: Urls.view_classification(id),
        class: 'hover-link',
        html: [
            icon,
            $('<span>', {text: `${org_name} / ${lab_name}`}),
            '<br/>',
            limitLengthSpan(lab_record_id, 30)
        ]
    });

    let dom;
    if (data.matches) {
        let searchTerm = data.search;
        let match_doms = [link];
        let ekeys = EKeys.cachedKeys;
        for (let [key, value] of Object.entries(data.matches)) {
            match_doms.push($('<div>', {class: 'search-result', html:[
                $('<span>', { text: ekeys.key(key).label + ' : '}),
                    highlightTextAsDom(searchTerm, value)
                ]
            }));
        }
        dom = $('<div>', {html: match_doms});
    } else {
        dom = link;
    }
    let indicatorClassName = `allele-origin-indicator allele-origin-horizontal allele-origin-${allele_origin_bucket}`;
    let fullDom = $('<div>', {style: 'display:flex; flex-direction:row', html:[
        alleleOriginDiv,
        dom
    ]});
    return fullDom.prop('outerHTML');
};

VCTable.sample = (sample_name, type, row) => {
    let dom = $('<span/>');
    if (sample_name && row.sample_id) {
        let link = $('<a>', {
            href: Urls.view_sample(row.sample_id),
            class: 'hover-link',
            html: [
                $('<span>', {text: sample_name}),
            ]
        });
        dom = link;
    }
    return dom.prop('outerHTML');
};