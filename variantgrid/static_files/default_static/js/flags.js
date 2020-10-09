function isScrolledToBottom(scrollDiv) {
    if (scrollDiv) {
        scrollDiv = $(scrollDiv).get(0);
        let scrolledUpBy = scrollDiv.scrollHeight - scrollDiv.scrollTop - scrollDiv.clientHeight;
        if (scrolledUpBy < 1) {
            return true;
        }
    }
    return false;
}

let Flags = (function () {

    let DataCollection = (function () {
        let DataCollection = function (initer) {
            this.initer = initer;
            this.dataMap = {};
        };
        DataCollection.prototype = {
            upsert(data) {
                let objy = this.dataMap[data.id];
                if (!objy) {
                    objy = this.initer();
                    this.dataMap[data.id] = objy;
                }
                Object.assign(objy, data);
                if (objy.updated) {
                    objy.updated();
                }

                return objy;
            },
            _remove(id) {
                delete this.dataMap[id];
            },
            get(id) {
                let existing = this.dataMap[id];
                if (!existing) {
                    throw new Error(`No record for id ${id}`);
                }
                return existing;
            },
            all() {
                return Object.values(this.dataMap);
            }
        };
        return DataCollection;
    })();

    let User = (function () {
        let User = function (db) {
            this.db = db;
        };
        User.prototype = {
            dom() {
                return $('<div>', {
                    class: 'user', html: [
                        $('<div>', { class: 'avatar' }).css(`background-image`, `url(${this.avatar})`).css('background-color', this.color),
                        $('<div>', { class: 'username', text: this.name })
                    ]
                });
            }
        };
        return User;
    })();
    
    let FlagResolution = (function() {
        let FlagResolution = function(db) {
            this.db = db;
        };
        FlagResolution.prototype = {};
        return FlagResolution;
    })();

    let FlagType = (function () {
    
        let FlagType = function (db) {
            this.db = db;
        };
        FlagType.prototype = {
            collectionObj() { return this.db.collections.get(this.collection) },

            dom(params) {
                params = params || {};
                let {flag, statusFlag} = params;

                let content = [];
                if (statusFlag) {
                    content.push($('<div>', { class: `flag ${statusFlag}` }));
                }
                let titleText = this.label;
                if (flag && !flag.creating) {
                    titleText += ` (${flag.ageText()})`;
                }
                let bonusClass = '';
                if (flag) {
                    bonusClass = `res-${flag.resolution}`
                }
                content = content.concat([
                    $('<div>', { class: `flag flag-${this.id} ${bonusClass}`}),
                    $('<span>', { class: 'title', text: titleText})
                ]);
                return $('<div>', {html: content});
            }
        };

        return FlagType;
    })();
    FlagType.EVENT = 'E';
    FlagType.TOGGLE = 'T';
    FlagType.REJECTABLE = 'R';
    FlagType.STANDARD = 'S';

    let FlagCollection = (function () {
        let FlagCollection = function (db) {
            this.db = db;
        };
        FlagCollection.prototype = {
            comments() {
                return this.db.comments.all().filter(c => c.flagObj().collection == this.id).sort((c1, c2) => c1.created - c2.created);
            },
            commentsImportant() {
                return this.db.comments.all().filter(c => {
                    let flag = c.flagObj();
                    if (flag.collection != this.id) {
                        return false;
                    }
                    return flag.flagTypeObj().importance >= 1;
                }).sort((c1, c2) => c1.created - c2.created);
            },
            flagTypes() {
                return this.db.flagTypes.all().filter(f1 => f1.context == this.context).sort((f1, f2) => {
                    return f1.label.localeCompare(f2.label);
                });
            },
            flags() {
                return this.db.flags.all().filter(f => f.collection == this.id && f.open)
                    .sort((f1, f2) => {
                        let localeCompare = -f2.flagTypeObj().label.localeCompare(f1.flagTypeObj().label);
                        return localeCompare ? localeCompare : f2.created - f1.created;
                    });
            },
            importantResolvedFlags() {
                return this.db.flags.all().filter(f => f.collection == this.id && !f.open && f.flagTypeObj().importance >= 2)
                    .sort((f1, f2) => {
                        let localeCompare = -f2.flagTypeObj().label.localeCompare(f1.flagTypeObj().label);
                        return localeCompare ? localeCompare : f2.created - f1.created;
                    });
            },
            raisableFlags() {
                flagTypes = this.flagTypes().filter(ft => {
                    let userPermission = this.user_permission;
                    let requiredPermission = ft.raise_permission;
                    return userPermission >= requiredPermission;
                });
                return flagTypes.map(ft => this.createFlagFor(ft)).filter(f => f != null);
            },
            createFlagFor(flagType) {
                if (flagType.only_one) {
                    let existing = this.db.flags.all().find(f => f.collection == this.id && f.flag_type == flagType.id);
                    if (existing) {
                        if (existing.status != 'O') {
                            return existing;
                        }
                        return null;
                    }
                }
                let newFlag = new Flag(this.db);
                Object.assign(newFlag, {
                    flag_type: flagType.id,
                    collection: this.id,
                    user: this.db.currentUser().id,
                    open: true,
                    creating: true
                });
                return newFlag;
            },
            remove() {
                for (let flag of this.flags()) {
                    flag.remove(flag.id);
                }
                this.db.collections._remove(this.id)
            },

            recent_activity() {
                // FIXME, filter to actually be recent
                return this.comments();
            },

            async watch_toggle() {
                let watching = this.watching === 0 || this.watching;
                await this.db.post(
                    {
                        watch: !watching
                    },
                    this.db.flagsUrl + this.id
                );
            }
        };
        return FlagCollection;
    })();

    let FlagTimeline = (function() {
        let FlagTimeline = function(flagCollection, dialog) {
            this.flagCollection = flagCollection;
            this.dialog = dialog;
            this.lastEntry = null;
            this.firstUpdate = true;
            this.scrollablePanel = null;
            this.title = "Timeline";
            this.body = null;
            this.footer = null;
        };
        FlagTimeline.prototype = {
            init() {
                let dom = $('<div>');
                this.body = dom;

                let content = $('<div>', { class: 'edit-flag-dialog container' }).appendTo(dom);
                let details = $('<div>', {class: 'flag-details'}).appendTo(content);
                this.scrollablePanel = details;
                
                $('<div>', { class: 'flag-comments' }).appendTo(details);
            },
            update() {
                // fixme, don't re-render comments already rendered
                let shouldScroll = (isScrolledToBottom(this.scrollablePanel) && !this.firstUpdate) || this.justPosted;

                let lastEntry = this.lastEntry;
                let reachedStart = lastEntry == null;
                let commentsDom = this.body.find('.flag-comments');
                let importantStuff = this.flagCollection.commentsImportant();

                if (!importantStuff.length) {
                    commentsDom.empty();
                    $('<p>', {text: 'There have been no important flag actions against this record.'}).appendTo(commentsDom);
                } else {
                    for (let comment of importantStuff) {
                        if (!reachedStart) {
                            reachedStart = comment.id == lastEntry.id;
                            continue;
                        }
                        comment.dom({ lastEntry: lastEntry, includeTitle: true, dialog: this.dialog }).appendTo(commentsDom);
                        lastEntry = comment;
                    }
                    this.lastEntry = lastEntry;
                    if (shouldScroll) {
                        this.scrollToBottom();
                    }
                }
                this.firstUpdate = false;
            },
            scrollToBottom() {
                window.setTimeout(() => {
                    $(this.scrollablePanel).get(0).scrollTop = 10000000;
                },1);
            }
        };
        return FlagTimeline;
    })();

    let FlagContent = (function() {
        let FlagContent = function(flag, dialog) {
            this.flag = flag;
            this.dialog = dialog;
            
            this.scrollablePanel = null;
            this.textarea = null;
            this.states = null;

            this.title = null;
            this.body = null;
            this.footer = null;

            this.cachedState = null;
            this.newActionDom = null;
            this.firstUpdate = true;
            this.justPosted = false;
        };
        FlagContent.prototype = {
            init() {
                this.title = this.flag.flagTypeObj().label;
                let dom = $('<div>');
                this.body = dom;

                let content = $('<div>', { class: 'edit-flag-dialog' }).appendTo(dom);
                let details = $('<div>', {class: 'flag-details'}).appendTo(content);
                this.scrollablePanel = details;
                
                $('<div>', { class: 'description border rounded p-4' }).appendTo(details);
                $('<div>', { class: 'flag-comments' }).appendTo(details);

                this.states = $('<div>');
                this.textarea = $('<textarea>', { class: 'form-control', placeholder: 'Enter a comment here', rows:5 }).hide();

                this.footer = $('<div>');

                this.newActionDom = $('<div>', {
                    class: 'new-action', html: [
                        this.states,
                        this.textarea
                    ]
                }).appendTo(content);

                this.footer = $('<button>', { class:'btn btn-primary w-100',  text: 'Save & Continue Editing', click: () => this.postAction() });
            },
            
            async postAction(params) {
                let flag = this.flag;
            
                params = params || {};
                let { close } = params;
                let comment = null;
                if (this.textarea) {
                    comment = this.textarea.val().trim();
                    if (comment.length == 0) {
                        comment = null;
                    }
                    this.textarea.val('');
                }
                let resolution = $(`select[name='flag-res']`).val() || $(`input[name='flag-res']:checked`).val();
                
                if (!comment && flag.creating) {
                    alert('Please enter a comment before raising this flag.');
                    // nothing to do, user probably has to add a comment
                    return;
                }
                this.setPosting();

                sendParams = {
                    comment: comment || null,
                };
                if (resolution) {
                    sendParams['resolution'] = resolution;
                }

                let url = null;
                if (flag.creating) {
                    sendParams.flag_type = flag.flag_type;
                    url = flag.db.flagsUrl + flag.collectionObj().id;
                } else {
                    url = flag.db.flagUrl + flag.id;
                }

                if (resolution || flag.creating) {
                    flag.db.addedOrClosedFlag = true;
                }

                let response = await flag.db.post(
                    sendParams,
                    url
                );
                await flag.db.refresh();
                
                this.dialog.update();
                if (flag.creating && response.created_flag_id) {
                    let newFlagId = response.created_flag_id;
                    let newFlag = flag.db.flags.get(newFlagId);
                    this.dialog.open(newFlag, {created: true});
                } else {
                    //this.dialog.open(flag);
                    this.update();
                }
                if (close) {
                    window.setTimeout(() => {
                        this.dialog.back();
                    },1);
                }
            },

            updateDescription() {
                let description = this.body.find('.description');
                description.empty();
                let flag = this.flag;

                description.append(flag.flagTypeObj().dom({flag: flag}));
                $('<div>', {class: 'flag-help', html: [
                    flag.helpFlagType(),
                    flag.helpSpecific(),
                    flag.statusDom()
                ]}).appendTo(description);
                
                if (flag.creating) {
                    $('<h5>', { text: 'Enter a comment and click save to raise this flag' }).appendTo(description);
                }
                return description;
            },

            scrollToBottom() {
                window.setTimeout(() => {
                    $(this.scrollablePanel).get(0).scrollTop = 10000000;
                },1);
            },

            updateComments() {
                // fixme, don't re-render comments already rendered
                shouldScroll = (isScrolledToBottom(this.scrollablePanel) && !this.firstUpdate) || this.justPosted;
                
                let lastEntry = this.lastEntry;
                let reachedStart = lastEntry == null;
                let commentsDom = this.body.find('.flag-comments');
                for (let comment of this.flag.comments()) {
                    if (!reachedStart) {
                        reachedStart = comment.id == lastEntry.id;
                        continue;
                    }
                    comment.dom({ lastEntry: lastEntry, includeTitle: false, dialog: this.dialog }).appendTo(commentsDom);
                    lastEntry = comment;
                }
                this.lastEntry = lastEntry;
                if (shouldScroll) {
                    this.scrollToBottom();
                }
            },

            setPosting() {
                this.footer.LoadingOverlay('show');
                this.justPosted = true;
            },
            
            currentResolution() {
                if (this.flag.creating) {
                    return 'XX';
                } else {
                    return this.flag.resolution;
                }
            },

            updateAction() {
                let flag = this.flag;
                let flagType = flag.flagTypeObj();
                let creating = flag.creating;                
                let newRes = this.currentResolution();
                
                this.footer.LoadingOverlay('hide');
                
                if (newRes == this.cachedState) {
                    // State has not changed, not updating
                    return;
                } else {
                    this.cachedState = newRes;
                }
                let resolutions = flag.resolutions();
                if (flag.creating && resolutions.length > 0) {
                    newRes = resolutions[0].id;
                }
                
                let commentsEnabled = flagType.comments_enabled;
                this.textarea.hide();
                if (commentsEnabled) {
                    this.textarea.show();
                }

                this.states.empty();
                if (resolutions.length > 0) {
                    let content = $('<div>', {class:'inline-form'}).appendTo(this.states);
                    for (let resolution of resolutions) {
                        let did = `res-${resolution.id}`;
                        let input = $('<input>', {class:"form-check-input", type:"radio", name:"flag-res", id:did, value:resolution.id});
                        if (resolution.id === newRes) {
                            input.attr('checked', true);
                        }
                        $('<div>', {class:'form-check form-check-inline', html:[
                            input,
                            $('<label>', {class:"form-check-label", 'for':did, text:resolution.label})
                        ]}).appendTo(content);
                    }

                    /*
                    $('<legend>', {class: 'status', text: flag.statusText()}).appendTo(content);
                    if (resolutions.length <= 4) {
                        let inputs = [];
                        for (let resolution of resolutions) {
                            let did = `res-${resolution.id}`;
                            $('<label>', {'for': did, text: resolution.label}).appendTo(content);
                            let input = $('<input>', {type: 'radio', name: 'flag-res', value: resolution.id, id: did}).appendTo(content);
                            inputs.push(input);
                            if (resolution.id == newRes) {
                                input.attr('checked', true);
                            }
                        }
                        $(inputs).checkboxradio();


                    } else {
                        let select = $('<select>', {name: 'flag-res'}).appendTo(content);
                        for (let resolution of resolutions) {
                            let opt = $('<option>', {value: resolution.id, text: resolution.label}).appendTo(select);
                            if (resolution.id == newRes) {
                                opt.attr('selected', true);
                            }
                        }
                    }

                     */
                }
                
                let canDoSomething = resolutions.length > 1 || commentsEnabled;
                if (canDoSomething) {
                    this.newActionDom.show();
                    this.footer.show();
                } else {
                    this.newActionDom.hide();
                    this.footer.hide();
                }
            },

            update() {
                this.updateDescription();
                this.updateComments();
                this.updateAction();
                this.firstUpdate = false;
                this.justPosted = false;
            }
        };
        return FlagContent;
    })();

    let FlagCollectionSummaryContent = (function() {
        let FlagCollectionSummaryContent = function (collection, dialog) {
            this.collection = collection;
            this.dialog = dialog;
            
            this.title = null;
            this.body = null;
            this.footer = null;
        };
        
        FlagCollectionSummaryContent.prototype = {
            flagIcon(flag, status) {
                let flagType = flag.flagTypeObj();
                let bigIcon = $('<div>', {class:`big-icon`, click: () => {
                    this.dialog.open(flag);
                }});
                $('<div>', {class:`flag flag-${flagType.id} res-${flag.resolution}`}).appendTo(bigIcon);
                $('<label>', {text: flagType.label}).appendTo(bigIcon);
                $('<div>', {class: 'indicator'}).appendTo(bigIcon);
                if (status) {
                    bigIcon.addClass(status);
                } else {
                    if (flag.creating) {
                        bigIcon.addClass('new');
                    } else if (flag.status == 'O') {
                        bigIcon.addClass('open');
                    } else if (flag.status == 'C') {
                        bigIcon.addClass('closed');
                    } else if (flag.status == 'R') {
                        bigIcon.addClass('rejected');
                    }
                }
                if (!flag.creating) {
                    let ageDom = createTimestampDom(flag.created, true);
                    ageDom.addClass('age');
                    ageDom.appendTo(bigIcon);
                }
                
                return bigIcon;
            },
            
            flagChunk(flags, title, status, showDetailed) {
                let chunk = $('<div>', {class: 'containe container-flags'});
                $('<h5>', {text: title}).appendTo(chunk);
                
                if (showDetailed) {
                    let rows = chunk;
                    
                    let lastType = null;
                    
                    for (let flag of flags) {
                        let row = $('<div>', {class: 'd-flex flex-row align-content-start mt-4 flag-detail-row'}).appendTo(rows);
                        this.flagIcon(flag, status).appendTo(row);
                        let help = $('<div>', {style: 'flex-grow:1', class: 'ml-2'}).appendTo(row);
                        
                        let firstOfType = false;
                        if (flag.flag_type !== lastType) {
                            lastType = flag.flag_type;
                            // provide help for all flags of this type
                            flag.helpFlagType().appendTo(help);
                        }
                        flag.helpSpecific().appendTo(help);
                        if (status == 'open' && flag.resolution) {
                            let res = flag.resolutionObj();
                            $('<div>', {html: [
                                $('<br>'),
                                $('<span>', {text: 'Current status : '}),
                                $('<b>', {text: res.label})
                            ]}).appendTo(help);
                        }
                    }
                    return chunk;
                } else {
                    for (let flag of flags) {
                        this.flagIcon(flag, status).appendTo(chunk);
                    }
                }
                return chunk;
            },
        
            init() {
                let dom = $('<div>');
                let owner = this.collection.user_permission >= 1;

                let openFlags = this.collection.flags();
                if (openFlags.length) {
                    this.flagChunk(openFlags, 'In Progress Flags', 'open', owner).appendTo(dom);
                }
                
                let resolvedFlags = this.collection.importantResolvedFlags();
                if (resolvedFlags.length) {
                    this.flagChunk(resolvedFlags, 'Resolved Flags', false).appendTo(dom);
                }

                if (!openFlags.length && !resolvedFlags.length) {
                    let chunk = $('<div>', {class: 'flag-chunk'}).appendTo(dom);
                    $('<h3>', {text: 'No Flags'}).appendTo(chunk);
                    $('<p>', {text: 'There are no In Progress Flags or important Resolved Flags attached to this record.'}).appendTo(chunk);
                }

                let raisableFlags = this.collection.raisableFlags();
                if (raisableFlags.length) {
                    this.flagChunk(raisableFlags, 'Raise New Flags', 'new', owner).appendTo(dom);
                }

                // this is mostly redundant to Variant Classification History
                let timelineA = $('<a>', {text: 'See Full Timeline', class:'btn btn-outline-secondary'});
                timelineA.click(() => {
                    this.dialog.timeline(); 
                });
                this.footer = timelineA;

                this.body = dom;
            },
            update() {}
        };
        
        return FlagCollectionSummaryContent;
    })();
    
    let FlagCollectionDialog = (function() {
        let FlagCollectionDialog = function(collection) {
            this.collection = collection;
            this.dialog = null;
            this.content = null;
            this.activeContent = null;
        };
        FlagCollectionDialog.prototype = {
            update(params) {
                if (this.activeContent) {
                    this.activeContent.update();
                }
            },
            async init(params) {
                if (this.collection.db.hasOpenDialog) {
                    console.log('A dialog is already open/opening');
                    return;
                }
                params = params || {};
                if (params.triggerDom) {
                    $(params.triggerDom).LoadingOverlay('show');
                }
                this.collection.db.hasOpenDialog = true;
                
                params = params || {};
                await this.collection.db.refresh({ collection: this.collection });
                let activeFlag = params.activeFlag;
                
                let modalContent = $(`
                    <div class="modal fade" id="FlagModal" tabindex="-1" role="dialog" aria-labelledby="FlagModalLabel" aria-hidden="true">
                        <div class="modal-dialog modal-lg" role="document">
                            <div class="modal-content">
                                <div class="modal-header">
                                </div>
                                <div class="modal-body">
                                </div>
                                <div class="modal-footer">
                                </div>
                            </div>
                        </div>
                    </div>
                `);
                modalContent.find('.modal-header').html(`
                    <nav aria-label="breadcrumb">
                      <ol class="breadcrumb"></ol>
                    </nav>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                        <span aria-hidden="true">&times;</span>
                    </button>
                `);

                this.content = modalContent;

                let db = this.collection.db;
                let dialog = this;
                modalContent.on('shown.bs.modal', function() {
                    db.setDialogState(dialog);
                });
                modalContent.on('hidden.bs.modal', function() {
                    db.setDialogState(null);
                    if (db.onClose) {
                        db.onClose({addedOrClosedFlag: db.addedOrClosedFlag});
                    }
                    modalContent.modal('dispose');
                    modalContent.remove();
                });

                if (activeFlag) {
                    this.open(activeFlag);
                } else {
                    this.applyActiveContent(new FlagCollectionSummaryContent(this.collection, this))
                }
                modalDialog = modalContent.modal({focus:true, show:true});
                if (params.triggerDom) {
                    $(params.triggerDom).LoadingOverlay('hide');
                }
            },
            open(flag) {
                this.applyActiveContent(new FlagContent(flag, this));
            },

            applyActiveContent(activeContent) {
                activeContent.init();
                activeContent.update();
                this.activeContent = activeContent;
                let parentTitle = (this.collection.label || 'Flags').replace('/', '-'); // fixme trim title
                let title = this.activeContent.title;
                let breadcrumbs = this.content.find('ol.breadcrumb');
                if (!title) {
                    breadcrumbs.html([
                        $('<li>', {class: 'breadcrumb-item', text: parentTitle})
                    ]);
                } else {
                    breadcrumbs.html([
                        $('<li>', {class: 'breadcrumb-item', html: $('<a>', {text: parentTitle, click: () => { this.back() }})}),
                        $('<li>', {class: 'breadcrumb-item active', text: title})
                    ]);
                }
                this.content.find('.modal-body').html(this.activeContent.body);
                this.content.find('.modal-footer').html(this.activeContent.footer);
            },

            timeline() {
                this.applyActiveContent(new FlagTimeline(this.collection, this));
            },
            
            back() {
                this.applyActiveContent(new FlagCollectionSummaryContent(this.collection, this));
            },
            
        };
        return FlagCollectionDialog;
    })();

    let Flag = (function () {
        let Flag = function (db, data) {
            this.db = db;
        };
        Flag.prototype = {
            updated() { this.open = this.status == 'O'; },
            userObj() { return this.db.users.get(this.user); },
            flagTypeObj() { return this.db.flagTypes.get(this.flag_type); },
            collectionObj() { return this.db.collections.get(this.collection); },
            title() { return `${this.flagTypeObj().label}` }, //  (${this.id})
            resolutionObj() { return this.db.flagResolutions.get(this.resolution); },
            canEdit() {
                let userPermission = this.collectionObj().user_permission;
                let requiredPermission = this.flagTypeObj().permission;
                return userPermission >= requiredPermission;
            },
            subFlags() {
                let extraData = this.extra_data;
                if (extraData) {
                    let subFlags = extraData.sub_flags;
                    if (subFlags) {
                        return subFlags;
                    } else {
                        return [];
                    }
                }
            },
            
            statusText() {
                if (this.flag_type == 'classification_significance_change') {
                    return 'Primary Reason for Change';
                }
                if (this.creating) {
                    return 'Opening Status';
                }
                return 'Status';
            },
            
            resolutions() {
                let resolutions = this.flagTypeObj().resolutions.map(r => this.db.flagResolutions.get(r) );
                
                if (this.creating) {
                    let excludeStatuses = {'R': true};
                    if (!this.canEdit()) {
                        excludeStatuses['C'] = true;
                    }
                    return resolutions.filter(r => !excludeStatuses[r.status]); 
                } else if (!this.canEdit()) {
                    return [];
                } else {
                    return resolutions;
                }
            },
            
            helpFlagType() {
                if (this.creating) {
                    if (this.flag_type === 'classification_suggestion') {
                        return $(`<div>
                        If you have found some extra information that you think should be incorporated into this classification,
                        you can raise a suggestion for the classification owner to accept or reject.
                        </div>`);
                    } else if (this.flag_type === 'classification_internal_review') {
                        return $(`<div>
                        You can raise this flag to let people know the classification is currently in review, or raise it
                        as "Completed" to record the fact that a review has recently taken place.<br/>
                        Please record any internal reviews while a classification is marked as discordant.
                        </div>`);
                    }
                    return $('<div>');
                }
            
                if (this.flag_type === 'classification_suggestion') {
                    return $(`<div>
                    Someone has raised suggestion(s) against this classification.
                    <ol><li>Review the contents of each suggestion.</li>
                    <li>If appropriate, make changes in your curation system and mark the suggestion as Complete.</li>
                    <li>If you decline the suggestion, mark it as Rejected.</li>
                    </ol></div>
                    `);
                } else if (this.flag_type === 'classification_outstanding_edits') {
                    return $(`<div>
                    Edits have been made to this classification that are not included in a published version.
                    <ol><li>From the classification form, ensure there are no validation errors stopping this record from being published.</li>
                    <li>At the bottom of the form, click the tick to submit the outstanding changes.</li></ol></div>`
                    );
                } else if (this.flag_type === 'classification_matching_variant') {
                    return $(`<div>
                    This classification is not yet linked to a variant
                    <ol><li>If this has a status of In Progress we should match it to a variant shortly, no action required.</li>
                    <li>If this has a status of Matching Failed we were unable to normalise the variant provided based on the c.hgvs and genome build values.
                    Please contact Shariant support for help in resolving this.</li></ol></div>
                    `);
                } else if (this.flag_type === 'classification_matching_variant_warning') {
                    return $(`<div>
                    This classification has been matched to a variant, but requires a manual check to ensure it was matched correctly
                    <ol><li>If you believe it was matched correctly, select a status of Variant Confirmed and Save</li>
                    <li>If you believe that this is not the variant intended, select a status of Variant Rejected and Save. An admin will then attempt to fix the problem.</li>
                    </ol></div>
                    `);
                 } else if (this.flag_type === 'classification_transcript_version_change') {
                    return $(`<div>
                    This classification has been matched to a variant, but we did not have the requested transcript version on file and had to migrate it to another version of the same transcript.
                    <ol><li>If you believe the change of transcript version in this case has no impact on the variant, select a status of Change Accepted and Save.</li>
                    <li>If you believe this change is significant and that this classification should not apply to the variant on the new transcript, select a status of Change Rejected and Save. An admin will then attempt to fix the problem.</li>
                    </ol></div>
                    `);
                } else if (this.flag_type === 'classification_internal_review') {
                    return $(`<div>
                    This classification is marked as currently being internally reviewed.
                    <ol><li>Once the internal review is complete, ensure you update the classification in your curation system.</li>
                    <li>Mark the internal review as Completed</li></ol>
                    </div>
                    `);
                } else if (this.flag_type === 'classification_withdrawn') {
                    return $(`<div>
                    This classification has been marked as withdrawn. It will be hidden from almost all searches and exports.
                    <ol><li>If the classification is not of high enough quality or in error, you may leave it as "withdrawn" indefinitely.</li>
                    <li>If you wish to un-withdraw the classification, click the open bin icon in actions from the variant classification form</li></ol></div>
                    `);
                } else if (this.flag_type === 'classification_significance_change') {
                    return $(`<div>
                    This classification has changed its clinical significance compared to a previously published version.
                    <ol><li>Set the status of this flag to reflect the primary reason behind the change in classification</li>
                    <div><ul>
                    <li>Discordance Discussion - Data was changed as a result of talking to other labs when this classification was in discordance.
                    <li>Summation of Data - Data was changed as a result of combining information from multiple labs.
                    <li>Internal Review - New data or errors in old data were found during an internal review.
                    </ul>                 
                    </div>
                    <li>Please also add a comment providing some context.</li></ol></div>`);
                } else if (this.flag_type === 'classification_discordant') {
                    return $(`<div>
                    This classification is in discordance with one or more classifications.
                    <ol><li>Ensure that you have completed an internal review of your lab's classification recently (within the last 12 months is recommended). If not, raise the internal review flag and complete an internal review of your lab's classification.
                    <li>Review any outstanding suggestions against your lab's classification.
                    <li>View the other classifications in the discordance report and view the evidence differing between multiple records via the diff page. If appropriate, raise suggestions against other lab classifications.
                    <li>This Discordance flag will automatically be closed when concordance is reached.
                    </ol></div>
                    `);
                } else if (this.flag_type === 'classification_unshared') {
                    return $(`<div>
                    This classification is not yet shared outside of your lab or institution.
                    <ol><li>From the classification form, ensure there are no validation errors stopping this record from being published.</li>
                    <li>Review the content of the classification to make sure it's ready to be shared.</li>
                    <li>At the bottom of the form, click the Share to submit at a higher share level.</li></ol></div>
                    `);
                } else if (this.flag_type === 'allele_missing_37') {
                    return $(`<div>
                    This allele is missing a variant for GRCh37 and can't be lifted over in exports.<br/>
                    Our admins are looking at this allele and working on a fix.
                    </div>`);
                } else if (this.flag_type === 'allele_missing_38') {
                    return $(`<div>
                    This allele is missing a variant for GRCh38 and can't be lifted over in exports.<br/>
                    Our admins are looking at this allele and working on a fix.
                    </div>`);
                } else if (this.flag_type === 'allele_37_not_38') {
                    return $(`<div>
                    For a given transcript, the c.hgvs produced for this allele produces a different value for GRCh37 to GRCh38.<br/>
                    Our admins are looking at this allele and working on a fix.
                    </div>`);
                }
                return $('<div>');
            },
            
            statusDom() {
                if (this.creating) {
                    return $('<div>');
                }
                return $('<div>', {html: [
                    `<br/>This ${this.flagTypeObj().label} has the status of `,
                    $('<b>', {text: this.resolutionObj().label})
                ]});
            },
            
            helpSpecific() {
                if (this.creating) {
                    return $('<div>');
                }
                let firstComment = this.comments()[0].text;
                if (firstComment.length > 100) {
                    firstComment = firstComment.substring(0, 97) + '...';
                }
                
                if (this.data && this.data.transcript) {
                    let transcript = this.data.transcript;
                    parts = /^([_A-Z0-9]+)(?:[.]([0-9]+))?$/i;
                    let match = parts.exec(transcript);
                    if (match) {
                        let transcriptUrl = null;
                        let transcriptNoVer = match[1];
                        transcriptUrl = `/genes/view_transcript/${transcriptNoVer}`;
                        return $('<div>', {html: [
                            'View details about the transcript ',
                            $('<a>', {href:transcriptUrl, class:'hover-link', text:`${transcriptNoVer}`, target:'_blank'}),
                            '.'
                        ]});
                    }
                } else if (this.flag_type === 'classification_suggestion') {
                    return $('<div>', {html: [
                        `Raised by `,
                        $('<span>', {class: 'username', text:this.userObj().name}),
                        ` in regards to `,
                        $('<span>', {class: 'quote', text: firstComment})
                    ]});
                } else if (this.flag_type === 'classification_significance_change') {
                    return $('<div>', {text: firstComment});
                } else if (this.flag_type === 'classification_matching_variant_warning') {
                    let variantId = this.collectionObj().variant;
                    return $('<div>', {html: [
                        `See more information about the linked variant `,
                        $('<a>', {class: 'hover-link', text: `here`, href:`/variantdetails/view_allele_from_variant/${variantId}`})
                    ]});
                } else if (this.flag_type === 'classification_discordant') {
                    let reportId = this.collectionObj().discordance_report;
                    let variantId = this.collectionObj().variant;
                    let clinicalContext = this.collectionObj().clinical_context;
                    if (reportId) {
                        return $('<div>', {html: [
                            `Go to the `,
                            $('<a>', {class: 'hover-link', text: `Discordance Report`, href:`/classification/classification/discordance_report/${reportId}`}),
                            ` | `,
                            $('<a>', {class: 'hover-link', text: `Diff with other Classifications`, href: `/classification/diff/?clinical_context=${clinicalContext}`})
                        ]});
                    }
                }
                return $('<div>');
            },
            
            ageText() {
                return jQuery.timeago(this.created * 1000) + ' ago';
            },
            comments() {
                return Object.values(this.db.comments.all()).filter(c => c.flag == this.id).sort((c1, c2) => c1.created - c2.created);
            },
            dom() {
                let flagType = this.flagTypeObj();
                let user = this.userObj();

                let titleText = `${flagType.label} (${this.ageText()})`;

                let collectionObj = this.collectionObj();
                let flagDiv = $('<div>', { class: `flag flag-${flagType.id} res-${this.resolution}`, title: titleText});
                flagDiv.click(() => {
                    flagDiv.tooltip('hide');
                    new FlagCollectionDialog(this.collectionObj()).init({activeFlag:this, triggerDom:flagDiv})
                });
                if (this.open == false) {
                    flagDiv.addClass('closed');
                }
                return flagDiv;
            },
            remove() {
                for (let comment of this.comments()) {
                    comment.remove();
                }
                this.db.flags._remove(this.id);
            },
        };
        return Flag;
    })();


    let FlagComment = (function () {
        let FlagComment = function (db) {
            this.db = db;
        };
        FlagComment.prototype = {
            flagObj() { return this.db.flags.get(this.flag); },
            userObj() { return this.db.users.get(this.user); },
            resolutionObj() { return this.resolution ? this.db.flagResolutions.get(this.resolution) : null; },
            action() {
                let resolution = this.resolutionObj();
                if (!resolution) {
                    return 'Commented'
                } else {
                    return resolution.label;
                }
            },
            title() {
                return `
                    \n${this.action()} ${moment(this.created * 1000).format('DD-MMM-YYYY HH:mm')}
                    \n${this.userObj().name} : ${this.text || this.summary || '&lt;no comment provided&gt;'}
                `;
            },
            dom(params) {
                let { includeTitle, lastEntry, dialog } = params || {};

                let flag = this.flagObj();
                let user = this.userObj();
                let flagType = flag.flagTypeObj();
                let title = `${flagType.label} (${flag.ageText()})`;
                let commentDom = $('<div>', { class: 'flag-comment' });

                let timestamp = moment(this.created * 1000).format('DD-MMM-YYYY');
                if (lastEntry) {
                    oldTimestamp = moment(lastEntry.created * 1000).format('DD-MMM-YYYY');
                    if (oldTimestamp == timestamp) {
                        timestamp = null;
                    }
                }
                if (timestamp) {
                    $('<div>', { class: 'timestamp flag-title', text: timestamp }).appendTo(commentDom);
                    commentDom.addClass('with-title');
                }

                let userActionDom = $('<div>', { class: 'user-action' }).appendTo(commentDom);

                let avatar = $('<div>', { class: 'avatar' }).css(`background-image`, `url(${user.avatar})`).css('background-color', user.color).appendTo(userActionDom);
                let content = $('<div>', { class: 'content' }).appendTo(userActionDom);
                let time = moment(this.created * 1000).format('hh:mm A');
                let userInfo = $('<div>', {
                    class: 'header', html: [
                        $('<span>', { class: 'username', text: user.name }),
                        $('<span>', { class: 'time', text: time })
                    ]
                }).appendTo(content);

                if (!timestamp && lastEntry && lastEntry.user == this.user && time == moment(lastEntry.created * 1000).format('hh:mm A')) {
                    avatar.css('height', '0').css('visibility', 'hidden');
                    userInfo.css('display', 'none');
                }

                let actionCommentDom = $('<div>', { class: 'action-text' }).appendTo(content);
                action = this.action();

                if (params.includeTitle) {
                    $('<div>', {
                        class: `flag-link hover-link ${flag.open ? '' : 'closed'}`,
                        click: () => dialog.open(flag),
                        html: [
                            $('<div>', { class: `flag flag-${flagType.id} res-${flag.resolution}` }),
                            $('<span>', { text: flag.title() })
                        ]
                    }).appendTo(actionCommentDom);
                }
                $('<span>', {class:'action', text: `${action}`}).appendTo(actionCommentDom);

                if (this.text && this.text.length) {
                    $('<div>', { class: 'text', text: this.text }).appendTo(content);
                }

                return commentDom;
            },
            remove() {
                this.db.comments._remove(this.id);
            }
        };
        return FlagComment;
    })();


    let Flags = function (
        props
    ) {
        this.flagsUrl = '/flags/api/flags/';
        this.flagUrl = '/flags/api/flag/';

        this.userId = null;
        this.users = new DataCollection(() => new User(this));
        this.flagResolutions = new DataCollection(() => new FlagResolution(this));
        this.flagTypes = new DataCollection(() => new FlagType(this));
        this.flags = new DataCollection(() => new Flag(this));
        this.collections = new DataCollection(() => new FlagCollection(this));
        this.comments = new DataCollection(() => new FlagComment(this));
        this.dialogState = null;
        this.hasOpenDialog = false;
        this.refreshToken = this.refreshCurrentDialog();
        this.summaryFilter = null;
        this.summaryFilterInclusions = null;
        this.filterValue = null;
        
        this.addedOrClosedFlag = false;
    };
    
    Flags.prototype = {

        updateFilter() {
            if (this.summaryFilter.length) {
                this.summaryFilter.empty();
                let counts = {};
                for (let flag of this.flags.all()) {
                    if (flag.open) {
                        let resolution_counts = counts[flag.flag_type] || {};
                        resolution_counts[flag.resolution] = (resolution_counts[flag.resolution] || 0) + 1;
                        counts[flag.flag_type] = resolution_counts;
                    }
                }
                let ordered_keys = Object.keys(counts).sort((k1, k2) => this.flagTypes.get(k1).label.localeCompare(this.flagTypes.get(k2).label));
                for (let key of ordered_keys) {
                    
                    if (this.summaryFilterInclusions && !this.summaryFilterInclusions[key]) {
                        continue;
                    }
                    
                    let resolution_counts = counts[key];
                    let flagType = this.flagTypes.get(key);
                    for (let [resolution, value] of Object.entries(resolution_counts)) {

                        let label =  $('<span>', {class: 'label', text: flagType.label});
                        if (resolution !== 'open') {
                            let resolutionObj = this.flagResolutions.get(resolution);
                            $('<span>', {text: ` (${resolutionObj.label})`}).appendTo(label);
                        }
                        
                        let cell_parts = [
                            $('<div>').addClass('flag').addClass(`flag-${flagType.id} res-${resolution}`),
                            $('<div>', {class: 'text-monospace mx-2', text: value.toString(), style:'min-width: 2rem; text-align:right'}),
                            label
                        ];
                        if (value > 1 && !flagType.label.endsWith('s')) {
                            $('<span>', {class: 'plural', text: 's'}).appendTo(label);
                        }
                        
                        let filterCell = $('<a/>', {
                            class: 'list-group-item list-group-item-action d-flex',
                            'data-flag-type': `${flagType.id}-${resolution}`,
                            html: cell_parts,
                            click: () => {
                                this.applyFilter(flagType.id, resolution);
                            },
                            //title: flagType.description
                        }).appendTo(this.summaryFilter);
                    }
                }
            }
        },
        
        applyFilter(flag_type_id, resolution) {
            let filterValue = `${flag_type_id}-${resolution}`;
            this.summaryFilter.find('.list-group-item').removeClass('list-group-item-primary');
            if (this.filterValue === filterValue || flag_type_id === null) {
                this.filterValue = null;
                this.collections.all().forEach(c => c.dom.closest('tr').show());
            } else {
                this.filterValue = filterValue;
                this.summaryFilter.find(`[data-flag-type=${flag_type_id}-${resolution}]`).addClass('list-group-item-primary');
                let matching_collection_ids = {};
                let flags = this.flags.all().filter(f => f.flag_type === flag_type_id && f.resolution === resolution).forEach(f => {
                    matching_collection_ids[f.collection] = true;
                });
                for (let fc of this.collections.all()) {
                    if (matching_collection_ids[fc.id]) {
                        fc.dom.closest('tr').show();
                    } else {
                        fc.dom.closest('tr').hide();
                    }
                }
            }
        },

        setDialogState(collectionDialog) {
            this.dialogState = collectionDialog;
            if (!collectionDialog) {
                this.hasOpenDialog = false;
            }
        },

        async refreshCurrentDialog() {
            if (this.dialogState && this.dialogState) {
                await this.refresh();
                this.render();
                if (this.dialogState) {
                    this.dialogState.update();
                }
            }
            return window.setTimeout(() => {
                this.refreshCurrentDialog();
            }, 5000);
        },

        init(props) {
            props = props || {};
            if (props.userId) {
                this.userId = props.userId;
            }
            
            this.onClose = props.onClose || false;
            this.summaryFilter = $('#flags-filter');
            if (this.summaryFilter) {
                summaryFilterInclusions = this.summaryFilter.attr('data-filter') || null;
                if (summaryFilterInclusions) {
                    this.summaryFilterInclusions = {};
                    for (let inclusion of summaryFilterInclusions.split(' ')) {
                        this.summaryFilterInclusions[inclusion.trim()] = true;
                    }
                }
            }

            let changes = false;
            let existingCollections = {};
            for (let collection of this.collections.all()) {
                existingCollections[collection.id] = collection;
            }
            let flag_collection_doms = {};
            $('[data-flags]').toArray().forEach(e => {
                let dom = $(e);
                dom.addClass('flags');
                let id = parseInt(dom.attr('data-flags'));

                if (!isNaN(id)) {
                    flag_collection_doms[id] = Object.assign({ dom: dom }, flag_collection_doms[id] || {});
                }
            });
            for (let entry of Object.entries(flag_collection_doms)) {
                let id = entry[0];
                let data = entry[1];
                if (existingCollections[id]) {
                    delete existingCollections[id];
                } else {
                    changes = true;
                }

                this.collections.upsert({
                    id: id,
                    dom: data.dom
                });
            }

            for (let collection of Object.values(existingCollections)) {
                collection.remove();
            }
            if (changes) {
                this.refresh();
            } else {
                this.render();
            }
        },

        post(data, url) {
            return new Promise((resolve, reject) => {
                $.ajax({
                    headers: {
                        'Accept': 'application/json',
                        'Content-Type': 'application/json'
                    },
                    data: JSON.stringify(data),
                    url: url,
                    type: 'POST',
                    error: (call, status, text) => {
                        reject(status);
                    },
                    success: async (record) => {
                        this.updateData(record);
                        resolve(record);
                    }
                });
            });
        },

        currentUser() {
            return this.users.get(this.userId);
        },

        refresh(params) {
            return new Promise((resolve, reject) => {
                let url = null;
                let rerender = false;
                params = params || {};

                rerender = true;
                let collections = this.collections.all();
                if (collections.length == 0) {
                    resolve();
                    return;
                }
                collection = params.collection;
                if (collection) {
                    params['history'] = collection.id;
                } else if (this.dialogState && this.dialogState.collection) {
                    collection = this.dialogState.collection;
                }
                if (collection) {
                    url = this.flagsUrl + collection.id;
                    if (collection.since) {
                        params.since = collection.since;
                    }
                } else {
                    url = this.flagsUrl + this.collections.all().filter(c => !!c.dom).map(c => c.id).join(',');
                }
                
                params = Object.assign({}, params);
                delete params['collection'];

                $.ajax({
                    headers: {
                        'Accept': 'application/json',
                        'Content-Type': 'application/json'
                    },
                    url: url,
                    data: params || null,
                    type: 'GET',
                    error: (call, status, text) => {
                        reject(status);
                    },
                    success: (record) => {
                        if (record.since && collection) {
                            collection.since = record.since;
                        }
                        this.updateData(record);
                        if (rerender) {
                            this.render();
                        }
                        this.updateFilter();
                        resolve();
                    }
                });
            });
        },

        updateData(record) {
            for (let user of (record.users || [])) {
                this.users.upsert(user);
            }
            for (let flag of (record.flags || [])) {
                this.flags.upsert(flag);
            }
            for (let comment of (record.comments || [])) {
                this.comments.upsert(comment);
            }
            for (let flagResolution of (record.flag_resolutions || [])) {
                this.flagResolutions.upsert(flagResolution);
            }
            for (let flagType of (record.flag_types || [])) {
                this.flagTypes.upsert(flagType);
            }
            for (let collection of (record.collections || [])) {
                this.collections.upsert(collection);
            }
        },

        render() {
            for (let collection of this.collections.all()) {
                let dom = collection.dom;
                if (!dom) {
                    continue;
                }
                dom.empty();
                let flags = collection.flags();
                let watch = collection.watching === 0 || collection.watching;
                let flagSummary = $('<div>', { class: `flag add`, title: 'Add or review flags for ' + collection.label}).appendTo(dom);
                flagSummary.click(() => { new FlagCollectionDialog(collection).init({triggerDom: flagSummary}) });

                if (collection.watching) {
                    $('<div>', { class: 'notifications', title: `${collection.watching} unseen activities` }).appendTo(flagSummary);
                } else if (watch) {
                    $('<div>', { class: `notifications watching`, title: `starred` }).appendTo(flagSummary);
                }
                
                for (let flag of collection.flags()) {
                    let flagDom = flag.dom().appendTo(dom);
                    dom.append(flagDom);
                    let subFlags = flag.subFlags();
                    if (subFlags && subFlags.length) {
                        let subFlagsDom = $('<div>', {class: 'sub-flags'}).appendTo(dom);
                        for (let subFlag of subFlags) {
                            let subFlagDom = $('<div>', {class: `sub-flag`, title: subFlag.label, text: subFlag.letter || 'X'});
                            subFlagDom.appendTo(subFlagsDom);
                        }
                    }
                }
            }
        }
    };

    Flags.instance = new Flags();

    return Flags;

})();