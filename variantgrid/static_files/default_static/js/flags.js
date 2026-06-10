function isScrolledToBottom(scrollDiv) {
    if (scrollDiv) {
        scrollDiv = $(scrollDiv).get(0);
        const scrolledUpBy = scrollDiv.scrollHeight - scrollDiv.scrollTop - scrollDiv.clientHeight;
        if (scrolledUpBy < 1) {
            return true;
        }
    }
    return false;
}

const Flags = (function () {

    const DataCollection = (function () {
        const DataCollection = function (initer) {
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
                const existing = this.dataMap[id];
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

    const User = (function () {
        const User = function (db) {
            this.db = db;
        };
        User.prototype = {
            dom() {
                return $('<div>', {
                    class: 'user', html: [
                        $('<div>', { class: 'avatar' }).css(`background-image`, `url(${this.avatar})`).css('background-color', this.color),
                        $('<div>', { class: 'username', text: this.name }),
                    ]
                });
            }
        };
        return User;
    })();
    
    const FlagResolution = (function() {
        const FlagResolution = function(db) {
            this.db = db;
        };
        FlagResolution.prototype = {};
        return FlagResolution;
    })();

    const FlagType = (function () {
    
        const FlagType = function (db) {
            this.db = db;
        };
        FlagType.prototype = {
            collectionObj() { return this.db.collections.get(this.collection) },

            dom(params) {
                params = params || {};
                const {flag, statusFlag} = params;

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

    const FlagCollection = (function () {
        const FlagCollection = function (db) {
            this.db = db;
        };
        FlagCollection.prototype = {
            comments() {
                return this.db.comments.all().filter(c => c.flagObj().collection == this.id).sort((c1, c2) => c1.created - c2.created);
            },
            commentsImportant() {
                return this.db.comments.all().filter(c => {
                    const flag = c.flagObj();
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
                        const localeCompare = -f2.flagTypeObj().label.localeCompare(f1.flagTypeObj().label);
                        return localeCompare ? localeCompare : f2.created - f1.created;
                    });
            },
            importantResolvedFlags() {
                return this.db.flags.all().filter(f => f.collection == this.id && !f.open && f.flagTypeObj().importance >= 2)
                    .sort((f1, f2) => {
                        const localeCompare = -f2.flagTypeObj().label.localeCompare(f1.flagTypeObj().label);
                        return localeCompare ? localeCompare : f2.created - f1.created;
                    });
            },
            raisableFlags() {
                flagTypes = this.flagTypes().filter(ft => {
                    const userPermission = this.user_permission;
                    const requiredPermission = ft.raise_permission;
                    return userPermission >= requiredPermission;
                });
                return flagTypes.map(ft => this.createFlagFor(ft)).filter(f => f != null);
            },
            createFlagFor(flagType) {
                if (flagType.only_one) {
                    const existing = this.db.flags.all().find(f => f.collection == this.id && f.flag_type == flagType.id);
                    if (existing) {
                        if (existing.status != 'O') {
                            return existing;
                        }
                        return null;
                    }
                }
                const newFlag = new Flag(this.db);
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
                for (const flag of this.flags()) {
                    flag.remove(flag.id);
                }
                this.db.collections._remove(this.id)
            },

            recent_activity() {
                // FIXME, filter to actually be recent
                return this.comments();
            },

            async watch_toggle() {
                const watching = this.watching === 0 || this.watching;
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

    const FlagTimeline = (function() {
        const FlagTimeline = function(flagCollection, dialog) {
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
                const dom = $('<div>');
                this.body = dom;

                const content = $('<div>', { class: 'edit-flag-dialog container' }).appendTo(dom);
                const details = $('<div>', {class: 'flag-details'}).appendTo(content);
                this.scrollablePanel = details;
                
                $('<div>', { class: 'flag-comments' }).appendTo(details);
            },
            update() {
                // fixme, don't re-render comments already rendered
                const shouldScroll = (isScrolledToBottom(this.scrollablePanel) && !this.firstUpdate) || this.justPosted;

                let lastEntry = this.lastEntry;
                let reachedStart = lastEntry == null;
                const commentsDom = this.body.find('.flag-comments');
                const importantStuff = this.flagCollection.commentsImportant();

                if (!importantStuff.length) {
                    commentsDom.empty();
                    $('<p>', {text: 'There have been no important flag actions against this record.'}).appendTo(commentsDom);
                } else {
                    for (const comment of importantStuff) {
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

    const FlagContent = (function() {
        const FlagContent = function(flag, dialog) {
            this.flag = flag;
            this.dialog = dialog;
            
            this.scrollablePanel = null;
            this.textarea = null;
            this.states = null;

            this.title = null;
            this.titleDom = null;
            this.body = null;
            this.footer = null;

            this.cachedState = null;
            this.newActionDom = null;
            this.firstUpdate = true;
            this.justPosted = false;
        };
        FlagContent.prototype = {
            init() {
                const flag = this.flag;
                const flagType = flag.flagTypeObj();
                const flagDiv = $('<div>', { class: `flag flag-${flagType.id} res-${this.resolution} mr-1`});
                this.titleDom = [flagDiv, flagType.label];
                if (flag.open === false) {
                    flagDiv.addClass('closed');
                }
                if (flag.creating) {
                    this.titleDom.push(' (New)');
                } else {
                    this.titleDom.push(` (Raised: ${flag.ageText()})`);
                }
                const dom = $('<div>');
                this.body = dom;

                const content = $('<div>', { class: 'edit-flag-dialog' }).appendTo(dom);
                const details = $('<div>', {class: 'flag-details'}).appendTo(content);
                this.scrollablePanel = details;
                
                $('<div>', { class: 'description' }).appendTo(details);
                $('<hr/>').appendTo(details);
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
                const flag = this.flag;
            
                params = params || {};
                const { close } = params;
                let comment = null;
                if (this.textarea) {
                    comment = this.textarea.val().trim();
                    if (comment.length === 0) {
                        comment = null;
                    }
                    this.textarea.val('');
                }
                const resolution = $(`select[name='flag-res']`).val() || $(`input[name='flag-res']:checked`).val();
                
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

                const response = await flag.db.post(
                    sendParams,
                    url
                );
                await flag.db.refresh();
                
                this.dialog.update();
                if (flag.creating && response.created_flag_id) {
                    const newFlagId = response.created_flag_id;
                    const newFlag = flag.db.flags.get(newFlagId);
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
                const description = this.body.find('.description');
                description.empty();
                const flag = this.flag;

                $('<div>', {class: 'flag-help', html: [
                    flag.helpFlagType(),
                    flag.helpSpecific(),
                    flag.statusDom()
                ]}).appendTo(description);
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
                const commentsDom = this.body.find('.flag-comments');
                for (const comment of this.flag.comments()) {
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
                const flag = this.flag;
                const flagType = flag.flagTypeObj();
                const creating = flag.creating;                
                let newRes = this.currentResolution();
                
                this.footer.LoadingOverlay('hide');
                
                if (newRes == this.cachedState) {
                    // State has not changed, not updating
                    return;
                } else {
                    this.cachedState = newRes;
                }
                const resolutions = flag.resolutions();
                if (flag.creating && resolutions.length > 0) {
                    newRes = resolutions[0].id;
                }
                
                const commentsEnabled = flagType.comments_enabled;
                this.textarea.hide();
                if (commentsEnabled) {
                    this.textarea.show();
                }

                this.states.empty();
                if (resolutions.length > 0) {
                    if (flag.creating) {
                        $('<label>', { class:'mb-2', text: 'Enter a comment and click save to raise this flag' }).appendTo(this.states);
                    }
                    const content = $('<div>', {class:'d-flex mb-2'}).css('align-items','center').appendTo(this.states);
                    $('<label/>', {text: 'Status:', class:'mr-2 align-center'}).appendTo(content);
                    const statusButtons = $('<div>', {class:'btn-group btn-group-toggle', 'data-toggle':'buttons'}).appendTo(content)
                    for (const resolution of resolutions) {
                        const did = `res-${resolution.id}`;
                        const input = $('<input>', {type:"radio", name:"flag-res", id:did, value:resolution.id});
                        const label = $('<label>', {class:'btn btn-outline-secondary', html:[
                            input, resolution.label
                        ]});
                        if (resolution.id === newRes) {
                            input.attr('checked', true);
                            label.addClass('active');
                        }
                        label.appendTo(statusButtons);
                    }

                    /*
                    // old code for providing drop down if too many options
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
                
                const canDoSomething = resolutions.length > 1 || commentsEnabled;
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

    const FlagCollectionSummaryContent = (function() {
        const FlagCollectionSummaryContent = function (collection, dialog) {
            this.collection = collection;
            this.dialog = dialog;
            
            this.title = null;
            this.body = null;
            this.footer = null;
        };
        
        FlagCollectionSummaryContent.prototype = {
            flagIcon(flag, status) {
                const flagType = flag.flagTypeObj();
                const bigIcon = $('<div>', {class:`big-icon`, click: () => {
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
                    } else if (flag.status === 'O') {
                        bigIcon.addClass('open');
                    } else if (flag.status === 'C') {
                        bigIcon.addClass('closed');
                    } else if (flag.status === 'R') {
                        bigIcon.addClass('rejected');
                    }
                }
                if (!flag.creating) {
                    const ageDom = createTimestampDom(flag.created, true);
                    ageDom.addClass('age');
                    ageDom.appendTo(bigIcon);
                }
                
                return bigIcon;
            },
            
            flagChunk(flags, title, status, showDetailed) {
                const chunk = $('<div>', {class: 'containe container-flags'});
                $('<h5>', {text: title}).appendTo(chunk);
                
                if (showDetailed) {
                    const rows = chunk;
                    
                    let lastType = null;
                    
                    for (const flag of flags) {
                        const row = $('<div>', {class: 'd-flex flex-row align-content-start mt-4 flag-detail-row'}).appendTo(rows);
                        this.flagIcon(flag, status).appendTo(row);
                        const help = $('<div>', {style: 'flex-grow:1', class: 'ml-2'}).appendTo(row);
                        
                        const firstOfType = false;
                        if (flag.flag_type !== lastType) {
                            lastType = flag.flag_type;
                            // provide help for all flags of this type
                            flag.helpFlagType().appendTo(help);
                        }
                        flag.helpSpecific().appendTo(help);
                        if (status == 'open' && flag.resolution) {
                            const res = flag.resolutionObj();
                            $('<div>', {html: [
                                $('<span>', {text: 'Current status : '}),
                                $('<b>', {text: res.label})
                            ]}).appendTo(help);
                        }
                    }
                    return chunk;
                } else {
                    for (const flag of flags) {
                        this.flagIcon(flag, status).appendTo(chunk);
                    }
                }
                return chunk;
            },
        
            init() {
                const dom = $('<div>');
                const owner = this.collection.user_permission >= 1;

                const openFlags = this.collection.flags();
                if (openFlags.length) {
                    this.flagChunk(openFlags, 'In Progress Flags', 'open', owner).appendTo(dom);
                }
                
                const resolvedFlags = this.collection.importantResolvedFlags();
                const seenOnlyOneFlagTypes = {};
                if (resolvedFlags.length) {
                    this.flagChunk(resolvedFlags, 'Resolved Flags', false).appendTo(dom);
                    for (const rf of resolvedFlags) {
                        if (rf.db.flagTypes.get(rf.flag_type).only_one) {
                            seenOnlyOneFlagTypes[rf.flag_type] = true;
                        }
                    }
                }

                if (!openFlags.length && !resolvedFlags.length) {
                    const chunk = $('<div>', {class: 'flag-chunk'}).appendTo(dom);
                    $('<h3>', {text: 'No Flags'}).appendTo(chunk);
                    $('<p>', {text: 'There are no In Progress Flags or important Resolved Flags attached to this record.'}).appendTo(chunk);
                }

                const raisableFlags = this.collection.raisableFlags();
                const notSeenRaisableFlags = [];
                for (const raisableFlag of raisableFlags) {
                    if (!seenOnlyOneFlagTypes[raisableFlag.flag_type]) {
                        notSeenRaisableFlags.push(raisableFlag);
                    }
                }

                if (notSeenRaisableFlags.length) {
                    this.flagChunk(notSeenRaisableFlags, 'Raise New Flags', 'new', owner).appendTo(dom);
                }

                // this is mostly redundant to Variant Classification History
                const timelineA = $('<a>', {text: 'See Full Timeline', class:'btn btn-outline-secondary'});
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
    
    const FlagCollectionDialog = (function() {
        const FlagCollectionDialog = function(collection) {
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
                const activeFlag = params.activeFlag;
                
                const modalContent = $(`
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
                setupModalAnimationForWebTesting(modalContent);
                modalContent.find('.modal-header').html(`
                    <nav aria-label="breadcrumb">
                      <ol class="breadcrumb"></ol>
                    </nav>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                        <span aria-hidden="true">&times;</span>
                    </button>
                `);

                this.content = modalContent;

                const db = this.collection.db;
                const dialog = this;
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
                const parentTitle = (this.collection.label || 'Flags').replace('/', '-'); // fixme trim title
                const title = this.activeContent.title;
                let titleDom = this.activeContent.titleDom;
                if (!titleDom && title) {
                    titleDom = $('<span/>', {text: title});
                }
                const breadcrumbs = this.content.find('ol.breadcrumb');
                if (!title && !titleDom) {
                    breadcrumbs.html([
                        $('<li>', {class: 'breadcrumb-item', text: parentTitle})
                    ]);
                } else {
                    // not showing parent title as it can get a bit out of control
                    breadcrumbs.html([
                        $('<li>', {class: 'breadcrumb-item font-weight-bold', html: $('<a>', {html: '<i class="fas fa-angle-left"></i> Back to All Flags', click: () => { this.back(); }})}),
                        $('<li>', {class: 'breadcrumb-item active', html: titleDom})
                    ]);
                }
                this.content.find('.modal-body').html(this.activeContent.body);
                const footer = this.content.find('.modal-footer');
                footer.html(this.activeContent.footer);
                if (!this.activeContent.footer || this.activeContent.footer.css('display') === 'none') {
                    footer.addClass('d-none');
                } else {
                    footer.removeClass('d-none');
                }

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

    const Flag = (function () {
        const Flag = function (db, data) {
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
                const userPermission = this.collectionObj().user_permission;
                const requiredPermission = this.flagTypeObj().permission;
                return userPermission >= requiredPermission;
            },
            subFlags() {
                const extraData = this.extra_data;
                if (extraData) {
                    const subFlags = extraData.sub_flags;
                    if (subFlags) {
                        return subFlags;
                    } else {
                        return [];
                    }
                }
            },
            
            statusText() {
                if (this.flag_type === 'classification_significance_change') {
                    return 'Primary Reason for Change';
                }
                if (this.creating) {
                    return 'Opening Status';
                }
                return 'Status';
            },
            
            resolutions() {
                const resolutions = this.flagTypeObj().resolutions.map(r => this.db.flagResolutions.get(r) );
                
                if (this.creating) {
                    const excludeStatuses = {'R': true};
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
                        If you have found some extra information that you think should be incorporated into this classification record,
                        you can raise a suggestion for the classification record owner to accept or reject.
                        </div>`);
                    } else if (this.flag_type === 'classification_internal_review') {
                        return $(`<div>
                        You can raise this flag to let people know the classification record is currently in review, or raise it
                        as "Completed" to record the fact that a review has recently taken place.<br/>
                        Please record any internal reviews while a classification is marked as discordant.
                        </div>`);
                    } else if (this.flag_type === 'classification_pending_changes') {
                        return $(`<div>
                        Mark this classification as having important changes not yet reflected in this system.
                        </div>`);
                    } else if (this.flag_type === 'classification_not_public') {
                        return $(`<div>
                        Raise this flag to stop this specific record from being sent to ClinVar.<br/>
                        ClinVar exclusion patterns can be setup by administrators if required.
                        </div>`);
                    }
                    return $('<div>');
                }
            
                if (this.flag_type === 'classification_suggestion') {
                    return $(`<div>
                    Someone has raised suggestion(s) against this classification record.
                    <ol><li>Review the contents of each suggestion.</li>
                    <li>If appropriate, make changes in your curation system and mark the suggestion as Complete.</li>
                    <li>If you decline the suggestion, mark it as Rejected.</li>
                    </ol></div>
                    `);
                } else if (this.flag_type === 'classification_outstanding_edits') {
                    return $(`<div>
                    Edits have been made to this classification record that are not included in a published version.
                    <ol><li>From the classification record form, ensure there are no validation errors stopping this record from being published.</li>
                    <li>At the bottom of the form, click the tick to submit the outstanding changes.</li></ol></div>`
                    );
                } else if (this.flag_type === 'classification_internal_review') {
                    return $(`<div>
                    This classification record is marked as currently being internally reviewed.
                    <ol><li>Once the internal review is complete, ensure you update the classification in your curation system.</li>
                    <li>Mark the internal review as Completed</li></ol>
                    </div>
                    `);
                } else if (this.flag_type === 'classification_withdrawn') {
                    return $(`<div>
                    This classification record has been marked as withdrawn. It will be hidden from almost all searches and exports.
                    <ol><li>If the classification is not of high enough quality or in error, you may leave it as "withdrawn" indefinitely.</li>
                    <li>If you wish to un-withdraw the classification record, click the open bin icon in actions from the variant classification record form</li></ol></div>
                    `);
                } else if (this.flag_type === 'classification_significance_change') {
                    return $(`<div>
                    This classification record has changed its classification or clinical significance compared to a previously published version.
                    <ol><li>Set the status of this flag to reflect the primary reason behind the change in assertion</li>
                    <div><ul>
                    <li>Discordance Discussion - Data was changed as a result of talking to other labs when this classification was in discordance.
                    <li>Summation of Data - Data was changed as a result of combining information from multiple labs.
                    <li>Internal Review - New data or errors in old data were found during an internal review.
                    </ul>                 
                    </div>
                    <li>Please also add a comment providing some context.</li></ol></div>`);
                } else if (this.flag_type === 'classification_discordant') {
                    return $(`<div>
                    This classification record is in discordance with one or more classification records.
                    <ol><li>Ensure that you have completed an internal review of your lab's classification recently (within the last 12 months is recommended). If not, raise the internal review flag and complete an internal review of your lab's classification.
                    <li>Review any outstanding suggestions against your lab's classification.
                    <li>View the other classifications in the discordance report and view the evidence differing between multiple records via the diff page. If appropriate, raise suggestions against other lab classifications.
                    <li>This Discordance flag will automatically be closed when concordance is reached.
                    </ol></div>
                    `);
                } else if (this.flag_type === 'classification_pending_changes') {
                    return $(`<div>
                    This classification record has outstanding changes to apply.<br/>
                    If a subsequent sync to this system provides a new clinical significance this flag will automatically close.<br/>
                    Otherwise it can be closed manually if the flag was raised in error or circumstances have changed.
                    </div>`);
                } else if (this.flag_type === 'classification_unshared') {
                    return $(`<div>
                    This classification record is not yet shared outside of your lab or institution.
                    <ol><li>From the classification record form, ensure there are no validation errors stopping this record from being published.</li>
                    <li>Review the content of the classification record to make sure it's ready to be shared.</li>
                    <li>At the side of the form, click the Share to submit at a higher share level.</li></ol></div>
                    `);
                }
                return $('<div>');
            },
            
            statusDom() {
                if (this.creating) {
                    return $('<div>');
                }
                return $('<div>', {html: [
                    `This ${this.flagTypeObj().label} has the status of `,
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
                    const transcript = this.data.transcript;
                    parts = /^([_A-Z0-9]+)(?:[.]([0-9]+))?$/i;
                    const match = parts.exec(transcript);
                    if (match) {
                        let transcriptUrl = null;
                        const transcriptNoVer = match[1];
                        transcriptUrl = Urls.view_transcript(transcriptNoVer);
                        return $('<div>', {html: [
                            'View details about the transcript ',
                            $('<a>', {href:transcriptUrl, class:'hover-link', text:`${transcriptNoVer}`, target:'_blank'}),
                            '.'
                        ]});
                    }
                } else if (this.flag_type === 'classification_suggestion') {
                    const user = this.userObj();
                    return $('<div>', {html: [
                        `Raised by `,
                        $('<span>', {class: 'username', text:user.name}),
                        $('<span>', { class: 'text-secondary text-italic font-italic d-inline-block mx-2', text: user.lab}),
                        ` in regards to `,
                        $('<span>', {class: 'quote', text: firstComment})
                    ]});
                } else if (this.flag_type === 'classification_significance_change') {
                    return $('<div>', {text: firstComment});
                } else if (this.flag_type === 'classification_matching_variant_warning') {
                    const variantId = this.collectionObj().variant;
                    return $('<div>', {html: [
                        `See more information about the linked variant `,
                        $('<a>', {class: 'hover-link', text: `here`, href: Urls.view_allele_from_variant(variantId)})
                    ]});
                } else if (this.flag_type === 'classification_discordant') {
                    const reportId = this.collectionObj().discordance_report;
                    const variantId = this.collectionObj().variant;
                    const clinicalContext = this.collectionObj().clinical_context;
                    if (reportId) {
                        return $('<div>', {html: [
                            `Go to the `,
                            $('<a>', {class: 'hover-link', text: `Discordance Report`, href:`/classification/classification/discordance_report/${reportId}`}),
                            ` | `,
                            $('<a>', {class: 'hover-link', text: `Diff with other Classification Records`, href: `/classification/diff/?clinical_context=${clinicalContext}`})
                        ]});
                    }
                }
                return $('<div>');
            },
            
            ageText() {
                return jQuery.timeago(this.created * 1000);
            },
            comments() {
                return Object.values(this.db.comments.all()).filter(c => c.flag == this.id).sort((c1, c2) => c1.created - c2.created);
            },
            dom() {
                const flagType = this.flagTypeObj();
                const user = this.userObj();

                const titleText = `${flagType.label} (${this.ageText()})`;

                const collectionObj = this.collectionObj();
                const flagDiv = $('<div>', { class: `flag flag-${flagType.id} res-${this.resolution}`, title: titleText});
                flagDiv.click(() => {
                    flagDiv.tooltip('hide');
                    new FlagCollectionDialog(this.collectionObj()).init({activeFlag:this, triggerDom:flagDiv})
                });
                if (this.open === false) {
                    flagDiv.addClass('closed');
                }
                return flagDiv;
            },
            remove() {
                for (const comment of this.comments()) {
                    comment.remove();
                }
                this.db.flags._remove(this.id);
            },
        };
        return Flag;
    })();


    const FlagComment = (function () {
        const FlagComment = function (db) {
            this.db = db;
        };
        FlagComment.prototype = {
            flagObj() { return this.db.flags.get(this.flag); },
            userObj() { return this.db.users.get(this.user); },
            resolutionObj() { return this.resolution ? this.db.flagResolutions.get(this.resolution) : null; },
            action() {
                const resolution = this.resolutionObj();
                if (!resolution) {
                    return 'Commented';
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
                const { includeTitle, lastEntry, dialog } = params || {};

                const flag = this.flagObj();
                const user = this.userObj();
                const flagType = flag.flagTypeObj();
                const title = `${flagType.label} (${flag.ageText()})`;
                const commentDom = $('<div>', { class: 'flag-comment' });

                let timestamp = moment(this.created * 1000).format('DD-MMM-YYYY');
                if (lastEntry) {
                    oldTimestamp = moment(lastEntry.created * 1000).format('DD-MMM-YYYY');
                    if (oldTimestamp === timestamp) {
                        timestamp = null;
                    }
                }
                if (timestamp) {
                    $('<div>', { class: 'timestamp flag-title', text: timestamp }).appendTo(commentDom);
                    commentDom.addClass('with-title');
                }

                const userActionDom = $('<div>', { class: 'user-action' }).appendTo(commentDom);

                const avatar = $('<div>', { class: 'avatar' }).css(`background-image`, `url(${user.avatar})`).css('background-color', user.color).appendTo(userActionDom);
                const content = $('<div>', { class: 'content' }).appendTo(userActionDom);
                const time = moment(this.created * 1000).format('hh:mm A');
                const userInfo = $('<div>', {
                    class: 'header', html: [
                        $('<span>', { class: 'username', text: user.name }),
                        $('<span>', { class: 'text-secondary text-italic font-italic d-inline-block mr-2', text: user.lab}),
                        $('<span>', { class: 'time', text: time })
                    ]
                }).appendTo(content);

                if (!timestamp && lastEntry && lastEntry.user == this.user && time == moment(lastEntry.created * 1000).format('hh:mm A')) {
                    avatar.css('height', '0').css('visibility', 'hidden');
                    userInfo.css('display', 'none');
                }

                const actionCommentDom = $('<div>', { class: 'action-text' }).appendTo(content);
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


    const Flags = function (
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
                const counts = {};
                for (const flag of this.flags.all()) {
                    if (flag.open) {
                        const resolution_counts = counts[flag.flag_type] || {};
                        resolution_counts[flag.resolution] = (resolution_counts[flag.resolution] || 0) + 1;
                        counts[flag.flag_type] = resolution_counts;
                    }
                }
                const ordered_keys = Object.keys(counts).sort((k1, k2) => this.flagTypes.get(k1).label.localeCompare(this.flagTypes.get(k2).label));
                for (const key of ordered_keys) {
                    
                    if (this.summaryFilterInclusions && !this.summaryFilterInclusions[key]) {
                        continue;
                    }
                    
                    const resolution_counts = counts[key];
                    const flagType = this.flagTypes.get(key);
                    for (const [resolution, value] of Object.entries(resolution_counts)) {

                        const label =  $('<span>', {class: 'label', text: flagType.label});
                        if (resolution !== 'open') {
                            const resolutionObj = this.flagResolutions.get(resolution);
                            $('<span>', {text: ` (${resolutionObj.label})`}).appendTo(label);
                        }
                        
                        const cell_parts = [
                            $('<div>').addClass('flag').addClass(`flag-${flagType.id} res-${resolution}`),
                            $('<div>', {class: 'text-monospace mx-2', text: value.toString(), style:'min-width: 2rem; text-align:right'}),
                            label
                        ];
                        if (value > 1 && !flagType.label.endsWith('s')) {
                            $('<span>', {class: 'plural', text: 's'}).appendTo(label);
                        }
                        
                        const filterCell = $('<a/>', {
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
            const filterValue = `${flag_type_id}-${resolution}`;
            this.summaryFilter.find('.list-group-item').removeClass('list-group-item-primary');
            if (this.filterValue === filterValue || flag_type_id === null) {
                this.filterValue = null;
                this.collections.all().forEach(c => c.dom.closest('tr').show());
            } else {
                this.filterValue = filterValue;
                this.summaryFilter.find(`[data-flag-type=${flag_type_id}-${resolution}]`).addClass('list-group-item-primary');
                const matching_collection_ids = {};
                const flags = this.flags.all().filter(f => f.flag_type === flag_type_id && f.resolution === resolution).forEach(f => {
                    matching_collection_ids[f.collection] = true;
                });
                for (const fc of this.collections.all()) {
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
            let flagDoms = null;
            const flagGroup = props.flagGroup || 'default';
            if (props.filter) {
                flagDoms = $(`${props.filter} *[data-flags]`);
            } else {
                flagDoms = $('[data-flags]');
            }
            
            this.onClose = props.onClose || false;
            this.summaryFilter = $('#flags-filter');
            if (this.summaryFilter) {
                const summaryFilterInclusions = this.summaryFilter.attr('data-filter') || null;
                if (summaryFilterInclusions) {
                    this.summaryFilterInclusions = {};
                    for (const inclusion of summaryFilterInclusions.split(' ')) {
                        this.summaryFilterInclusions[inclusion.trim()] = true;
                    }
                }
            }

            let changes = false;
            const existingCollections = {};
            for (const collection of this.collections.all()) {
                existingCollections[collection.id] = collection;
            }
            const flag_collection_doms = {};
            flagDoms.toArray().forEach(e => {
                const dom = $(e);
                dom.addClass('flags');
                const id = parseInt(dom.attr('data-flags'));

                if (!isNaN(id)) {
                    flag_collection_doms[id] = Object.assign({ dom: dom }, flag_collection_doms[id] || {});
                }
            });
            for (const entry of Object.entries(flag_collection_doms)) {
                const id = entry[0];
                const data = entry[1];
                if (existingCollections[id]) {
                    delete existingCollections[id];
                } else {
                    changes = true;
                }

                this.collections.upsert({
                    id: id,
                    dom: data.dom,
                    flagGroup: flagGroup
                });
            }

            for (const collection of Object.values(existingCollections)) {
                if (collection.flagGroup === flagGroup) {
                    collection.remove();
                }
            }
            if (changes || props.forceRender) {
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
                const collections = this.collections.all();
                if (collections.length === 0) {
                    resolve();
                    return;
                }
                let collection = params.collection;
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
            for (const user of (record.users || [])) {
                this.users.upsert(user);
            }
            for (const flag of (record.flags || [])) {
                this.flags.upsert(flag);
            }
            for (const comment of (record.comments || [])) {
                this.comments.upsert(comment);
            }
            for (const flagResolution of (record.flag_resolutions || [])) {
                this.flagResolutions.upsert(flagResolution);
            }
            for (const flagType of (record.flag_types || [])) {
                this.flagTypes.upsert(flagType);
            }
            for (const collection of (record.collections || [])) {
                this.collections.upsert(collection);
            }
        },

        render() {
            for (const collection of this.collections.all()) {
                const dom = collection.dom;
                if (!dom) {
                    continue;
                }
                dom.empty();
                const flags = collection.flags();
                const watch = collection.watching === 0 || collection.watching;
                const flagSummary = $('<div>', { class: `flag add`, title: 'Add or review flags for ' + collection.label}).appendTo(dom);
                flagSummary.click(() => { new FlagCollectionDialog(collection).init({triggerDom: flagSummary}) });

                if (collection.watching) {
                    $('<div>', { class: 'notifications', title: `${collection.watching} unseen activities` }).appendTo(flagSummary);
                } else if (watch) {
                    $('<div>', { class: `notifications watching`, title: `starred` }).appendTo(flagSummary);
                }
                
                for (const flag of collection.flags()) {
                    const flagDom = flag.dom().appendTo(dom);
                    dom.append(flagDom);
                    const subFlags = flag.subFlags();
                    if (subFlags && subFlags.length) {
                        const subFlagsDom = $('<div>', {class: 'sub-flags'}).appendTo(dom);
                        for (const subFlag of subFlags) {
                            const subFlagDom = $('<div>', {class: `sub-flag`, title: subFlag.label, text: subFlag.letter || 'X'});
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