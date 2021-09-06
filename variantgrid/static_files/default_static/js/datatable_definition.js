let DataTableDefinition = (function() {
    "use strict";

    let DataTableDefinition = function(params) {
        this.dom = params.dom;
        this.url = params.url;
        this.data = params.data;
        this.filterCount = params.filterCount;
        this.waitOn = Promise.resolve();

        this.tableId = null;
        this.lengthKey = null;
        this.serverParams = null;
        this.dtParams = null;
        this.dataTable = null;
        this.expandData = {
            expandedTr: null,
            expandedRow: null,
            hoverTimeout: null
        };
    };

    DataTableDefinition.definitions = {};

    DataTableDefinition.prototype = {

        ensureState: function() {
            if (!this.dom) {
                throw "DatatTableDefinition missing parameter `dom`";
            }
            if (!this.url) {
                throw "DatatTableDefinition missing parameter `url`";
            }

            let tableId = this.dom.attr('id');
            if (!tableId) {
                tableId = "tid" + _.random(0,100000);
                this.dom.attr('id', tableId);
            }
            // this.dom.style('width', '100%');
            this.tableId = tableId;
            this.lengthKey = `datatable_length_${tableId}`;

            let dom = this.dom;
            dom.addClass('table');
            dom.addClass('stripe');
            dom.addClass('dataTable');
        },

        loadDefinition: function() {
            let definitionData = DataTableDefinition.definitions[this.url];
            if (!definitionData) {
                definitionData = $.getJSON(this.url + '?dataTableDefinition=1');
                DataTableDefinition.definitions[this.url] = definitionData;
            }
            return definitionData.then(data => {this.serverParams = data});
        },

        convertDefinition: function() {
            let defn = this.serverParams;
            let tableId = this.tableId;
            let lengthKey = this.lengthKey;

            let lengthValue = 10;
            if (tableId) {
                lengthValue = parseInt(localStorage.getItem(lengthKey)) || 10;
            }

            let domString = `<"top"><"toolbar"<"custom">${ defn.searchBoxEnabled ? 'f' : ''}>rt<"bottom"ilp><"clear">`;

            let dtParams = {
                processing: true,
                serverSide: true,
                pageLength: lengthValue,
                dom: domString,
                order: defn.order,
                ajax: {
                    url: this.url,
                    type: 'POST',
                    data: this.data ? eval(this.data) : null
                },
                bFilter: defn.searchBoxEnabled,
                bAutoWidth: false
            };
            if (this.filterCount === 'hide') {
                dtParams.langauge = {infoFiltered: ""};
            }
            if (defn.responsive) {
                // might need to have a calculated with for this to hide columns automatically?
                dtParams.responsive = {
                    details: {
                        type: 'column',
                        target: 'td.dt-preview',
                        renderer: TableFormat.detailRendererHtml
                    }
                };
            }

            let columnDefs = [];
            dtParams.columnDefs = columnDefs;

            let waitOnEKeys = null;

            // GENERATE COLUMNS
            for (let col of defn.columns) {
                let columnDef = Object.assign({}, col);
                let target = columnDefs.length;
                columnDef.targets = target;
                columnDefs.push(columnDef);
                if (col.render) {
                    columnDef.render = eval(col.render);
                    if (col.render.includes('VCTable')) {
                        waitOnEKeys = true;
                    }
                }
                if (col.createdCell) {
                    columnDef.createdCell = eval(col.createdCell);
                }
                if (target === 0 && defn.expandClientRenderer) {
                    // if we're expanding rows with ajax, make first column have the toggle
                    columnDef.className = (columnDef.className || "") + " toggle-link";
                }
            }

            if (waitOnEKeys) {
                this.waitOn = EKeys.load();
            }

            this.dtParams = dtParams;
            return dtParams;
        },

        setupDom: function() {
            let dom = this.dom;
            let dtParams = this.dtParams;

            let tHead = $('<thead/>').appendTo(dom);
            let tHeadTr = $('<tr/>').appendTo(tHead);

            // GENERATE COlumns
            for (let columnDef of dtParams.columnDefs) {
                $('<th/>', {class: columnDef.classNames, html: columnDef.label}).appendTo(tHeadTr);
            }

            let dataTable = dom.DataTable(dtParams);
            this.dataTable = dataTable;

            dom.on('error.dt', function ( e, settings, techNote, message ) {
                Rollbar.warning("DataTables error " + message);
                    console.log( 'An error has been reported by DataTables: ', message );
                }
            );

            let tableId = this.tableId;
            $(`select[name=${tableId}_length]`).change(function() {
                localStorage.setItem(lengthKey, $(this).val());
            });
            // move any externally defined toolbar elements onto it
            $(`[data-toolbar="#${tableId}"]`).detach().appendTo(dom.closest('.dataTables_wrapper').find('.toolbar .custom'));
        },

        setupClientExpend: function() {
            if (!this.serverParams.expandClientRenderer) {
                return;
            }

            let dom = this.dom;
            dom.addClass('expandable');

            let dataTable = this.dataTable;
            let expandFn = eval(this.serverParams.expandClientRenderer);
            let expandData = this.expandData;

            dom.on('click', 'tr', function() {
                let tr = $(this); //.closest('tr');
                if (!tr.hasClass('odd') && !tr.hasClass('even')) {
                    // not a regular row
                    return;
                }
                let row = dataTable.row( tr );
                if ( row.child.isShown() ) {
                    // This row is already open - close it
                    row.child.hide();
                    tr.removeClass('shown');

                    expandData.expandedRow = null;
                    expandData.expandedTr = null;
                } else {
                    // Close previous row (if there is one)
                    try {
                        if (expandData.expandedRow && expandData.expandedTr) {
                            expandData.expandedRow.child.hide();
                            expandData.expandedTr.removeClass('shown');
                        }
                    } catch (ex) {
                        // no-op
                    }
                    expandData.expandedRow = row;
                    expandData.expandedTr = tr;

                    if (!tr.hasClass('loaded')) {
                        // loading hasn't started yet, load the row
                        let childHtml = expandFn(row.data());
                        row.child( childHtml );
                        tr.addClass('loaded');
                    }
                    // show the row
                    row.child.show();
                    tr.addClass('shown');
                }
            });
            // PRE-FETCH data
            // if hovering over a single row for 500ms, pre-fetch the client data ready to display
            dom.on('mouseenter', 'tr', function() {
                let tr = $(this);
                if (!tr.hasClass('odd') && !tr.hasClass('even') || tr.hasClass('loaded')) {
                    return; // either not an odd or even row, or already
                }
                window.clearTimeout(expandData.hoverTimeout);
                expandData.hoverTimeout = window.setTimeout(() => {
                    if (!tr.hasClass('loaded')) { // could have been clicked on
                        let row = dataTable.row( tr );
                        let childHtml = expandFn(row.data());
                        row.child( childHtml );
                        tr.addClass('loaded');
                        tr.addClass('pre-fetched'); // in case we ever want to do stats on it
                    }
                }, 500);
            });
            dom.on('mouseleave', 'tr', function() {
                window.clearTimeout(expandData.hoverTimeout);
            });
        },

        setup: function() {
            if (this.dom.hasClass('dataTable')) {
                return;
            }
            this.ensureState();
            this.loadDefinition().then(() => {
                this.convertDefinition();
                this.waitOn.then(() => {
                    this.setupDom();
                    this.setupClientExpend();
                    this.dataTable.columns.adjust().draw();
                });
            });
        }
    };

    return DataTableDefinition;
})();


// ******************************************************************************************
// COMMON TABLE FORMATTERS
// ******************************************************************************************


let TableFormat = (function() {
    "use strict";
    let TableFormat = function() {};
    TableFormat.prototype = {};
    return TableFormat;
})();

TableFormat.timestamp = (data, type, row) => {
    if (data) {
        let timestampStr = convertTimestamp(data);
        return $('<span>', {class:'timestamp', text: timestampStr}).prop('outerHTML');
    } else {
        return '';
    }
};

TableFormat.timeAgo = (data, type, row) => {
    if (data) {
        return $('<data>', {class:'convert-timestamp time-ago', 'data-timestamp':data, text:data}).prop('outerHTML');

        let timestampStr = convertTimestamp(data);
        return $('<span>', {class:'timestamp', text: timestampStr}).prop('outerHTML');
    } else {
        return '';
    }
};

TableFormat.choices = (choices, data, type, row) => {
    return $('<span>', {class:`val-${data}`, text:choices[data] || data}).prop('outerHTML');
};

TableFormat.flags = (data, type, row) => {
    if (data) {
        return $('<div>', {'data-flags': data, class:'flags', text:'...'}).prop('outerHTML');
    }
};

TableFormat.limit = (limit, data, type, row) => {
    if (typeof(data) === 'string' && data.length > limit) {
        return $('<span>', {class:'hover-detail', text:data.substring(0, limit) + '...', title:data}).prop('outerHTML');
    }
    return data;
};

TableFormat.text = (data, type, row) => {
    if (data === '' || data === null) {
        return $('<span/>', {class:'no-value', text:'-'}).prop('outerHTML');
    } else {
        return data;
    }
};

TableFormat.preview = (columns, data, type, row) => {
    let dom = $('<div>');
    let hasValue = false;
    for (let col of columns) {
        let value = row[col];
        if (value && value.length) {
            hasValue = true;
            if (value.length > 80) {
                value = value.substring(0, 80) + '...';
            }
            $('<div>', {text: value}).appendTo(dom);
        }
    }
    if (!hasValue) {
        dom = $('<div>', {class:'no-value', text:'-'});
    }
    return dom.prop('outerHTML');
};

TableFormat.boolean = function(style, data, type, columns) {
    if (style === 'warning') {
        if (data) {
            return '<i class="fas fa-exclamation-circle"></i>';
        }
    } else {
        return data ? '<i class="fas fa-check-circle text-success"></i>' : '<i class="far fa-circle"></i>';
    }
    return null;
};

TableFormat.severeNumber = function(severity, data, type, columns) {
    if (data === 0) {
        return '<span class="no-value mono font-weight-bold">0</span>';
    } else {
        return `<span class="mono font-weight-bold text-${severity}">${data}</span>`;
    }
};

TableFormat.hgvs = function(data, type, columns) {
    if (!data) {
        return "?";
    }
    let genomeBuild = data.genomeBuild;
    let transcript = data.transcript;
    let geneSymbol = data.geneSymbol;
    let variant = data.variant;
    let allele = data.allele;
    // also turn into a link
    let dom = $('<div>');
    if (allele) {
        dom.append($('<div>', {text: allele, class:'font-weight-bold'}));
    }

    if (genomeBuild) {
        dom.append($('<div>', {text: genomeBuild, class:'text-info'}));
        // <span style="white-space: nowrap"><span>{{ c_hgvs.transcript }}</span>{% if c_hgvs.gene_symbol %}(<span class="text-secondary" style="letter-spacing: 0.5px">{{ c_hgvs.gene_symbol }}</span>){% endif %}:</span><span style="display:inline-block;word-break: break-all">{{ c_hgvs.raw_c }}</span>
    }

    let cDom = $('<span>', {style:'white-space:nowrap'});
    if (transcript && variant) {
        cDom.append($('<span>', {text: transcript}));
        if (geneSymbol) {
            cDom.append("(");
            cDom.append($('<span>', {class: 'text-secondary', style: 'letter-spacing: 0.5px', text: geneSymbol}));
            cDom.append(")");
        }
        cDom.append(":");
        // used to be display:inline-block; but that doesn't underline
        cDom.append($('<span>', {style: 'word-break:break-all', text: variant}));
    } else {
        cDom.append(data.full);
    }
    let variantId = data.variantId;
    if (variantId) {
        cDom = $('<a>', {href: Urls.view_allele_from_variant(variantId), class:'hover-link', html: cDom});
    }
    dom.append(cDom);

    return dom.prop('outerHTML');
};

TableFormat.expandAjax = function(url, param, expectedHeight, data) {
    if (data) {
        let dataId = data[param];
        if (!dataId) {
            return `<i class="fas fa-bomb text-danger"></i> No value for "${param}" in this ${JSON.stringify(data)} : DEBUG - is ${param} a column in this table, visible or otherwise?`;
        }
        let ajaxId = `ajax_${dataId}`;
        let reverseUrl = Urls[url];
        if (!reverseUrl) {
            return `<i class="fas fa-bomb text-danger"></i> URL not configured for "${url} : Developer may need to run<br/>
            <div class="code">manage.py collectstatic_js_reverse</div>`;
        }
        if (param) {
            reverseUrl = reverseUrl(dataId);
        }

        let ajaxDom =
            $('<div>', {html:[
                $('<div>', {style:`text-align: center;color: #888; min-height:${expectedHeight}`, text:'Loading...'})
            ]});

        loadAjaxBlock(ajaxDom, reverseUrl);

        if (!expectedHeight) {
            // put a small div in the row to show that we're thinking
            expectedHeight = "50px";
        }
        return ajaxDom;
    } else {
        return '';
    }
};

TableFormat.detailRendererHtml = function ( api, rowIdx, columns ) {
    let fieldset = $('<div>', {class:'mt-3'});
    for (let col of columns) {
        if (col.hidden) {
            if (col === null || col.data.length === 0) {
                // pass
            } else {
                $('<div>', {
                    class: 'row mt-2', style:'align-items:center', html: [
                        $('<div>', {
                            class: 'col-2 text-right', html:
                                $('<label>', {html: col.title})
                        }),
                        $('<div>', {
                            class: 'col-10', html:
                                $('<span>', {class: 'dt-detail', html: col.data})
                        }),
                    ]
                }).appendTo(fieldset);
            }
        }
    }
    return fieldset;
};