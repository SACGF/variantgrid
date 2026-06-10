const oldExportAction = function (self, e, dt, button, config) {
    if (button[0].className.indexOf('buttons-csv') >= 0) {
        if ($.fn.dataTable.ext.buttons.csvHtml5.available(dt, config)) {
            $.fn.dataTable.ext.buttons.csvHtml5.action.call(self, e, dt, button, config);
        }
    }
};

const newExportAction = function (e, dt, button, config) {
    /* Pull all data from Ajax - but don't do the draw
       Code taken from https://stackoverflow.com/a/44635032/295724 */
    const self = this;
    const oldStart = dt.settings()[0]._iDisplayStart;

    dt.one('preXhr', function (e, s, data) {
        // Just this once, load all data from the server...
        data.start = 0;
        data.length = 2147483647;

        dt.one('preDraw', function (e, settings) {
            // Call the original action function
            oldExportAction(self, e, dt, button, config);

            dt.one('preXhr', function (e, s, data) {
                // DataTables thinks the first item displayed is index 0, but we're not drawing that.
                // Set the property to what it was before exporting.
                settings._iDisplayStart = oldStart;
                data.start = oldStart;
            });

            // Reload the grid with the original page. Otherwise, API functions like table.cell(this) don't work properly.
            setTimeout(dt.ajax.reload, 0);

            // Prevent rendering of the full data to the DOM
            return false;
        });
    });

    // Requery the server with the new one-time export settings
    dt.ajax.reload();
};


const DataTableDefinition = (function() {
    "use strict";

    const DataTableDefinition = function(params) {
        this.dom = params.dom;
        this.url = params.url;
        this.data = params.data;
        this.filterCount = params.filterCount;
        this.waitOn = Promise.resolve();
        this.adjustColumns = params.adjustColumns !== false;

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
            // this seems to happen after convertDefinition, so it's relatively useless
            if (!this.dom) {
                throw new Error("DatatTableDefinition missing parameter `dom`");
            }
            if (!this.url) {
                throw new Error("DatatTableDefinition missing parameter `url`");
            }

            let tableId = this.dom.attr('id');
            if (!tableId) {
                tableId = "tid" + _.random(0,100000);
                this.dom.attr('id', tableId);
            }
            // this.dom.style('width', '100%');
            this.tableId = tableId;
            this.lengthKey = `datatable_length_${tableId}`;

            const dom = this.dom;
            dom.addClass('table');
            dom.addClass('stripe');
            dom.addClass('dataTable');
        },

        loadDefinition: function() {
            let definitionData = DataTableDefinition.definitions[this.url];
            if (!definitionData) {
                let sep = '?';
                if (this.url.indexOf('?') !== -1) {
                    sep = '&';
                }
                definitionData = $.getJSON(this.url + sep + 'dataTableDefinition=1');
                DataTableDefinition.definitions[this.url] = definitionData;
            }
            return definitionData.then(data => {this.serverParams = data});
        },

        convertDefinition: function() {
            const defn = this.serverParams;
            const tableId = this.tableId;
            const lengthKey = this.lengthKey;

            let lengthValue = 10;
            if (tableId) {
                lengthValue = parseInt(localStorage.getItem(lengthKey)) || 10;
            }

            const domString = `<"top"><"toolbar"<"custom">${ defn.searchBoxEnabled ? 'f' : ''}>rt${ defn.downloadCsvButtonEnabled ? 'B' : ''}<"bottom"<"showing"il>p><"clear">`;

            const dtParams = {
                processing: true,
                serverSide: true,
                pageLength: lengthValue,
                dom: domString,
                order: defn.order,
                fixedOrder: defn.order,
                pagingType: "input",
                classes: {
                    'sPageButton': 'btn btn-outline-primary btn-rnd-rect',
                    'sPageButtonDisabled': 'disabled'
                },
                ajax: {
                    url: this.url,
                    type: 'POST',
                    data: this.data ? eval(this.data) : null
                },
                bFilter: defn.searchBoxEnabled,
                bAutoWidth: false,
                scrollX: defn.scrollX,
                initComplete:  function( settings, json ) {
                    $('th.toggle-link').removeClass('toggle-link');
                }
            };
            if (defn.order) {
                dtParams.orderSequence = defn.orderSequence
            }

            if (defn.downloadCsvButtonEnabled) {
                const csvName = defn.csvName || 'export';
                const dateStr = new Date().toISOString().slice(0, 10);
                dtParams.buttons = [
                    {
                        extend: 'csvHtml5',
                        action: newExportAction,
                        text: "Download as CSV",
                        filename: csvName + '_' + dateStr,
                    }
                ]
            }

            if (this.filterCount === 'hide') {
                dtParams.langauge = {infoFiltered: ""};
            }
            if (defn.responsive) {
                // might need to have a calculated with for this to hide columns automatically?
                dtParams.responsive = {
                    details: {
                        type: 'column',
                        target: 'tr',
                        renderer: TableFormat.detailRendererHtml
                    }
                };
            }

            const columnDefs = [];
            dtParams.columnDefs = columnDefs;

            let waitOnEKeys = null;

            if (!Array.isArray(defn.columns)) {
                console.log("Invalid TableDefinition")
                console.log(defn);
                throw new Error("Received invalid datatable definition");
            }

            // GENERATE COLUMNS
            for (const col of defn.columns) {
                const columnDef = Object.assign({}, col);
                const target = columnDefs.length;
                columnDef.targets = target;
                columnDefs.push(columnDef);
                if (col.render) {
                    const rawRenderer = eval(col.render);
                    const renderer = (data, type, row) => {
                        const output = rawRenderer(data, type, row);
                        if (output instanceof jQuery) {
                            return output.prop("outerHTML");
                        }
                        return output;
                    }
                    columnDef.render = renderer;
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
            const dom = this.dom;
            dom.empty();
            const dtParams = this.dtParams;

            const tHead = $('<thead/>').appendTo(dom);
            const tHeadTr = $('<tr/>').appendTo(tHead);

            // GENERATE COLUMNS
            for (const columnDef of dtParams.columnDefs) {
                $('<th/>', {class: columnDef.classNames, html: columnDef.label}).appendTo(tHeadTr);
            }

            dtParams.createdRow = (row, data, dataIndex) => {
                const row_css = data.row_css;
                if (row_css) {
                    $(row).addClass(row_css);
                }
            }

            const dataTable = dom.DataTable(dtParams);
            this.dataTable = dataTable;

            dom.on('error.dt', function (e, settings, techNote, message ) {
                Rollbar.warning("DataTables error " + message);
                dom.replaceWith($(`
                    <div class="w-100 m-4 p-4 rounded border border-danger text-center">
                        <p class="text-danger"><i class="fas fa-bomb"></i> There was an error loading this table.<br/><br/><a class="font-weight-bold" href="javascript:window.location.reload(true)">Click Here to reload the page.</a></p>
                    </div>`)
                );
            });

            const tableId = this.tableId;
            const lengthKey = this.lengthKey;

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

            const dom = this.dom;
            dom.addClass('expandable');

            const dataTable = this.dataTable;
            const expandFn = eval(this.serverParams.expandClientRenderer);
            const expandData = this.expandData;

            dom.on('click', 'tr', function() {
                const tr = $(this); //.closest('tr');
                if (!tr.hasClass('odd') && !tr.hasClass('even')) {
                    // not a regular row
                    return;
                }
                const row = dataTable.row( tr );
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
                        const childHtml = expandFn(row.data());
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
                const tr = $(this);
                if (!tr.hasClass('odd') && !tr.hasClass('even') || tr.hasClass('loaded')) {
                    return; // either not an odd or even row, or already
                }
                window.clearTimeout(expandData.hoverTimeout);
                expandData.hoverTimeout = window.setTimeout(() => {
                    if (!tr.hasClass('loaded')) { // could have been clicked on
                        const row = dataTable.row( tr );
                        const childHtml = expandFn(row.data());
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

        setupResponsiveExpand: function() {
            if (!this.serverParams.responsive) {
                return;
            }
            const dt = this.dataTable;
            let lastShown = null;

            dt.on('responsive-display', function (e, datatable, row, showHide, update ) {
                if (showHide) {
                    if (lastShown) {
                        // just calling hide keeps the 'parent' class
                        lastShown.nodes().to$().trigger('click');
                    }
                    lastShown = row;
                } else {
                    lastShown = null;
                }
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
                    this.setupResponsiveExpand();
                    // note the below causes the dataTable to re-download data from the server redundantly
                    // not sure under what circumstances it's actually required
                    if (this.adjustColumns) {
                        this.dataTable.columns.adjust().draw(false);
                    }
                });
            });
        }
    };

    return DataTableDefinition;
})();


// ******************************************************************************************
// COMMON TABLE FORMATTERS
// ******************************************************************************************


const TableFormat = (function() {
    "use strict";
    const TableFormat = function() {};
    TableFormat.prototype = {};
    return TableFormat;
})();

TableFormat.timestamp = (data, type, row) => {
    if (data) {
        const timestampStr = convertTimestamp(data);
        return $('<span>', {class:'timestamp', text: timestampStr}).prop('outerHTML');
    } else {
        return '';
    }
};

TableFormat.timestampSeconds = (data, type, row) => {
    if (data) {
        const timestampStr = moment(Number(data) * 1000).format(JS_DATE_FORMAT_SECONDS);
        return $('<span>', {class:'timestamp', 'text': timestampStr}).prop('outerHTML');
    } else {
        return '';
    }
};

TableFormat.timestampMilliseconds = (data, type, row) => {
    if (data) {
        const momentValue = moment(Number(data) * 1000)
        const timestampStr = momentValue.format(JS_DATE_FORMAT_SCIENTIFIC);
        const seconds = momentValue.format("ss")
        const milliseconds = momentValue.format("SSS")
        return $('<span>', {class:'timestamp', 'html': [
                timestampStr + ":",
                $('<span>', {class:'seconds', text:seconds}),
                ".",
                $('<span>', {class:'milliseconds', text:milliseconds})
            ]}).prop('outerHTML');
    } else {
        return '';
    }
};

TableFormat.list_codes = (data, type, row) => {
    if (!data) {
        return $("<span>", {text: "-"});
    }
    const elements = [];
    let isFirst = true;
    for (value of data) {
        if (!isFirst) {
            elements.push(", ");
        }
        isFirst = false;
        elements.push($('<span>', {class: 'text-monospace text-secondary', text: value}));
    }
    return $('<div>', {html: elements});
}

TableFormat.sizeBytes = (data, type, row) => {
    if (data) {
        let unit = 'bytes';
        let value = data;
        if (value > 1024) {
            value = value / 1024;
            unit = 'KB';
            if (value > 1024) {
                value = value / 1024;
                unit = 'MB'
            }
        }
        value = Math.round(value);

        return $('<span>', {class:'text-monospace', text: `${value} ${unit}`}).prop('outerHTML');
    } else {
        return '';
    }
};

TableFormat.timeAgo = (data, type, row) => {
    if (data) {
        return $('<data>', {class:'convert-timestamp time-ago', 'data-timestamp':data, text:data}).prop('outerHTML');
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

TableFormat.plain = (data, type, row) => {
    return data;
}

TableFormat.number = (data, type, row) => {
    if (data === '' || data === null) {
        return $('<span/>', {class:'no-value', text:'-'}).prop('outerHTML');
    } else {
        // TODO, put in format commas, monospace etc
        return `<span class="text-number">${data.toLocaleString('en-US')}</span>`;
    }
}

TableFormat.linkUrl = (data, type, row) => {
    if (!data) {
        return '<span class="no-value">-</span>';
    }
    const text = data.text;

    let textDom;
    if (!text) {
        textDom = $('<span>', {class: 'no-value', text:'<blank>'});
    } else {
        textDom = $('<span>', {text: text});
    }
    const aDom = $('<a>', {href:data.url, html: textDom, class: 'hover-link'});

    return aDom.prop('outerHTML');
};

TableFormat.preview = (columns, data, type, row) => {
    let dom = $('<div>');
    let hasValue = false;
    for (const col of columns) {
        let value = data[col] || row[col];
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
    } else if (style === 'false_is_error') {
        return data ? '<i class="fas fa-check-circle text-success"></i>' : '<i class="fas fa-times-circle text-danger"></i>';
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

TableFormat.combine = function(formatters, settings, data, type, columns) {
    if (settings === null) {
        settings = {};
    }
    const dom = $('<div>');
    formatters.forEach((formatter, index) => {
        let part;
        if (settings.dataMode === "combined") {
            part = eval(formatter)(data, type, columns);
        } else {
            part = eval(formatter)(data[index], type, columns);
        }
        if (settings.separator) {
            dom.append($(separator));
            dom.append(part);
        } else {
            dom.append($('<div>', {'html': part}));
        }
    });
    return dom;
};

TableFormat.repeat = function(settings, data, type, columns) {
    if (settings === null) {
        settings = {};
    }
    if (data === null) {
        return "";
    }
    const subFormatter = eval(settings.formatter);
    let cssClass = "repeat";
    if (settings.groupCSS) {
        cssClass = settings.groupCSS;
    }
    const dom = $('<div>', {"class": cssClass});
    data.forEach((subData, index) => {
        const subDom = subFormatter(subData, type, columns);
        dom.append(subDom);
    });
    return dom;
};

TableFormat.json = function(data, type, columns) {
    // TODO format JSON
    return JSON.stringify(data);
};

TableFormat.expandAjax = function(url_or_method, param, expectedHeight, data) {
    if (data) {
        let dataId = data[param];
        if (typeof(dataId) === "object") {
            dataId = dataId.id;
        }
        if (!dataId) {
            return `<i class="fas fa-bomb text-danger"></i> No value for "${param}" in this ${JSON.stringify(data)} : DEBUG - is ${param} a column in this table, visible or otherwise?`;
        }
        let reverseUrl = window[url_or_method] || Urls[url_or_method];
        if (!reverseUrl) {
            return `<i class="fas fa-bomb text-danger"></i> Method or URL not configured for "${url_or_method} : Developer may need to run<br/>
            <div class="code">python3 manage.py collectstatic_js_reverse</div>`;
        }
        if (param) {
            reverseUrl = reverseUrl(dataId);
            if (!reverseUrl) {
                throw `DataId ${dataId} did not make a URL from ${url_or_method}`;
            }
        }

        const ajaxDom =
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

TableFormat.deleteRow = function(data, type, row) {
    if (!data) {
        return '';
    }
    return $('<button>', {
        type: 'button',
        class: 'btn btn-sm btn-danger dt-delete-row',
        'data-url': data,
        html: $('<i>', {class: 'fas fa-trash'})
    }).prop('outerHTML');
};

$(document).on('click', '.dt-delete-row', function() {
    if (!confirm('Are you sure you want to delete this?')) {
        return;
    }
    const btn = $(this);
    const url = btn.data('url');
    $.ajax({
        type: 'POST',
        url: url,
        headers: {'X-CSRFToken': Cookies.get('csrftoken')},
        success: function() {
            btn.closest('table').DataTable().ajax.reload(null, false);
        },
        error: function(xhr) {
            alert('Error: ' + (xhr.responseText || 'Failed to delete'));
        }
    });
});

TableFormat.detailRendererHtml = function ( api, rowIdx, columns ) {
    const fieldset = $('<div>', {class:'mt-3'});
    for (const col of columns) {
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

// *******
// Jump to page
// https://github.com/DataTables/Plugins/blob/master/pagination/input.js
// *******

(function ($) {
	function calcDisableClasses(oSettings) {
		const start = oSettings._iDisplayStart;
		const length = oSettings._iDisplayLength;
		const visibleRecords = oSettings.fnRecordsDisplay();
		const all = length === -1;

		// Gordey Doronin: Re-used this code from main jQuery.dataTables source code. To be consistent.
		const page = all ? 0 : Math.ceil(start / length);
		const pages = all ? 1 : Math.ceil(visibleRecords / length);

		const disabledClass = oSettings.oClasses.sPageButtonDisabled;
		const disableFirstPrevClass = (page > 0 ? '' : disabledClass);
		const disableNextLastClass = (page < pages - 1 ? '' : disabledClass);

		return {
			'first': disableFirstPrevClass,
			'previous': disableFirstPrevClass,
			'next': disableNextLastClass,
			'last': disableNextLastClass
		};
	}

	function calcCurrentPage(oSettings) {
		return Math.ceil(oSettings._iDisplayStart / oSettings._iDisplayLength) + 1;
	}

	function calcPages(oSettings) {
		return Math.ceil(oSettings.fnRecordsDisplay() / oSettings._iDisplayLength);
	}

	const firstClassName = 'first';
	const previousClassName = 'previous';
	const nextClassName = 'next';
	const lastClassName = 'last';

	const paginateClassName = 'paginate';
	const paginatePageClassName = 'paginate_page';
	const paginateInputClassName = 'paginate_input';
	const paginateTotalClassName = 'paginate_total';

	$.fn.dataTableExt.oPagination.input = {
		'fnInit': function (oSettings, nPaging, fnCallbackDraw) {
			const nFirst = document.createElement('span');
			const nPrevious = document.createElement('span');
			const nNext = document.createElement('span');
			const nLast = document.createElement('span');
			const nInput = document.createElement('input');
			const nTotal = document.createElement('span');
			const nInfo = document.createElement('span');

			const language = oSettings.oLanguage.oPaginate;
			const classes = oSettings.oClasses;
			let info = language.info || 'Page _INPUT_ of _TOTAL_';

			nFirst.innerHTML = '<i class="fas fa-fast-backward"></i>'; // language.sFirst;
			nPrevious.innerHTML = '<i class="fas fa-step-backward"></i>'; //language.sPrevious;
			nNext.innerHTML = '<i class="fas fa-step-forward"></i>'; // language.sNext;
			nLast.innerHTML = '<i class="fas fa-fast-forward"></i>'; //language.sLast;

			nFirst.className = firstClassName + ' ' + classes.sPageButton;
			nPrevious.className = previousClassName + ' ' + classes.sPageButton;
			nNext.className = nextClassName + ' ' + classes.sPageButton;
			nLast.className = lastClassName + ' ' + classes.sPageButton;

			nInput.className = paginateInputClassName + " form-control d-inline-block";
			nInput.style.width = "50px";
			nTotal.className = paginateTotalClassName;

			if (oSettings.sTableId !== '') {
				nPaging.setAttribute('id', oSettings.sTableId + '_' + paginateClassName);
				nFirst.setAttribute('id', oSettings.sTableId + '_' + firstClassName);
				nPrevious.setAttribute('id', oSettings.sTableId + '_' + previousClassName);
				nNext.setAttribute('id', oSettings.sTableId + '_' + nextClassName);
				nLast.setAttribute('id', oSettings.sTableId + '_' + lastClassName);
			}

			nInput.type = 'text';

			info = info.replace(/_INPUT_/g, '</span>' + nInput.outerHTML + '<span>');
			info = info.replace(/_TOTAL_/g, '</span>' + nTotal.outerHTML + '<span>');
			nInfo.innerHTML = '<label>' + info + '</label>';

			nPaging.appendChild(nFirst);
			nPaging.appendChild(nPrevious);
			$(nInfo).children().each(function (i, n) {
			    nPaging.appendChild(n);
			});
			nPaging.appendChild(nNext);
			nPaging.appendChild(nLast);

			$(nFirst).click(function() {
				const iCurrentPage = calcCurrentPage(oSettings);
				if (iCurrentPage !== 1) {
					oSettings.oApi._fnPageChange(oSettings, 'first');
					fnCallbackDraw(oSettings);
				}
			});

			$(nPrevious).click(function() {
				const iCurrentPage = calcCurrentPage(oSettings);
				if (iCurrentPage !== 1) {
					oSettings.oApi._fnPageChange(oSettings, 'previous');
					fnCallbackDraw(oSettings);
				}
			});

			$(nNext).click(function() {
				const iCurrentPage = calcCurrentPage(oSettings);
				if (iCurrentPage !== calcPages(oSettings)) {
					oSettings.oApi._fnPageChange(oSettings, 'next');
					fnCallbackDraw(oSettings);
				}
			});

			$(nLast).click(function() {
				const iCurrentPage = calcCurrentPage(oSettings);
				if (iCurrentPage !== calcPages(oSettings)) {
					oSettings.oApi._fnPageChange(oSettings, 'last');
					fnCallbackDraw(oSettings);
				}
			});

			$(nPaging).find('.' + paginateInputClassName).keyup(function (e) {
				// 38 = up arrow, 39 = right arrow
				if (e.which === 38 || e.which === 39) {
					this.value++;
				}
				// 37 = left arrow, 40 = down arrow
				else if ((e.which === 37 || e.which === 40) && this.value > 1) {
					this.value--;
				}

				if (this.value === '' || this.value.match(/[^0-9]/)) {
					/* Nothing entered or non-numeric character */
					this.value = this.value.replace(/[^\d]/g, ''); // don't even allow anything but digits
					return;
				}

				let iNewStart = oSettings._iDisplayLength * (this.value - 1);
				if (iNewStart < 0) {
					iNewStart = 0;
				}
				if (iNewStart >= oSettings.fnRecordsDisplay()) {
					iNewStart = (Math.ceil((oSettings.fnRecordsDisplay()) / oSettings._iDisplayLength) - 1) * oSettings._iDisplayLength;
				}

				oSettings._iDisplayStart = iNewStart;
				oSettings.oInstance.trigger("page.dt", oSettings);
				fnCallbackDraw(oSettings);
			});

			// Take the brutal approach to cancelling text selection.
			$('span', nPaging).bind('mousedown', function () { return false; });
			$('span', nPaging).bind('selectstart', function() { return false; });

			// If we can't page anyway, might as well not show it.
			const iPages = calcPages(oSettings);
			if (iPages <= 1) {
				$(nPaging).hide();
			}
		},

		'fnUpdate': function (oSettings) {
			if (!oSettings.aanFeatures.p) {
				return;
			}

			const iPages = calcPages(oSettings);
			const iCurrentPage = calcCurrentPage(oSettings);

			const an = oSettings.aanFeatures.p;
			if (iPages <= 1) // hide paging when we can't page
			{
				$(an).hide();
				return;
			}

			const disableClasses = calcDisableClasses(oSettings);

			$(an).show();

			// Enable/Disable `first` button.
			$(an).children('.' + firstClassName)
				.removeClass(oSettings.oClasses.sPageButtonDisabled)
				.addClass(disableClasses[firstClassName]);

			// Enable/Disable `prev` button.
			$(an).children('.' + previousClassName)
				.removeClass(oSettings.oClasses.sPageButtonDisabled)
				.addClass(disableClasses[previousClassName]);

			// Enable/Disable `next` button.
			$(an).children('.' + nextClassName)
				.removeClass(oSettings.oClasses.sPageButtonDisabled)
				.addClass(disableClasses[nextClassName]);

			// Enable/Disable `last` button.
			$(an).children('.' + lastClassName)
				.removeClass(oSettings.oClasses.sPageButtonDisabled)
				.addClass(disableClasses[lastClassName]);

			// Paginate of N pages text
			$(an).find('.' + paginateTotalClassName).html(iPages);

			// Current page number input value
			$(an).find('.' + paginateInputClassName).val(iCurrentPage);
		}
	};
})(jQuery);