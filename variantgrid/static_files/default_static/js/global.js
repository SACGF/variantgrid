/*jshint esversion: 6 */
/*globals $:false, jQuery:false, console:false, moment:false, Rollbar:false */
/*globals toMarkdown:false */

function tweakAjax() {
    // automatically send X-CSRFToken on ajax calls
    $.ajaxSetup({
        beforeSend: function(xhr, settings) {
            //if (settings.type == 'POST' || settings.type == 'PUT' || settings.type == 'DELETE' || settings.type == 'PATCH') {
                if (!(/^http:.*/.test(settings.url) || /^https:.*/.test(settings.url))) {
                    // Only send the token to relative URLs i.e. locally.
                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                } else {
                    // console.log("Didn't add CSRF! to external url '" + settings.url + "'");
                }
            //}
        }
    });

    $(document).ajaxError(function(event, jqxhr, settings, thrownError) {
        console.log(settings);
        if (settings.suppressErrors) {
            return;
        }
        // only relevant when using ODIC https://mozilla-django-oidc.readthedocs.io/en/stable/xhr.html
        if (jqxhr.status === 403 && jqxhr.getResponseHeader('refresh_url')) {
                // shouldn't happen anymore if avoid sending cookie refresh on ajax call
                alert('Session has timed out, refreshing');
                //document.location.href = jqxhr.headers.get("refresh_url");
                location.reload();
            //}
        }
        console.log("Ajax error");
        console.log(event);
        console.log(jqxhr);
        // messagePoller.stop_polling();
        // $.blockUI({ message: $('#ajax-error') });
    });
}

function enhanceAndMonitor() {
    let popoverOpts = {
        html: true,
        trigger: 'hover click',
        title: function() {
            return 'Help';
        },
        content: function() {
            return $(this).attr('title');
        }
    };

    // setup a list of processors
    // they have a test (a selector for a node to meet)
    // then an action to perform over that node if the test is met
    // done like this so can respond to dynamically added elements through ajax or complicated js
    let processors = [
        // these elements don't like being screwed with
        {test: `[role="gridcell"]`, func: null},
        {test: `.select2-selection__clear`, func: null},

        // load ajax blocks as soon as we see them
        {test: '[data-toggle="ajax"]', func: (node) => {loadAjaxBlock(node);}},

        {test: '[data-toggle="ajax-modal"]', func: (node) => {
            console.log("FOUND AJAX MODAL");
            node.click(function() {
                console.log("CLIKCED AJAX MODAL");
                loadAjaxModal($(this));
                return false;
            });
        }},

        // setup popovers
        {test: '[data-content]', func: (node) => {node.addClass('hover-detail'); node.popover(popoverOpts);}},

        // everything with a title (that isn't data-content aka popover) give a tooltip
        {test: '[title]:not([data-content])',
            func: (node) => {
                node.tooltip({html:true, trigger : 'hover'});
                node.click(function(e) {$(this).tooltip('hide');});
            }
        },
        // load the active ajax tab now
        {test: '.nav-tabs a.active[data-href][data-toggle="tab"]',
            func: (node) => {loadAjaxTab(node);}
        },
        // load ajax tabes when they become shown
        {test: '.nav-tabs a[data-href][data-toggle="tab"]',
            func: (node) => {node.on('shown.bs.tab', function(e) {loadAjaxTab($(this));});}
        },
        // input with button at the end, have it so if you hit enter in the input, hte button activates
        {test: '.input-group-append',
            func: (node) => {
                node.closest('.input-group').find('input').keyup(function (event) {
                    if (event.which === 13) {
                        let button = $(this).closest('.input-group').find('.input-group-append .btn');
                        if (button.length) {
                            event.preventDefault();
                            event.stopPropagation();
                            button.click();
                        }
                    }
                });
            }
        },
        // timestamps
        {test: '.convert-timestamp', func: (node) => { convertTimestampDom(node); }},

        // if have a wide checkbox row, make it so clicking anywhere on the row activates the checkbox
        {test: '.list-group-checkbox', func: (node) => {
            node.click(function(event) {
                $(this).find('input[type=radio]').prop("checked", true);
            });
        }},
        // similar but for radio buttons
        {test: '.radio-row', func: (node) => {node.click(event => {$(event.currentTarget).find('[type="radio"]').prop('checked', 'checked').change();});}},
        // we don't generally allow future dates
        {test: '.date-picker', func: (node) => {node.datepicker({changeYear: true, yearRange: "-120:+0"});}},

        // is this still used?? Would like to get rid of
        {test: '#id_import_status',
            func: (node) => {
                let importStatus = node.val();
                if (importStatus !== 'S') {
                    let className = null;
                    if (importStatus === 'E') {
                        className = "error";
                    } else {
                        className = "warning";
                    }
                    node.addClass(className);
                }
            }
        },

        // fixing external link styles
        {test: `a[href^="http"], a[href^="ftp"], a[href*="google.com"]`,
            func: (node) => {
                if (!node.attr('target')) {
                    node.attr('target', '_blank');
                }
                if (!node.hasClass('external-link')) {
                    node.addClass('external-link');
                }
            }
        },

        {test: 'input[name=csrfmiddlewaretoken]',
            func: (node) => {
                let form = node.closest('form');
                form.submit(function (e) {
                    // if a page load has changed the cookie, change it in the form
                    let cookie = getCookie('csrftoken');
                    node.val(cookie);
                    return true;
                });
            }
        }
    ];

    // run the processors, and check recursively
    function checkNode(node, recursive) {
        for (let processor of processors) {
            if (node.is(processor.test)) {
                if (processor.func) {
                    processor.func(node);
                } else {
                    return; // no function means it's some weird element type we should stay away from
                }
            }
        }
        // check recursively
        if (recursive) {
            for (let child of node.children()) {
                checkNode($(child), true);
            }
        }
    }

    // update the initial state
    for (let processor of processors) {
        if (processor.func) {
            $(processor.test).each((index, node) => {
                processor.func($(node));
            });
        }
    }

    // would be best only to do this if added element has data-toggle, but still too much code that doens't provide that
    const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((mutation) => {
            for (let node of mutation.addedNodes) {
                checkNode($(node), true);
            }
        });
    });
    observer.observe($('body')[0], { attributes: false, childList: true, subtree: true });
}

function loadAjaxModal(linkDom) {
    let url = linkDom.attr('href');
    let useId = url.replace('/', '_');
    let modalContent = createModalShell(useId, linkDom.attr('data-title') || 'Dialog');
    let body = modalContent.find('.modal-body');
    modalContent.find('.modal-footer').remove();
    let content = $('<div>').appendTo(body);
    let modalDialog = modalContent.modal({focus:true, show:false});
    loadAjaxBlock(content, url);

    modalContent.on('hidden.bs.modal', function() {
        modalContent.modal('dispose');
        modalContent.remove();
    });
    modalDialog.modal('show');
}

function loadAjaxBlock(dom, url) {
    if (!url) {
        url = dom.attr('href');
    }

    let showingOverlay = false;
    // give ajax 300 ms to load before we start showing the spinner
    let spinnerTimeout = window.setTimeout(() => {
        showingOverlay = true;
        dom.LoadingOverlay('show');
    }, 300);

    $.ajax({
        type: "GET",
        url: url,
        async: true,
        success: (results) => {
            dom.html(results);
        },
        error: (call, status, text) => {
            dom.replaceWith("Error Loading Data");
        },
        complete: (jqXHR, textStatus) => {
            if (showingOverlay) {
                dom.LoadingOverlay('hide');
            }
            window.clearTimeout(spinnerTimeout);
        }
    });
}

function loadAjaxTab(tab) {
    if (tab.data('loaded') !== '1') {
        tab.data('loaded', '1');
        if (!tab.attr('data-href')) {
            console.log(tab);
            alert(`Tab doesn't have a data-href`);
            return;
        }

        // is this all really required?
        let firstSlash = tab.attr('data-href').indexOf('/');
        let url = '';
        if (firstSlash !== 0) {
            url = window.location.href + '/' + tab.attr('data-href');
        } else {
            url = tab.attr('data-href');
        }
        let tabContent = $(tab.attr('href'));
        loadAjaxBlock(tabContent, url);
    }
}

function globalSetup() {
    tweakAjax();
    configureTimestamps();
    // stops there being a popup to the user
    $.fn.DataTable.ext.errMode = 'none';

    // applies many tweaks and functionality (such as ajax blocks)
    // as well as applying them to dynamically added elements
    enhanceAndMonitor();
}

// fix copied from https://github.com/allpro/layout/issues/17
(function ($){
    $.fn.selector = { split: function() { return ""; }};
})(jQuery);

function getCookie(name) {
    let match = document.cookie.match(new RegExp('(^| )' + name + '=([^;]+)'));
    if (match) {
        return match[2];
    }
    return null;
}

function setCrossLink(link_selector, urlFunc, pk) {
	if (pk) {
		link_selector.attr("href", urlFunc(pk));
        link_selector.show();
	} else {
		link_selector.hide();
	}
}

function EncodeQueryData(data, skipNulls) {
	var ret = [];
	for (let d in data) {
	    let value = data[d];
	    if (skipNulls && value == null) {
	        continue;
        }
		ret.push(encodeURIComponent(d) + "=" + encodeURIComponent(value));
	}
	return ret.join("&");
}

function dictFromLabelsAndValues(labels, values) {
    var dict = {};
    for (var i=0; i<labels.length ;i++) {
        dict[labels[i]] = values[i];
    }
    return dict;
}

function getVariantTagHtml(variantId, tag) {
    return "<span class='grid-tag tagged-" + tag + "' title='Tagged as " + tag + "' variant_id='" + variantId + "' tag_id='" + tag + "'><span class='user-tag-colored'>" + tag + "</span></span>";
}

function deleteItemClickHandler(outerElement, innerSpan, deleteClickHandler) {
    var isExpanded = innerSpan.attr("original_width");
    var completeFunc;
    var desiredWidth, desiredHeight;

    if (isExpanded) {
        desiredWidth = innerSpan.attr("original_width");
        desiredHeight = innerSpan.attr("original_height");
        completeFunc = function() {
            innerSpan.removeAttr("original_width");
            innerSpan.removeAttr("original_height");
            $(".click-to-delete-button", $(outerElement)).remove();
        };

    } else {
        innerSpan.attr("original_width", innerSpan.outerWidth());
        innerSpan.attr("original_height", innerSpan.outerHeight());

        desiredWidth = innerSpan.outerWidth() * 2;
        desiredHeight = innerSpan.outerHeight(); // * 2;
        completeFunc = function() {
            // Add [X] button to call delete
            var deleteButton = $('<span />').attr({ 'style' : "display: inline-block;",
                                                    'class' : "click-to-delete-button",
                                                    'title' : 'Delete'});
            deleteButton.click(deleteClickHandler);
            deleteButton.appendTo(innerSpan);
        };
    }
    var params = {'width' : desiredWidth, 'height' : desiredHeight};
    innerSpan.animate(params, 200, 'swing', completeFunc);
}

function getValue(val, defaultValue) {
    return (typeof val !== 'undefined') ?  val : defaultValue;
}

function markdownize(content) {
    let html = content.split("\n").map($.trim).filter(function(line) {
        return line !== "";
    }).join("\n");
    return toMarkdown(html);
}

function format(str, col) {
    col = typeof col === 'object' ? col : Array.prototype.slice.call(arguments, 1);

    return str.replace(/\{\{|\}\}|\{(\w+)\}/g, function (m, n) {
        if (m === "{{") { return "{"; }
        if (m === "}}") { return "}"; }
        return col[n];
    });
}


// From https://jsfiddle.net/salman/f9Re3/
function invertColor(hexTripletColor) {
    var color = hexTripletColor;
    color = color.substring(1); // remove #
    color = parseInt(color, 16); // convert to integer
    color = 0xFFFFFF ^ color; // invert three bytes
    color = color.toString(16); // convert to hex
    color = ("000000" + color).slice(-6); // pad with leading zeros
    color = "#" + color; // prepend #
    return color;
}

function removeItemFromArray(item, array) {
    var index = array.indexOf(item);
    if (index > -1) {
        array.splice(index, 1);
    }
}

function dynamicSort(property, caseSensitive) {
    caseSensitive = getValue(caseSensitive, true);
    let sortOrder = 1;
    if(property[0] === "-") {
        sortOrder = -1;
        property = property.substr(1);
    }
    return function (a,b) {
        let valA = a[property];
        let valB = b[property];
        if (!caseSensitive) {
            valA = valA.toUpperCase();
            valB = valB.toUpperCase();
        }

        let result = (valA < valB) ? -1 : (valA > valB) ? 1 : 0;
        return result * sortOrder;
    };
}

function zero_pad(num, size) {
    // deprecated use _.pad(num + "", size, '0');
    let s = num + "";
    while (s.length < size) s = "0" + s;
    return s;
}

function clearAutocompleteChoice(selector) {
    $(selector).val(null).trigger('change');
}

function setAutocompleteValue(selector, value, label) {
    // https://select2.org/programmatic-control/add-select-clear-items#selecting-options
    var option = new Option(label, value, true, true);
    $(selector).append(option).trigger('change');
}

// called by update_django_messages tag
function update_django_messages(messages) {
    if (messages.length === 0) {
        return;
    }
    let messagesDom = $("#django-messages");
    messagesDom.empty();
    for (let message of messages) {
        let text = message.text;
        let timestamp = "";
        if (text.indexOf("saved successfully") !== -1) {
            let m = moment();
            timestamp= $('<time>', {'class': 'float-right', 'datetime': m.toISOString(), text: m.format(JS_DATE_FORMAT_SECONDS)});
        }
        $('<div>', {class: `alert ${message.tags}`, role:"alert", html:[
            "<i class=\"fas fa-exclamation-circle\"></i>",
            $('<span>', {text: text}),
            timestamp
        ]}).appendTo(messagesDom);
    }
    messagesDom.hide();
    messagesDom.fadeIn('slow');
    $('body').scrollTop(0); // needed when in phone size
    $('.main-content').scrollTop(0); // needed when in desktop size
}

function createMessage(className, message) {
    let errorMessageUl = $("<ul/>", {class: "messages"});
    let errorMessageLi = $("<li/>", {class: "save-message"});
    errorMessageUl.append(errorMessageLi);
    errorMessageLi.addClass(className);
    errorMessageLi.html(message);
    return errorMessageUl;
}


function checkLoggedIn(loggedInHandler, loggedOutHandler) {
    // call server to check if authenticated & call appropriate handler
    $.ajax({
        type: "GET",
        url: "/authenticated",
        success: function(data) {
            if (data["authenticated"]) {
                if (loggedInHandler) {
                    loggedInHandler();
                }
            } else {
                if (loggedOutHandler) {
                    loggedOutHandler();
                }
            }
        }
    });
}

const JS_DATE_FORMAT_DETAILED = 'YYYY-MM-DD HH:mm:ss ZZ';
const JS_DATE_FORMAT_SECONDS = 'YYYY-MM-DD HH:mm:ss';
const JS_DATE_FORMAT_SCIENTIFIC = 'YYYY-MM-DD HH:mm';
const JS_DATE_FORMAT = 'lll';
function configureTimestamps() {
    $.timeago.settings.allowFuture = true;
    $.timeago.settings.strings = {
        prefixAgo: null,
        prefixFromNow: null,
        suffixAgo: "ago",
        suffixFromNow: "from now",
        seconds: "<1 min",
        minute: "1 min",
        minutes: "%d mins",
        hour: "1 hour",
        hours: "%d hours",
        day: "1 day",
        days: "%d days",
        month: "1 month",
        months: "%d months",
        year: "1 year",
        years: "%d years",
        wordSeparator: " ",
        numbers: []
    };
}
function convertTimestampDom(elem) {
    elem = $(elem);
    let unix = Number(elem.attr('data-timestamp'));
    let m = moment();
    if (unix) {
        m = moment(unix * 1000);
    }
    if (elem.hasClass('time-ago')) {
        let isFuture = moment().diff(m) < 0;
        let newElement = $('<time>', {class: 'ago', datetime: m.toISOString(), text: m.format(JS_DATE_FORMAT_DETAILED)});
        if (isFuture) {
            newElement.addClass('future');
        }
        elem.replaceWith(newElement);
        newElement.timeago();
        newElement.tooltip();
    } else if (elem.hasClass('seconds')) {
        let newElement = $('<time>', {'datetime': m.toISOString(), text: m.format(JS_DATE_FORMAT_SECONDS)});
        elem.replaceWith(newElement);
    } else {
        let newElement = $('<time>', {'datetime': m.toISOString(), text: m.format(JS_DATE_FORMAT)});
        elem.replaceWith(newElement);
    }
}
function createTimestampDom(unix, timeAgo) {
    let jsTime = unix * 1000;
    let m = moment(jsTime);
    if (timeAgo) {
        let timeAgoText = jQuery.timeago(jsTime);
        return $('<time>', {class: 'ago', datetime: m.toISOString(), text: timeAgoText, title: m.format(JS_DATE_FORMAT_DETAILED)});
    } else {
        return $('<time>', {'datetime': m.toISOString(), text: m.format(JS_DATE_FORMAT)});
    }
}
function convertTimestamp(unixDate) {
    return moment(Number(unixDate) * 1000).format(JS_DATE_FORMAT_SCIENTIFIC);
}

function blankToNull(val) {
    if (!val) {
        return null;
    } else {
        return val;
    }
}

function limitLengthSpan(text, limit) {
    if (text && text.length > limit) {
        show_text = text.substring(0, limit) + '...';
        return $('<span>', {text: show_text, title: text});
    } else {
        return $('<span>', {text: text});
    }
}
function limitLength(text, limit) {
    if (text && text.length > limit) {
        return text.substring(0, limit) + '...';
    } else {
        return text;
    }
}

function debounce( func , timeout ) {
    var timeoutID , timeout = timeout || 200;
    return function () {
        var scope = this , args = arguments;
        clearTimeout( timeoutID );
        timeoutID = setTimeout( function () {
            func.apply( scope , Array.prototype.slice.call( args ) );
        }, timeout );
   }
}

function highlightTextAsDom(value, full_text) {
    let startsAt = full_text.toLowerCase().indexOf(value.toLowerCase());
    if (startsAt != -1) {
        return $('<span>', {html: [
            $('<span>', {text: full_text.substring(0, startsAt)}),
            $('<span>', {text: full_text.substring(startsAt, startsAt + value.length), 'class': 'match'}),
            $('<span>', {text: full_text.substring(startsAt + value.length)})
        ]});
    } else {
        return $('<span>', {text: full_text});
    }
}

let TableFormat = (function() {
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
}
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
TableFormat.detailRenderer = function ( api, rowIdx, columns ) {
    console.log(api);
    let fieldset = $('<div>', {class:'mt-3'});
    for (let col of columns) {
        if (col.hidden) {
            if (col === null || col.data.length === 0) {
                // pass
            } else {
                $('<div>', {class:'row mt-2', html:[
                    $('<div>', {class: 'col-2 text-right', html:
                        $('<label>', {text: col.title})
                    }),
                    $('<div>', {class: 'col-10', html:
                        $('<span>', {class: 'dt-detail', text: col.data})
                    }),
                ]}).appendTo(fieldset);
            }
        }
    }
    return fieldset;
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
TableFormat.boolean = function(style, data, type, columns) {
    if (style == 'warning') {
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
TableFormat.expandAjax = function ajaxBlock(url, param, data) {
    if (data) {
        let dataId = data[param];
        let ajaxId = `ajax_${dataId}`
        let reverseUrl = Urls[url];
        if (!reverseUrl) {
            return `<i class="fas fa-bomb text-danger"></i> URL not configured for "${url} : Developer may need to run<br/>
            <div class="code">manage.py collectstatic_js_reverse</div>`;
        }
        if (!dataId) {
            return `<i class="fas fa-bomb text-danger"></i> No value for "${param}" in this row : Developer, is ${param} a column in this table, visible or otherwise?`;
        }
        if (param) {
            reverseUrl = reverseUrl(dataId)
        }
        window.setTimeout(() => {
            let ajaxDom = $(`#${ajaxId}`);
            $.ajax({
                type: "GET",
                url: reverseUrl,
                success: (results) => {
                    ajaxDom.LoadingOverlay('hide');
                    ajaxDom.html(results);
                },
                error: (call, status, text) => {
                    ajaxDom.LoadingOverlay('hide');
                    ajaxDom.html("Error Loading Data");
                }
            });
        }, 0);
        return `<div id="${ajaxId}">Loading</div>`;
    } else {
        return '';
    }
};

// Dialogs
function createModalShell(id, title) {
    return $(`
        <div class="modal fade" id="${id}" tabindex="-1" role="dialog" aria-labelledby="${id}Label" aria-hidden="true">
            <div class="modal-dialog" role="document">
                <div class="modal-content">
                    <div class="modal-header">
                        <h5 class="modal-title" id="${id}Label">${title}</h5>
                        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                            <span aria-hidden="true">&times;</span>
                        </button>
                    </div>
                    <div class="modal-body">
                    </div>
                    <div class="modal-footer">
                    </div>
                </div>
            </div>
        </div>
    `);
}

// Suggestions

function suggestionDialog(userName) {
    let modalDialog = window.MODAL_SUGGESTION;
    if (!modalDialog) {
        let siteName = window.SITE_NAME || 'Variant Grid';
        // FIXME need to escape username, siteName, location etc
        let modalContent = createModalShell('suggestionModal', 'Suggestion / Bug Report')
        modalContent.find('.modal-body').html(
            `<p>Thank you for taking the time to report a bug or raise a suggestion to help us improve this product.</p>
            <form>
                <div class="form-group row">
                    <label class="col-form-label col-12 col-md-3 text-md-right">Reporter</label>
                    <div class="col-12 col-md-9 text-left align-self-center">${userName}</div>
                </div>
                <div class="form-group row">
                    <label class="col-form-label col-12 col-md-3 text-md-right">Assignee</label>
                    <div class="col-12 col-md-9 text-left align-self-center">The ${siteName} Team</div>
                </div>
                <div class="form-group row">
                    <label class="col-form-label col-12 col-md-3 text-md-right">URL</label>
                    <div class="col-12 col-md-9 text-left align-self-center text-break">${window.location}</div>
                </div>
                <div class="form-group row">
                    <label class="col-form-label col-12 col-md-3 text-md-right" for="suggestion-subject">Subject</label>
                    <div class="col-12 col-md-9 text-left align-self-center">
                        <input class="form-control" id="suggestion-subject" name="form-subject">
                    </div>
                </div>
                <div class="form-group row">
                    <label class="col-form-label col-12 col-md-3 text-md-right" for="suggestion-description">Description</label>
                    <div class="col-12 col-md-9 text-left align-self-center">
                        <textarea class="form-control" id="suggestion-description" name="suggestion-description" rows="5"></textarea>
                    </div>
                </div>
            </form>`
        );
        modalContent.find('.modal-footer').html(`
            <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
            <button type="button" class="btn btn-primary" onclick="suggestionDialogSaved()">Send</button>
        `);
        modalContent.on('shown.bs.modal', function(e) {
            $('#suggestion-subject').focus();
        });
        modalDialog = modalContent.modal({focus:true, show:false});

        window.MODAL_SUGGESTION = modalDialog;
    }
    modalDialog.modal('show');
}
function suggestionDialogSaved() {
    let subjectInput = $('#suggestion-subject');
    let descriptionInput = $('#suggestion-description');

    let subjectText = subjectInput.val();
    let contentText = descriptionInput.val();

    let loadOverlayMe = $('#suggestionModal .modal-content');
    loadOverlayMe.LoadingOverlay('show');

    Rollbar.info(`User Feedback : ${subjectText}`, {subject: subjectText, content: contentText}, (err, data) => {
        loadOverlayMe.LoadingOverlay('hide');
        window.MODAL_SUGGESTION.modal('hide');

        if (err) {
            console.log(err);
            window.alert('Unable to send suggestion.');
        } else {
            subjectInput.val(null);
            descriptionInput.val(null);
            window.alert('Suggestion sent, thank you');
        }
    });
}

// FIXME turn into Bootstrap modal
function showReloadPageErrorDialog(selector, message, allowClose) {
    let buttons = [
        {   text: "Reload Page",
            class: "btn",
            click: function() {
                $(this).dialog("close");
                location.reload();
            },
        },
    ];
    if (allowClose) {
        let closeButton = {
            text: "Close and continue (not recommended)",
            class: "btn btn-outline-danger",
            click: function () {
                $(this).dialog("close");
            },
        }
        buttons.push(closeButton);
    }

    $(selector).html(message).dialog({
        dialogClass: "no-close",
        minWidth: 500,
        buttons: buttons,
    });
}

function severityIcon(severity) {
    let first = severity.toUpperCase()[0];
    switch (first) {
        case 'C': return $('<i class="fas fa-bomb text-danger"></i>'); // critical
        case 'E':  // error
        case 'D':  // danger (BootStrap)
            return $('<i class="fas fa-exclamation-circle text-danger"></i>'); // error
        case 'W': return $('<i class="fas fa-exclamation-triangle text-warning"></i>'); // warning
        case 'I': return $('<i class="fas fa-info-circle text-info"></i>'); // info
        case 'S': return $('<i class="fas fa-check-circle text-success"></i>'); // success
        default: return $(severity);
    }
}