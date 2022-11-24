/*jshint esversion: 6 */
/*globals $:false, jQuery:false, console:false, moment:false, Rollbar:false */
/*globals toMarkdown:false */

// Required until we can up our minimum browser suppert to this https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/String/replaceAll
if (!String.prototype.replaceAll) {
	String.prototype.replaceAll = function(str, newStr){
		// If a regex pattern
		if (Object.prototype.toString.call(str).toLowerCase() === '[object regexp]') {
			return this.replace(str, newStr);
		}
		// If a string
		return this.replace(new RegExp(str, 'g'), newStr);
	};
}

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
            return $(this).attr('popover-header') || 'Help';
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
        // by providing no function it means anything tht is one of these will be skipped
        {test: `[role="gridcell"]`, func: null},
        {test: `.select2-selection__clear`, func: null},

        // load ajax blocks as soon as we see them
        {test: '[data-toggle="ajax"]', func: (node) => {loadAjaxBlock(node);}},

        {test: '[data-toggle="ajax-modal"]', func: (node) => {
            node.addClass('modal-link');
            node.click(function() {
                loadAjaxModal($(this));
                return false;
            });
        }},

        {test: '[data-toggle="ajax-collapse"]', func: (node) => {
            let $node = $(node);
            let href = $node.attr('href');
            let dataId = $node.attr('data-id');
            let title = $node.attr('title');
            let toggleLink = $(`<a data-toggle="collapse" class="toggle-link" href="#${dataId}">Toggle ${title}</a>`);
            let ajaxBlob = $(`<div class="collapse mt-2" id="${dataId}"><div class="loading-message">Loading ${title}</div></div>`);
            $node.replaceWith($('<div>', {html: [toggleLink, ajaxBlob]}));

            window.setTimeout(() => {
                ajaxBlob.on('show.bs.collapse', () => {
                    if ($node.attr('loading')) {
                        return;
                    }
                    $node.attr('loading', 1);
                    loadAjaxBlock(ajaxBlob, href);
                });
            });
        }},

        {test: '[data-help]', func: ($node) => {
            let title = $node.text();

            // wrap everything in an inline-block span
            // this is because many elements will take up 100% of horizontal space
            // so help is show far off to the right after a lot of whitespace

            // TODO, have a solution for when on a device without a mouse
            $node.wrapInner('<span></span>');
            let $target = $node.children('span').first();

            $target.css('display', 'inline-block');
            $target.addClass('hover-detail');
            $target.addClass('popover-hover-stay');
            $target.addClass('helpful');
            $target.attr('data-toggle', 'popover');
            $target.attr('title', title);
            $target.attr('data-content', $node.attr('data-help'))
            $target.attr('data-html', true);
            $target.attr('data-placement', 'left'); // top & left are preferred as most help are labels with data to the right
            $target.on("mouseenter", function () {
                let _this = this;
                $(this).popover("show");
                $(".popover").on("mouseleave", function () {
                    $(_this).popover('hide');
                });
            }).on("mouseleave", function () {
                let _this = this;
                setTimeout(function () {
                    if (!$(".popover:hover").length) {
                        $(_this).popover("hide");
                    }
                }, 300);
            });

            // remove attributes from parent element as to not get overlapping help
            $node.removeAttr('title');
            $node.removeAttr('data-help');
        }},

        // setup popovers
        {test: '[data-content]', func: (node) => {
                node.addClass('hover-detail');
                let poOpts = Object.assign({}, popoverOpts);  // clone
                if (node.hasClass("popover-hover-stay")) {
                    node.on("mouseenter", function () {
                        let _this = this;
                        $(this).popover("show");
                        $(".popover").on("mouseleave", function () {
                            $(_this).popover('hide');
                        });
                    }).on("mouseleave", function () {
                        let _this = this;
                        setTimeout(function () {
                            if (!$(".popover:hover").length) {
                                $(_this).popover("hide");
                            }
                        }, 300);
                    });
                    poOpts["trigger"] = "manual";
                }
                node.popover(poOpts);
            }
        },

        // everything with a title (that isn't data-content aka popover) give a tooltip
        {test: '[title]:not([data-content])',
            func: (node) => {
                node.tooltip({html:true, trigger : 'hover'});
                node.click(function(e) {$(this).tooltip('hide');});
            }
        },

        {test: '.nav-tabs a',
            func: (node) => {node.on('shown.bs.tab', function(e) {
                let $this = $(this);
                let url = new URL(window.location);
                let id = $this.attr('id');
                if (id.endsWith('-tab')) {
                    id = id.substring(0, id.length - 4);
                }
                url.searchParams.set('activeTab', $this.attr('data-tab-set') + ":" + id);

                window.history.replaceState({}, $this.innerHTML, url);
            })}
        },
        // load the active ajax tab now
        {test: '.nav-tabs a.active[data-href][data-toggle="tab"]',
            func: (node) => {loadAjaxTab(node);}
        },
        // load ajax tabes when they become shown
        {test: '.nav-tabs a[data-href][data-toggle="tab"]',
            func: (node) => {node.on('shown.bs.tab', function(e) {loadAjaxTab($(this));});}
        },
        // input with button at the end, have it so if you hit enter in the input, the button activates
        {test: '.input-group-append',
            func: (node) => {
                node.closest('.input-group').find('input').keydown(function (event) {
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

        {test: '.format-json', func: (node) => {
            let text = node.text().trim();
            try {
                // if not valid JSON, just print as is
                let textJson = JSON.parse(text);
                let prettyHtml = formatJson(textJson);
                prettyHtml.attr('class', node.attr('class') + ' ' + prettyHtml.attr('class'));
                prettyHtml.removeClass('format-json');
                prettyHtml.attr('data-p', 1);
                node.replaceWith(prettyHtml);
            } catch (e) {}
        }},

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
                    // putting external icon next to images breaks a few layouts
                    if (node.find('div').length === 0 && node.find('img').length === 0) {
                        node.addClass('external-link');
                    }
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
        },

        {test: 'table[data-datatable-url]',
            func: (node) => {
                let dataTableUrl = node.attr('data-datatable-url');
                let data = node.attr('data-datatable-data');
                new DataTableDefinition({
                    url: dataTableUrl,
                    data: data, // as in a function that filters the data displayed
                    filterCount: node.attr('data-datatable-filter-count'),
                    dom: node
                }).setup();
            }
        },

        {test: 'form[data-loadscreen]',
            // unsure if this kicks in in time, will have to test and see
            func: (node) => {
                let $node = $(node);
                $node.submit((event) => {
                    let selector = $node.attr('data-loadscreen');
                    let dom = $node.closest(selector) || $(selector);
                    dom.LoadingOverlay('show', {fade:false});
                    return true;
                })
            }
        },

        {test: '.current-record-menu-item',
            func: (node) => {
                let $node = $(node);
                let $moveTo = $('#current-record-spot');
                if ($moveTo.length == 0) {
                    $moveTo = $('#current-record-spot-fallback');
                }
                $node.addClass('active').detach().appendTo($moveTo);
            }
        },

        // use to have a checkbox synced to cookie (no save to database required)
        // checkbox will start as its default, but if a cookie has been set to "true" or "false" and the prop
        // checked is the opposite of that, the checkbox will be toggled to the other state (firing any change listeners)
        {test: 'input[type=checkbox][data-cookie]',
            func: (node) => {
                let $node = $(node);
                let cookieName = $node.attr('data-cookie') || $node.attr('id') || $node.attr('name');

                $node.change(() => {
                   let checked = !!$node.prop('checked');
                   Cookies.set(cookieName, checked ? 'true' : 'false', {sameSite: 'strict'});
                });

                let checked = !!$node.prop('checked') ? 'true' : 'false';
                let existingCookie = Cookies.get(cookieName);

                if (existingCookie && existingCookie != checked) {
                    $node.click();
                }
            }
        },

        {test: '[data-citation-id]',
            func: (node) => {
                CitationsManager.defaultManager.populate(node);
            }
        }

        /*
        // makes .main-icon icons in divs with the same data-group-id glow when one is highlighted
        {test: '[data-group-id]',
            func: (node) => {
                let $node = $(node);
                let dataGroupId = $node.attr('data-group-id');
                $node.mouseenter(() => {
                    console.log($(`[data-group-id='${dataGroupId}'`));
                   $(`[data-group-id='${dataGroupId}'`).addClass('group-selected');
                });
                $node.mouseleave(() => {
                    $(`[data-group-id='${dataGroupId}'`).removeClass('group-selected');
                });
            }
        }
         */
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

    let badElementTests = [];
    for (let processor of processors) {
        if (!processor.func) {
            badElementTests.push(processor.test);
        }
    }

    // update the initial state
    for (let processor of processors) {
        if (processor.func) {
            $(processor.test).each((index, node) => {
                node = $(node);
                let isGood = true;
                for (let badTest in badElementTests) {
                    if (node.is(badTest)) {
                        isGood = false;
                        break;
                    }
                }
                if (isGood) {
                    processor.func($(node));
                }
            });
        }
    }

    // would be best only to do this if added element has data-toggle, but still too much code that doens't provide that
    const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((mutation) => {
            for (let node of mutation.addedNodes) {
                node = $(node);
                if (!node.attr('data-p')) {
                    checkNode($(node), true);
                    // add an attribute at the top level to stop items being processed multiple times
                    // only do at top level so we don't mark everything with data-p, but does
                    // run the risk of adding a child of something that's already been added
                    node.attr('data-p', '1');
                }
            }
        });
    });
    observer.observe($('body')[0], { attributes: false, childList: true, subtree: true });
}

function loadAjaxModal(linkDom) {
    let url = linkDom.attr('data-href') || linkDom.attr('href');
    let useId = url.replace('/', '_');
    let modalContent = createModalShell(useId, linkDom.attr('data-title') || linkDom.text());
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
    if (!url) {
        dom.html([
            severityIcon('C'),
            "No href provided for AJAX block"
        ]);
        throw new Error("Ajax Block not given href");
    }

    let showingOverlay = false;
    // give ajax 300 ms to load before we start showing the spinner
    let spinnerTimeout = window.setTimeout(() => {
        showingOverlay = true;
        dom.LoadingOverlay('show', {zIndex: 100000});
    }, 300);

    $.ajax({
        type: "GET",
        url: url,
        async: true,
        success: (results, textStatus, jqXHR) => {
            // TODO provide the ability for a cache token, so we only reload data if something's changed
            dom.html(results);
            let autoRefreshTime = jqXHR.getResponseHeader('Auto-refresh');
            if (autoRefreshTime) {
                window.setTimeout(() => {
                    loadAjaxBlock(dom, url);
                }, parseInt(autoRefreshTime));
            }
        },
        error: (call, status, text) => {
            dom.replaceWith($('<div>', {class: 'ajax-error', html:[severityIcon('C'), "Error Loading Data"]}));
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
const JS_DATE_FORMAT = 'YYYY-MM-DD HH:mm'; //'lll';
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
        let newElement = $('<time>', {'datetime': m.toISOString(), class:'timestamp', text: m.format(JS_DATE_FORMAT)});
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

// Dialogs
function createModalShell(id, title) {
    return $(`
        <div class="modal fade" id="${id}" tabindex="-1" role="dialog" aria-labelledby="${id}Label" aria-hidden="true">
            <div class="modal-dialog modal-xl" role="document">
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

function createModal(id, title, body) {
    let modalContent = createModalShell(id, title);
    modalContent.find('.modal-body').addClass('modal-body-scroll').html(body);
    modalContent.find('.modal-footer').html(`
        <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
    `);
    let modalDialog = modalContent.modal({focus:true, show:false});
    modalContent.on('hidden.bs.modal', function() {
        modalContent.modal('dispose');
        modalContent.remove();
    });
    modalDialog.modal('show');
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

function formatJson(jsonObj) {
    return $('<div>', {class:'json', html:_formatJson(jsonObj)});
}

function _formatJson(jsonObj) {
    if (_.isNumber(jsonObj)) {
        return $('<span>', {class: 'js-num', text: jsonObj});
    } else if (_.isBoolean(jsonObj)) {
        return $('<span>', {class: 'js-bool', text: jsonObj});
    } else if (_.isString(jsonObj)) {
        let text = JSON.stringify(jsonObj);
        text = text.substring(1, text.length-1);
        let html = [];
        html.push($('<span>', {class: 'js-qt', text:"\""}));
        // TODO format "/" as special escape character
        html.push($('<span>', {class: 'js-str', text: text}));
        html.push($('<span>', {class: 'js-qt', text:"\""}));
        return $('<span>', {html:html});
        // return $('<span>', {class: 'js-str', text: JSON.stringify(jsonObj)});
    } else if (jsonObj === null) {
        return $('<span>', {class: 'js-null', text: 'null'});
    } else if (_.isArray(jsonObj)) {
        let html = [];
        html.push($('<span>', {class: 'js-br', text: '['}));
        let first = true;
        for (let elem of jsonObj) {
            if (!first) {
                html.push($('<span>', {class: 'js-comma', text: ','}));
            } else {
                first = false;
            }
            html.push(_formatJson(elem));
        }
        html.push($('<span>', {class: 'js-br', text: ']'}));
        return ($('<span>', {html: html}));
    } else if  (jsonObj && jsonObj["*wrapper$"] === "VJ") {
        let messages = jsonObj.messages;
        if (messages.length) {
            let html = [];
            let items = [];
            for (let message of messages) {
                let bsSeverity = "info";
                switch (message.severity) {
                    case "error": bsSeverity = "danger"; break;
                    case "warning": bsSeverity = "warning"; break;
                }
                items.push($('<li>', {class: `list-group-item list-group-item-${bsSeverity}`, html: [severityIcon(message.severity), message.text]}));
            }
            html.push($('<ul>', {class: 'list-group', html: items}));

            if (jsonObj.void) {
                html.push($('<div>', {class: 'js-valid-body', html: $('<span class="no-value">- entry omitted -</span>')}));
            } else {
                html.push($('<div>', {class: 'js-valid-body', html: _formatJson(jsonObj.wrap)}));
            }
            return $('<div>', {class:'js-valid', html:html});
        } else {
            // could check for void but should never have a void without messages
            return _formatJson(jsonObj.wrap);
        }
    } else {
        let html = [];
        let content = [];
        html.push($('<span>', {class: 'js-pr', text: '{'}));
        let first = true;

        for (let [key, value] of Object.entries(jsonObj)) {
            if (!first) {
                content.push($('<span>', {class: 'js-comma', text: ','}));
                content.push($('<br/>'));
            } else {
                first = false;
            }
            content.push(_formatJson(key));
            content.push($('<span>', {class: 'js-colon', text: ':'}));
            content.push(_formatJson(value));
        }
        html.push($('<span>', {class: 'js-block', html: content}));
        html.push($('<span>', {class: 'js-pr', text: '}'}));
        return ($('<span>', {class:'js-obj', html: html}));
    }
}

function diffToggle(e) {
    let diffBox = $(e).closest('.diff-box');
    let mode = diffBox.find('[name=diff-mode]:checked').val();
    diffBox.find('[data-diff-mode]').hide();
    diffBox.find(`[data-diff-mode=${mode}]`).show();
}