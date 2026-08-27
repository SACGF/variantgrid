function renderGeneSymbol(geneSymbol, type, row) {
    let link = "";
    if (geneSymbol) {
        const linkObj = $('<a>', {
            href: Urls.view_gene_symbol(geneSymbol),
            class: 'hover-link',
            html: [
                $('<span>', {text: geneSymbol}),
            ]
        });
        link = linkObj.prop('outerHTML');
    }
    return link;
}

function renderAnalysisAuditLogSummary(summary, type, row) {
    const summaryText = summary['summary_text'];
    if (summaryText) {
        return "<span>" + summaryText + "</span>";
    }
    const changes = summary['changes'];
    if (changes) {
        const hideValues = new Set(['valid', 'status', 'version', 'shadow_color', 'appearance_version']);
        let changesSummary = "<table class='table'>";
        changesSummary += "<tr><th>field</th><th>old</th><th>new</th></tr>";
        for (const [key, value] of Object.entries(changes)) {
            if (!hideValues.has(key)) {
                changesSummary += `<tr><td>${key}</td><td>${value[0]}</td><td>${value[1]}</td>`;
            }
        }
        changesSummary += "</table>";
        return changesSummary;
    }

}


function renderExpandAnalysisAuditLogEntry(x) {
    return formatJson({
        "changes": x["changes"],
        "additional_data": x["additional_data"],
    });
}


// Server sends {text, url} - url is null where there's nothing to link to
function renderOptionalLink(data, type, row) {
    if (!data || !data.text) {
        return '<span class="no-value">-</span>';
    }
    if (!data.url) {
        return $('<span>', {text: data.text}).prop('outerHTML');
    }
    return $('<a>', {href: data.url, class: 'hover-link', text: data.text}).prop('outerHTML');
}


function renderSampleGeneListCount(data, type, row) {
    const dom = $('<span>');
    if (data.active) {
        $('<div>', {class: 'left icon16 check-mark-green', title: 'Active GeneList'}).appendTo(dom);
    }
    $('<span>', {text: data.count}).appendTo(dom);
    return dom.prop('outerHTML');
}


// Server sends the URL to POST to, or null where the user can't clone
function renderCloneRow(data, type, row) {
    if (!data) {
        return '';
    }
    return $('<button>', {
        type: 'button',
        class: 'btn btn-sm btn-outline-primary dt-clone-row',
        title: 'Clone',
        'data-url': data,
        html: $('<i>', {class: 'fas fa-copy'})
    }).prop('outerHTML');
}


$(document).on('click', '.dt-clone-row', function() {
    const btn = $(this);
    $.ajax({
        type: 'POST',
        url: btn.data('url'),
        headers: {'X-CSRFToken': Cookies.get('csrftoken')},
        success: function() {
            btn.closest('table').DataTable().ajax.reload(null, false);
        },
        error: function(xhr) {
            alert('Error: ' + (xhr.responseText || 'Failed to clone'));
        }
    });
});


// Server sends a list of {label, icon, url, css_class} - hidden behind a per-row toggle so the grid
// isn't cluttered with buttons. The menu is positioned fixed so it escapes the scrollX overflow clip.
function renderRowActions(data, type, row) {
    if (!data || !data.length) {
        return '';
    }
    const menu = $('<div>', {class: 'dropdown-menu dt-row-actions-menu'});
    for (const action of data) {
        $('<a>', {
            href: 'javascript:void(0)',
            class: `dropdown-item ${action.css_class || ''}`,
            'data-url': action.url,
            html: [
                $('<i>', {class: `${action.icon} fa-fw mr-2`}),
                $('<span>', {text: action.label})
            ]
        }).appendTo(menu);
    }
    return $('<div>', {
        class: 'dt-row-actions', html: [
            $('<button>', {
                type: 'button',
                class: 'btn btn-sm btn-outline-secondary dt-row-actions-toggle',
                title: 'Actions',
                html: $('<i>', {class: 'fas fa-ellipsis-v'})
            }),
            menu
        ]
    }).prop('outerHTML');
}

function closeRowActionsMenus() {
    $('.dt-row-actions-menu.show').removeClass('show');
}

$(document).on('click', '.dt-row-actions-toggle', function(event) {
    event.stopPropagation();
    const menu = $(this).siblings('.dt-row-actions-menu');
    const wasOpen = menu.hasClass('show');
    closeRowActionsMenus();
    if (!wasOpen) {
        const toggleRect = this.getBoundingClientRect();
        menu.css({position: 'fixed', margin: 0, top: `${toggleRect.bottom}px`, left: 'auto',
                  right: `${window.innerWidth - toggleRect.right}px`}).addClass('show');
    }
});

$(document).on('click', closeRowActionsMenus);
$(window).on('resize', closeRowActionsMenus);
// scroll doesn't bubble, so capture it to catch the DataTables scroll body as well as the window
document.addEventListener('scroll', closeRowActionsMenus, true);


// Server sends a boolean - only draw an icon when it's set, so the grid stays sparse
function renderIconFlag(settings, data, type, row) {
    if (!data) {
        return '';
    }
    return $('<div>', {class: settings.cssClass, title: settings.title}).prop('outerHTML');
}


// Server sends a list of {label, url}
function renderExternalLinks(data, type, row) {
    if (!data || !data.length) {
        return '';
    }
    const dom = $('<div>');
    data.forEach((link, index) => {
        if (index) {
            dom.append(' | ');
        }
        $('<a>', {href: link.url, target: '_blank', rel: 'noopener', text: link.label}).appendTo(dom);
    });
    return dom.prop('outerHTML');
}


// Server sends one entry per VCF loaded from the run - {id, url, import_status, variant_caller}
function renderSequencingRunVCFs(data, type, row) {
    if (!data || !data.length) {
        return '';
    }
    const dom = $('<div>', {class: 'sequencing-run-vcfs'});
    for (const vcf of data) {
        let icon;
        if (vcf.import_status === 'S') {
            icon = $('<div>', {class: 'grid-link-icon vcf-icon', title: vcf.variant_caller});
        } else if (vcf.import_status === 'E') {
            icon = $('<i>', {class: 'fas fa-times-circle text-danger', title: vcf.variant_caller});
        } else if (vcf.import_status === 'C' || vcf.import_status === 'I') {
            icon = $('<i>', {class: 'fas fa-spinner fa-spin', title: vcf.variant_caller});
        } else {
            icon = $('<span>', {text: `VCF ${vcf.id} (${vcf.import_status})`});
        }
        $('<a>', {href: vcf.url, html: icon}).appendTo(dom);
    }
    return dom.prop('outerHTML');
}


// Variant tags grid (@see VariantTagsColumns) - the tag colour CSS comes from render_tag_styles_and_formatter
function renderVariantTagVariant(data, type, row) {
    if (!data || !data.variant_string) {
        return '';
    }
    const dom = $('<div>');
    if (data.url) {
        $('<a>', {href: data.url, target: '_blank', text: data.variant_string}).appendTo(dom);
    } else {
        $('<span>', {text: data.variant_string}).appendTo(dom);
    }
    if (data.classify_url) {
        $('<a>', {
            href: data.classify_url,
            target: '_blank',
            class: 'btn btn-primary new-classification-button',
            html: [$('<i>', {class: 'fas fa-plus-circle'}), ' New Classification']
        }).appendTo(dom);
    }
    return dom.prop('outerHTML');
}


function renderVariantTagAnalysis(data, type, row) {
    if (!data || !data.text) {
        return '';
    }
    if (!data.url) {
        return $('<span>', {text: data.text}).prop('outerHTML');
    }
    return $('<a>', {href: data.url, target: '_blank', text: data.text}).prop('outerHTML');
}


function renderVariantTagPill(data, type, row) {
    if (!data || !data.tag) {
        return '';
    }
    const tag = data.tag;
    return $('<span>', {
        class: `grid-tag tagged-${tag}`,
        title: `Tagged as ${tag}`,
        variant_id: data.variant_id,
        tag_id: tag,
        html: $('<span>', {class: 'user-tag-colored', text: tag})
    }).prop('outerHTML');
}


// Comma separated list of symbols, each opening in a new window
function renderGeneSymbolNewWindow(data, type, row) {
    if (!data) {
        return '';
    }
    const dom = $('<div>');
    data.split(",").forEach((geneSymbol, index) => {
        if (index) {
            dom.append(', ');
        }
        $('<a>', {
            href: Urls.view_gene_symbol(geneSymbol),
            target: '_blank',
            title: 'View gene in new window',
            text: geneSymbol
        }).appendTo(dom);
    });
    return dom.prop('outerHTML');
}
