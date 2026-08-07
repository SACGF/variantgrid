/* Tracks the async node CSV/VCF exports for an analysis (issue #1257).

   #node-data-container is replaced wholesale every time the user clicks a node, so a widget rendered
   next to the grid would be destroyed long before a big export finishes. The state therefore lives on
   the analysis toolbar (which lasts as long as the page) plus sessionStorage, and download buttons ask
   the tracker how to render themselves whenever the grid is rebuilt. */

class AnalysisDownloadTracker {
    constructor(analysisId, toolbarSelector) {
        this.storageKey = "analysis-downloads-" + analysisId;
        this.toolbarSelector = toolbarSelector;
        this.entries = this._loadFromSession();
        this.buttons = [];
        Object.values(this.entries).forEach(entry => this._pruneSupersededVersions(entry.nodeId, entry.nodeVersion));
        this.render();
        // Polling died with the page that started these, so pick them back up
        Object.values(this.entries).filter(entry => entry.state === "pending").forEach(entry => this._poll(entry));
    }

    /* Identifies a download by everything that changes its content. Node version is in here, so editing
       a node gives its exports fresh entries rather than offering a stale file. */
    static keyFor(info) {
        const IGNORED = ["page", "rows", "sidx", "sord", "nd"];
        const gridParams = {};
        Object.keys(info.gridParams || {}).sort().forEach(function (k) {
            if (IGNORED.indexOf(k) === -1) {
                gridParams[k] = info.gridParams[k];
            }
        });
        return JSON.stringify([info.nodeId, info.nodeVersion, info.exportType,
                               Boolean(info.useCanonicalTranscripts), gridParams]);
    }

    /* Called by every download button click - the tracker decides whether that means launch, ignore,
       fetch the finished file, or retry a failure. */
    download(info) {
        const key = AnalysisDownloadTracker.keyFor(info);
        const entry = this.entries[key];
        if (entry) {
            if (entry.state === "done") {
                window.location.href = entry.url;
                return;
            }
            if (entry.state === "pending") {
                return;  // Already running - the toolbar is showing its progress
            }
            // Failed - the server caches that failure against the params hash, so clear it before relaunching
            const that = this;
            this._discard(entry, function () {
                that._launch(key, info);
            });
            return;
        }
        this._launch(key, info);
    }

    /* Buttons register the element they want kept in sync, along with a function producing the export
       info. That's re-evaluated on each render, so a button follows the grid's current filters, and
       elements are dropped as their grid detaches on a node switch. */
    registerButton(element, infoFunc) {
        const el = $(element);
        if (!el.length) {
            return;
        }
        const existing = this.buttons.find(b => b.container[0] === el[0]);
        if (existing) {
            // Re-registered on a grid reload - keep the original defaultHtml, which may since have
            // been swapped out for a spinner
            existing.infoFunc = infoFunc;
            this._renderButton(existing);
            return;
        }
        const inner = el.find(".ui-pg-div");  // jqGrid nav buttons keep their padding on an inner div
        const target = inner.length ? inner : el;
        const button = {
            infoFunc: infoFunc,
            target: target,
            defaultHtml: target.html(),
            defaultTitle: el.attr("title") || "",
            container: el,
        };
        this.buttons.push(button);
        this._renderButton(button);
    }

    _launch(key, info) {
        const that = this;
        const entry = {
            key: key,
            nodeId: info.nodeId,
            nodeVersion: info.nodeVersion,
            label: info.label,
            state: "pending",
            progress: 0,
        };
        this._pruneSupersededVersions(info.nodeId, info.nodeVersion);
        this.entries[key] = entry;
        this._changed();

        $.getJSON(info.url, function (data) {
            entry.cgfId = data.cgf_id;
            entry.pollUrl = Urls.cached_generated_file_check(data.cgf_id);
            that._changed();
            that._poll(entry);
        }).fail(function () {
            that._fail(entry, "Could not start the download");
        });
    }

    _poll(entry) {
        const that = this;
        if (!entry.pollUrl) {
            this._fail(entry, "Lost track of this download");
            return;
        }
        poll_cached_generated_file(entry.pollUrl,
            function (data) {
                entry.state = "done";
                entry.progress = 1;
                entry.url = data.url;
                that._changed(true);
                window.location.href = data.url;  // They asked for the file, so give them the file
            },
            function (data) {
                that._fail(entry, (data && data.exception) || "Download failed");
            },
            function (progress) {
                entry.progress = progress || 0;
                that._changed();
            });
    }

    _fail(entry, exception) {
        entry.state = "failed";
        entry.exception = exception;
        this._changed(true);
    }

    _discard(entry, doneFunc) {
        delete this.entries[entry.key];
        this._changed();
        if (entry.cgfId) {
            $.ajax({type: "POST", url: Urls.cached_generated_file_delete(),
                    data: {cgf_id: entry.cgfId}, complete: doneFunc});
        } else {
            doneFunc();
        }
    }

    /* A node's older versions can't produce a file anyone still wants - drop them rather than offer
       a download of a superseded view */
    _pruneSupersededVersions(nodeId, nodeVersion) {
        for (const key in this.entries) {
            const entry = this.entries[key];
            if (entry.nodeId === nodeId && entry.nodeVersion < nodeVersion) {
                delete this.entries[key];
            }
        }
    }

    _loadFromSession() {
        try {
            return JSON.parse(window.sessionStorage.getItem(this.storageKey)) || {};
        } catch (e) {
            return {};
        }
    }

    _saveToSession() {
        try {
            window.sessionStorage.setItem(this.storageKey, JSON.stringify(this.entries));
        } catch (e) {
            console.log("Could not store analysis download state: " + e);
        }
    }

    _changed(pulse) {
        this._saveToSession();
        this.render(pulse);
    }

    render(pulse) {
        this._renderToolbar(pulse);
        const that = this;
        this.buttons = this.buttons.filter(function (button) {
            if (!button.container.closest("body").length) {
                return false;  // Node switched - this button's grid is gone
            }
            that._renderButton(button);
            return true;
        });
    }

    _renderToolbar(pulse) {
        const toolbar = $(this.toolbarSelector);
        if (!toolbar.length) {
            return;
        }
        const entries = Object.values(this.entries);
        if (!entries.length) {
            toolbar.hide();
            return;
        }

        const numPending = entries.filter(e => e.state === "pending").length;
        const icon = numPending ? '<i class="fas fa-spinner fa-spin"></i>' : '<i class="fas fa-download"></i>';
        const titleLines = entries.map(function (entry) {
            if (entry.state === "pending") {
                return `${entry.label} - ${Math.floor(100 * (entry.progress || 0))}%`;
            }
            if (entry.state === "failed") {
                return `${entry.label} - failed: ${entry.exception}`;
            }
            return `${entry.label} - ready`;
        });

        toolbar.html(`${icon} Downloads: ${entries.length}`);
        toolbar.attr("title", titleLines.join("\n"));
        toolbar.show();

        if (pulse) {
            toolbar.addClass("red-pulse");
            setTimeout(function () {
                toolbar.removeClass("red-pulse");
            }, 2000);
        }
    }

    _renderButton(button) {
        const info = button.infoFunc();
        const entry = this.entries[AnalysisDownloadTracker.keyFor(info)];
        const target = button.target;
        if (!entry) {
            target.html(button.defaultHtml);
            button.container.attr("title", button.defaultTitle).css("opacity", "");
            return;
        }

        if (entry.state === "pending") {
            const percent = Math.floor(100 * (entry.progress || 0));
            target.html(`<i class="fas fa-spinner fa-spin"></i> ${percent}%`);
            button.container.attr("title", `${entry.label} - preparing download`).css("opacity", 0.6);
        } else if (entry.state === "done") {
            target.html(`<i class="fas fa-download"></i> ${info.caption}`);
            button.container.attr("title", `${entry.label} - click to download`).css("opacity", "");
        } else {
            target.html(`<i class="fas fa-xmark"></i> ${info.caption}`);
            button.container.attr("title", `${entry.label} failed: ${entry.exception} - click to try again`)
                            .css("opacity", "");
        }
    }
}
