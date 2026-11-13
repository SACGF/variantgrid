// @ts-check
/* global $ */

/*
 * CSV downloads built in the browser, for tables whose rows are all in the page already (a server
 * driven grid streams its CSV instead - DatatableConfig.server_csv_download). Cells are written the
 * way the server's exports write them (library/utils/text_utils.py:delimited_row): formula-like text
 * gets a leading quote, and the file is named <basename>_<date>.csv.
 *
 *   CsvExport.download(basename, header, rows)            arrays of cell values
 *   CsvExport.downloadDataTable(dt, basename, options)    a DataTable's rows (buttons.exportData options)
 *   CsvExport.dataTableButton(basename, options)          a DataTables button that does the above
 */
const CsvExport = (function() {
    "use strict";

    const SPREADSHEET_FORMULA_PREFIXES = ['=', '+', '-', '@', '\t', '\r'];

    /** @see library/utils/text_utils.py:spreadsheet_safe - numbers like "-1.5" are left alone */
    function spreadsheetSafe(value) {
        if (typeof value === 'string' && value.length > 1
                && SPREADSHEET_FORMULA_PREFIXES.includes(value[0]) && isNaN(Number(value))) {
            return "'" + value;
        }
        return value;
    }

    function csvCell(value, forceQuote) {
        const text = value === null || value === undefined ? '' : String(spreadsheetSafe(value));
        if (forceQuote || /[",\r\n]/.test(text)) {
            return '"' + text.replace(/"/g, '""') + '"';
        }
        return text;
    }

    function csvText(header, rows) {
        // Excel reads a file whose first letters are "ID" as a SYLK file - quoting stops that
        const quoteHeader = header.length > 0 && String(header[0]).toUpperCase().startsWith('ID');
        const lines = [header.map(value => csvCell(value, quoteHeader)).join(',')];
        for (const row of rows) {
            lines.push(row.map(value => csvCell(value, false)).join(','));
        }
        return lines.join('\r\n') + '\r\n';
    }

    function download(basename, header, rows) {
        const dateStr = new Date().toISOString().slice(0, 10);
        const blob = new Blob([csvText(header, rows)], {type: 'text/csv;charset=utf-8'});
        const url = URL.createObjectURL(blob);
        const link = $('<a>', {href: url, download: `${basename}_${dateStr}.csv`}).css('display', 'none');
        link.appendTo(document.body);
        link[0].click();
        link.remove();
        URL.revokeObjectURL(url);
    }

    /* options are DataTables' buttons.exportData() options - by default every column, every row
       matching the search, with the cell's display text */
    function downloadDataTable(dt, basename, options) {
        const data = dt.buttons.exportData(options || {});
        download(basename, data.header, data.body);
    }

    function dataTableButton(basename, options) {
        return {
            text: 'Download as CSV',
            action: function(e, dt) {
                downloadDataTable(dt, typeof basename === 'function' ? basename() : basename, options);
            },
        };
    }

    return {spreadsheetSafe, download, downloadDataTable, dataTableButton};
})();
window.CsvExport = CsvExport;
