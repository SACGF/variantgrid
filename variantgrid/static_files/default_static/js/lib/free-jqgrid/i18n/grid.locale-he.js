/**
 * @license jqGrid Hebrew Translation
 * Shuki Shukrun shukrun.shuki@gmail.com
 * http://trirand.com/blog/
 * Dual licensed under the MIT and GPL licenses:
 * http://www.opensource.org/licenses/mit-license.php
 * http://www.gnu.org/licenses/gpl.html
**/

/*jslint white: true */
/*global jQuery, module, require */
(function (factory) {
	"use strict";
	if (typeof define === "function" && define.amd) {
		// AMD. Register as an anonymous module.
		define(["jquery"], factory);
	} else if (typeof module === "object" && module.exports) {
		// Node/CommonJS
		module.exports = function (root, $) {
			if ($ === undefined) {
				// require("jquery") returns a factory that requires window to
				// build a jQuery instance, we normalize how we use modules
				// that require this pattern but the window provided is a noop
				// if it's defined (how jquery works)
				$ = typeof window !== "undefined" ?
						require("jquery") :
						require("jquery")(root || window);
			}
			factory($);
			return $;
		};
	} else {
		// Browser globals
		factory(jQuery);
	}
}(function ($) {
"use strict";
var locInfo = {
	isRTL: true,
	defaults: {
		recordtext: "מציג {0} - {1} מתוך {2}",
		emptyrecords: "אין רשומות להציג",
		loadtext: "טוען...",
		pgtext: "דף {0} מתוך {1}",
		pgfirst: "First Page",
		pglast: "Last Page",
		pgnext: "Next Page",
		pgprev: "Previous Page",
		pgrecs: "Records per Page",
		showhide: "Toggle Expand Collapse Grid",
		savetext: "שומר..."
	},
	search: {
		caption: "מחפש...",
		Find: "חפש",
		Reset: "התחל",
		odata: [
			{ oper: "eq", text: "שווה" },
			{ oper: "ne", text: "לא שווה" },
			{ oper: "lt", text: "קטן" },
			{ oper: "le", text: "קטן או שווה" },
			{ oper: "gt", text: "גדול" },
			{ oper: "ge", text: "גדול או שווה" },
			{ oper: "bw", text: "מתחיל ב" },
			{ oper: "bn", text: "לא מתחיל ב" },
			{ oper: "in", text: "נמצא ב" },
			{ oper: "ni", text: "לא נמצא ב" },
			{ oper: "ew", text: "מסתיים ב" },
			{ oper: "en", text: "לא מסתיים ב" },
			{ oper: "cn", text: "מכיל" },
			{ oper: "nc", text: "לא מכיל" },
			{ oper: "nu", text: "is null" },
			{ oper: "nn", text: "is not null" }
		],
		groupOps: [
			{ op: "AND", text: "הכל" },
			{ op: "OR",  text: "אחד מ" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Click to select search operation.",
		resetTitle: "Reset Search Value"
	},
	edit: {
		addCaption: "הוסף רשומה",
		editCaption: "ערוך רשומה",
		bSubmit: "שלח",
		bCancel: "בטל",
		bClose: "סגור",
		saveData: "נתונים השתנו! לשמור?",
		bYes: "כן",
		bNo: "לא",
		bExit: "בטל",
		msg: {
			required: "שדה חובה",
			number: "אנא, הכנס מספר תקין",
			minValue: "ערך צריך להיות גדול או שווה ל ",
			maxValue: "ערך צריך להיות קטן או שווה ל ",
			email: "היא לא כתובת איימל תקינה",
			integer: "אנא, הכנס מספר שלם",
			date: "אנא, הכנס תאריך תקין",
			url: "הכתובת אינה תקינה. דרושה תחילית ('http://' או 'https://')",
			nodefined: " is not defined!",
			novalue: " return value is required!",
			customarray: "Custom function should return array!",
			customfcheck: "Custom function should be present in case of custom checking!"
		}
	},
	view: {
		caption: "הצג רשומה",
		bClose: "סגור"
	},
	del: {
		caption: "מחק",
		msg: "האם למחוק את הרשומה/ות המסומנות?",
		bSubmit: "מחק",
		bCancel: "בטל"
	},
	nav: {
		edittext: "",
		edittitle: "ערוך שורה מסומנת",
		addtext: "",
		addtitle: "הוסף שורה חדשה",
		deltext: "",
		deltitle: "מחק שורה מסומנת",
		searchtext: "",
		searchtitle: "חפש רשומות",
		refreshtext: "",
		refreshtitle: "טען גריד מחדש",
		alertcap: "אזהרה",
		alerttext: "אנא, בחר שורה",
		viewtext: "",
		viewtitle: "הצג שורה מסומנת",
		savetext: "",
		savetitle: "Save row",
		canceltext: "",
		canceltitle: "Cancel row editing"
	},
	col: {
		caption: "הצג/הסתר עמודות",
		bSubmit: "שלח",
		bCancel: "בטל"
	},
	errors: {
		errcap: "שגיאה",
		nourl: "לא הוגדרה כתובת url",
		norecords: "אין רשומות לעבד",
		model: "אורך של colNames <> colModel!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0.00" },
		currency: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0.00" },
		date: {
			dayNames:   [
				"א", "ב", "ג", "ד", "ה", "ו", "ש",
				"ראשון", "שני", "שלישי", "רביעי", "חמישי", "שישי", "שבת"
			],
			monthNames: [
				"ינו", "פבר", "מרץ", "אפר", "מאי", "יונ", "יול", "אוג", "ספט", "אוק", "נוב", "דצמ",
				"ינואר", "פברואר", "מרץ", "אפריל", "מאי", "יוני", "יולי", "אוגוסט", "ספטמבר", "אוקטובר", "נובמבר", "דצמבר"
			],
			AmPm: ["לפני הצהרים", "אחר הצהרים", "לפני הצהרים", "אחר הצהרים"],
			S: function (j) {
				return j < 11 || j > 13 ? ["", "", "", ""][Math.min((j - 1) % 10, 3)] : "";
			},
			srcformat: "Y-m-d",
			newformat: "d/m/Y",
			masks: {
				ShortDate: "n/j/Y",
				LongDate: "l, F d, Y",
				FullDateTime: "l, F d, Y g:i:s A",
				MonthDay: "F d",
				ShortTime: "g:i A",
				LongTime: "g:i:s A",
				YearMonth: "F, Y"
			}
		}
	}
};
$.jgrid = $.jgrid || {};
$.extend(true, $.jgrid, {
	defaults: {
		locale: "he"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		he: $.extend({}, locInfo, { name: "עברית", nameEnglish: "Hebrew" }),
		"he-IL": $.extend({}, locInfo, { name: "עברית (ישראל)", nameEnglish: "Hebrew (Israel)" })
	}
});
}));
