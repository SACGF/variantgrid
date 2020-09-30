/**
 * @license jqGrid Icelandic Translation
 * jtm@hi.is Univercity of Iceland
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
	isRTL: false,
	defaults: {
		recordtext: "Skoða {0} - {1} af {2}",
		emptyrecords: "Engar færslur",
		loadtext: "Hleður...",
		pgtext: "Síða {0} af {1}",
		pgfirst: "First Page",
		pglast: "Last Page",
		pgnext: "Next Page",
		pgprev: "Previous Page",
		pgrecs: "Records per Page",
		showhide: "Toggle Expand Collapse Grid",
		savetext: "Vistar..."
	},
	search: {
		caption: "Leita...",
		Find: "Leita",
		Reset: "Endursetja",
		odata: [
			{ oper: "eq", text: "sama og" },
			{ oper: "ne", text: "ekki sama og" },
			{ oper: "lt", text: "minna en" },
			{ oper: "le", text: "minna eða jafnt og" },
			{ oper: "gt", text: "stærra en" },
			{ oper: "ge", text: "stærra eða jafnt og" },
			{ oper: "bw", text: "byrjar á" },
			{ oper: "bn", text: "byrjar ekki á" },
			{ oper: "in", text: "er í" },
			{ oper: "ni", text: "er ekki í" },
			{ oper: "ew", text: "endar á" },
			{ oper: "en", text: "endar ekki á" },
			{ oper: "cn", text: "inniheldur" },
			{ oper: "nc", text: "inniheldur ekki" },
			{ oper: "nu", text: "is null" },
			{ oper: "nn", text: "is not null" }
		],
		groupOps: [
			{ op: "AND", text: "allt" },
			{ op: "OR",  text: "eða" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Click to select search operation.",
		resetTitle: "Reset Search Value"
	},
	edit: {
		addCaption: "Bæta við færslu",
		editCaption: "Breyta færslu",
		bSubmit: "Vista",
		bCancel: "Hætta við",
		bClose: "Loka",
		saveData: "Gögn hafa breyst! Vista breytingar?",
		bYes: "Já",
		bNo: "Nei",
		bExit: "Hætta við",
		msg: {
			required: "Reitur er nauðsynlegur",
			number: "Vinsamlega settu inn tölu",
			minValue: "gildi verður að vera meira en eða jafnt og ",
			maxValue: "gildi verður að vera minna en eða jafnt og ",
			email: "er ekki löglegt email",
			integer: "Vinsamlega settu inn tölu",
			date: "Vinsamlega setti inn dagsetningu",
			url: "er ekki löglegt URL. Vantar ('http://' eða 'https://')",
			nodefined: " er ekki skilgreint!",
			novalue: " skilagildi nauðsynlegt!",
			customarray: "Fall skal skila fylki!",
			customfcheck: "Fall skal vera skilgreint!"
		}
	},
	view: {
		caption: "Skoða færslu",
		bClose: "Loka"
	},
	del: {
		caption: "Eyða",
		msg: "Eyða völdum færslum ?",
		bSubmit: "Eyða",
		bCancel: "Hætta við"
	},
	nav: {
		edittext: "",
		edittitle: "Breyta færslu",
		addtext: "",
		addtitle: "Ný færsla",
		deltext: "",
		deltitle: "Eyða færslu",
		searchtext: "",
		searchtitle: "Leita",
		refreshtext: "",
		refreshtitle: "Endurhlaða",
		alertcap: "Viðvörun",
		alerttext: "Vinsamlega veldu færslu",
		viewtext: "",
		viewtitle: "Skoða valda færslu",
		savetext: "",
		savetitle: "Save row",
		canceltext: "",
		canceltitle: "Cancel row editing"
	},
	col: {
		caption: "Sýna / fela dálka",
		bSubmit: "Vista",
		bCancel: "Hætta við"
	},
	errors: {
		errcap: "Villa",
		nourl: "Vantar slóð",
		norecords: "Engar færslur valdar",
		model: "Lengd colNames <> colModel!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0.00" },
		currency: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0.00" },
		date: {
			dayNames:   [
				"Sun", "Mán", "Þri", "Mið", "Fim", "Fös", "Lau",
				"Sunnudagur", "Mánudagur", "Þriðjudagur", "Miðvikudagur", "Fimmtudagur", "Föstudagur", "Laugardagur"
			],
			monthNames: [
				"Jan", "Feb", "Mar", "Apr", "Maí", "Jún", "Júl", "Ágú", "Sep", "Oct", "Nóv", "Des",
				"Janúar", "Febrúar", "Mars", "Apríl", "Maí", "Júný", "Júlý", "Ágúst", "September", "Október", "Nóvember", "Desember"
			],
			AmPm: ["am", "pm", "AM", "PM"],
			S: function (j) {
				return j < 11 || j > 13 ? ["st", "nd", "rd", "th"][Math.min((j - 1) % 10, 3)] : "th";
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
		locale: "is"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		is: $.extend({}, locInfo, { name: "íslenska", nameEnglish: "Icelandic" }),
		"is-IS": $.extend({}, locInfo, { name: "íslenska (Ísland)", nameEnglish: "Icelandic (Iceland)" })
	}
});
}));
