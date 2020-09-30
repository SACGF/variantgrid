/**
 * @license jqGrid Swedish Translation
 * Harald Normann harald.normann@wts.se, harald.normann@gmail.com
 * http://www.worldteamsoftware.com
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
		recordtext: "Visar {0} - {1} av {2}",
		emptyrecords: "Det finns inga poster att visa",
		loadtext: "Laddar...",
		pgtext: "Sida {0} av {1}",
		pgfirst: "First Page",
		pglast: "Last Page",
		pgnext: "Next Page",
		pgprev: "Previous Page",
		pgrecs: "Records per Page",
		showhide: "Toggle Expand Collapse Grid",
		savetext: "Sparar..."
	},
	search: {
		caption: "Sök Poster - Ange sökvillkor",
		Find: "Sök",
		Reset: "Nollställ Villkor",
		odata: [
			{ oper: "eq", text: "lika" },
			{ oper: "ne", text: "ej lika" },
			{ oper: "lt", text: "mindre" },
			{ oper: "le", text: "mindre eller lika" },
			{ oper: "gt", text: "större" },
			{ oper: "ge", text: "större eller lika" },
			{ oper: "bw", text: "börjar med" },
			{ oper: "bn", text: "börjar inte med" },
			{ oper: "in", text: "tillhör" },
			{ oper: "ni", text: "tillhör inte" },
			{ oper: "ew", text: "slutar med" },
			{ oper: "en", text: "slutar inte med" },
			{ oper: "cn", text: "innehåller" },
			{ oper: "nc", text: "innehåller inte" },
			{ oper: "nu", text: "is null" },
			{ oper: "nn", text: "is not null" }
		],
		groupOps: [
			{ op: "AND", text: "alla" },
			{ op: "OR",  text: "eller" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Click to select search operation.",
		resetTitle: "Reset Search Value"
	},
	edit: {
		addCaption: "Ny Post",
		editCaption: "Redigera Post",
		bSubmit: "Spara",
		bCancel: "Avbryt",
		bClose: "Stäng",
		saveData: "Data har ändrats! Spara förändringar?",
		bYes: "Ja",
		bNo: "Nej",
		bExit: "Avbryt",
		msg: {
			required: "Fältet är obligatoriskt",
			number: "Välj korrekt nummer",
			minValue: "värdet måste vara större än eller lika med",
			maxValue: "värdet måste vara mindre än eller lika med",
			email: "är inte korrekt e-post adress",
			integer: "Var god ange korrekt heltal",
			date: "Var god ange korrekt datum",
			url: "är inte en korrekt URL. Prefix måste anges ('http://' or 'https://')",
			nodefined: " är inte definierad!",
			novalue: " returvärde måste anges!",
			customarray: "Custom funktion måste returnera en vektor!",
			customfcheck: "Custom funktion måste finnas om Custom kontroll sker!"
		}
	},
	view: {
		caption: "Visa Post",
		bClose: "Stäng"
	},
	del: {
		caption: "Radera",
		msg: "Radera markerad(e) post(er)?",
		bSubmit: "Radera",
		bCancel: "Avbryt"
	},
	nav: {
		edittext: "",
		edittitle: "Redigera markerad rad",
		addtext: "",
		addtitle: "Skapa ny post",
		deltext: "",
		deltitle: "Radera markerad rad",
		searchtext: "",
		searchtitle: "Sök poster",
		refreshtext: "",
		refreshtitle: "Uppdatera data",
		alertcap: "Varning",
		alerttext: "Ingen rad är markerad",
		viewtext: "",
		viewtitle: "Visa markerad rad",
		savetext: "",
		savetitle: "Save row",
		canceltext: "",
		canceltitle: "Cancel row editing"
	},
	col: {
		caption: "Välj Kolumner",
		bSubmit: "OK",
		bCancel: "Avbryt"
	},
	errors: {
		errcap: "Fel",
		nourl: "URL saknas",
		norecords: "Det finns inga poster att bearbeta",
		model: "Antal colNames <> colModel!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ",", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0,00" },
		currency: { decimalSeparator: ",", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: "Kr", defaultValue: "0,00" },
		date: {
			dayNames:   [
				"Sön", "Mån", "Tis", "Ons", "Tor", "Fre", "Lör",
				"Söndag", "Måndag", "Tisdag", "Onsdag", "Torsdag", "Fredag", "Lördag"
			],
			monthNames: [
				"Jan", "Feb", "Mar", "Apr", "Maj", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dec",
				"Januari", "Februari", "Mars", "April", "Maj", "Juni", "Juli", "Augusti", "September", "Oktober", "November", "December"
			],
			AmPm: ["fm", "em", "FM", "EM"],
			S: function (j) {
				return j < 11 || j > 13 ? ["st", "nd", "rd", "th"][Math.min((j - 1) % 10, 3)] : "th";
			},
			srcformat: "Y-m-d",
			newformat: "Y-m-d",
			masks: {
				ShortDate:  "n/j/Y",
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
		locale: "sv"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		sv: $.extend({}, locInfo, { name: "svenska", nameEnglish: "Swedish" }),
		"sv-SE": $.extend({}, locInfo, { name: "svenska (Sverige)", nameEnglish: "Swedish (Sweden)" })
	}
});
}));
