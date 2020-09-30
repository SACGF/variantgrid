/**
 * @license jqGrid Hungarian Translation
 * Őrszigety Ádám udx6bs@freemail.hu
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
	isRTL: false,
	defaults: {
		recordtext: "Oldal {0} - {1} / {2}",
		emptyrecords: "Nincs találat",
		loadtext: "Betöltés...",
		pgtext: "Oldal {0} / {1}",
		pgfirst: "First Page",
		pglast: "Last Page",
		pgnext: "Next Page",
		pgprev: "Previous Page",
		pgrecs: "Records per Page",
		showhide: "Toggle Expand Collapse Grid",
		savetext: "Mentés..."
	},
	search: {
		caption: "Keresés...",
		Find: "Keres",
		Reset: "Alapértelmezett",
		odata: [
			{ oper: "eq", text: "egyenlő" },
			{ oper: "ne", text: "nem egyenlő" },
			{ oper: "lt", text: "kevesebb" },
			{ oper: "le", text: "kevesebb vagy egyenlő" },
			{ oper: "gt", text: "nagyobb" },
			{ oper: "ge", text: "nagyobb vagy egyenlő" },
			{ oper: "bw", text: "ezzel kezdődik" },
			{ oper: "bn", text: "nem ezzel kezdődik" },
			{ oper: "in", text: "tartalmaz" },
			{ oper: "ni", text: "nem tartalmaz" },
			{ oper: "ew", text: "végződik" },
			{ oper: "en", text: "nem végződik" },
			{ oper: "cn", text: "tartalmaz" },
			{ oper: "nc", text: "nem tartalmaz" },
			{ oper: "nu", text: "is null" },
			{ oper: "nn", text: "is not null" }
		],
		groupOps: [
			{ op: "AND", text: "all" },
			{ op: "OR",  text: "any" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Click to select search operation.",
		resetTitle: "Reset Search Value"
	},
	edit: {
		addCaption: "Új tétel",
		editCaption: "Tétel szerkesztése",
		bSubmit: "Mentés",
		bCancel: "Mégse",
		bClose: "Bezárás",
		saveData: "A tétel megváltozott! Tétel mentése?",
		bYes: "Igen",
		bNo: "Nem",
		bExit: "Mégse",
		msg: {
			required: "Kötelező mező",
			number: "Kérjük, adjon meg egy helyes számot",
			minValue: "Nagyobb vagy egyenlőnek kell lenni mint ",
			maxValue: "Kisebb vagy egyenlőnek kell lennie mint",
			email: "hibás emailcím",
			integer: "Kérjük adjon meg egy helyes egész számot",
			date: "Kérjük adjon meg egy helyes dátumot",
			url: "nem helyes cím. Előtag kötelező ('http://' vagy 'https://')",
			nodefined: " nem definiált!",
			novalue: " visszatérési érték kötelező!!",
			customarray: "Custom function should return array!",
			customfcheck: "Custom function should be present in case of custom checking!"
		}
	},
	view: {
		caption: "Tétel megtekintése",
		bClose: "Bezárás"
	},
	del: {
		caption: "Törlés",
		msg: "Kiválaztott tétel(ek) törlése?",
		bSubmit: "Törlés",
		bCancel: "Mégse"
	},
	nav: {
		edittext: "",
		edittitle: "Tétel szerkesztése",
		addtext: "",
		addtitle: "Új tétel hozzáadása",
		deltext: "",
		deltitle: "Tétel törlése",
		searchtext: "",
		searchtitle: "Keresés",
		refreshtext: "",
		refreshtitle: "Frissítés",
		alertcap: "Figyelmeztetés",
		alerttext: "Kérem válasszon tételt.",
		viewtext: "",
		viewtitle: "Tétel megtekintése",
		savetext: "",
		savetitle: "Save row",
		canceltext: "",
		canceltitle: "Cancel row editing"
	},
	col: {
		caption: "Oszlopok kiválasztása",
		bSubmit: "Ok",
		bCancel: "Mégse"
	},
	errors: {
		errcap: "Hiba",
		nourl: "Nincs URL beállítva",
		norecords: "Nincs feldolgozásra váró tétel",
		model: "colNames és colModel hossza nem egyenlő!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ",", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0,00" },
		currency: { decimalSeparator: ",", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0,00" },
		date: {
			dayNames:   [
				"Va", "Hé", "Ke", "Sze", "Csü", "Pé", "Szo",
				"Vasárnap", "Hétfő", "Kedd", "Szerda", "Csütörtök", "Péntek", "Szombat"
			],
			monthNames: [
				"Jan", "Feb", "Már", "Ápr", "Máj", "Jún", "Júl", "Aug", "Szep", "Okt", "Nov", "Dec",
				"Január", "Február", "Március", "Áprili", "Május", "Június", "Július", "Augusztus", "Szeptember", "Október", "November", "December"
			],
			AmPm: ["de", "du", "DE", "DU"],
			S: function () { return ".-ik"; },
			srcformat: "Y-m-d",
			newformat: "Y/m/d",
			masks: {
				ShortDate: "Y/j/n",
				LongDate: "Y. F hó d., l",
				FullDateTime: "l, F d, Y g:i:s A",
				MonthDay: "F d",
				ShortTime: "a g:i",
				LongTime: "a g:i:s",
				YearMonth: "Y, F"
			}
		}
	}
};
$.jgrid = $.jgrid || {};
$.extend(true, $.jgrid, {
	defaults: {
		locale: "hu"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		hu: $.extend({}, locInfo, { name: "magyar", nameEnglish: "Hungarian" }),
		"hu-HU": $.extend({}, locInfo, { name: "magyar (Magyarország))", nameEnglish: "Hungarian (Hungary)" })
	}
});
}));
