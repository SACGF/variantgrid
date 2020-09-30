/**
 * @license jqGrid Czech Translation
 * Pavel Jirak pavel.jirak@jipas.cz
 * doplnil Thomas Wagner xwagne01@stud.fit.vutbr.cz
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
		recordtext: "Zobrazeno {0} - {1} z {2} záznamů",
		emptyrecords: "Nenalezeny žádné záznamy",
		loadtext: "Načítám...",
		pgtext: "Strana {0} z {1}",
		pgfirst: "First Page",
		pglast: "Last Page",
		pgnext: "Next Page",
		pgprev: "Previous Page",
		pgrecs: "Records per Page",
		showhide: "Toggle Expand Collapse Grid",
		savetext: "Ukládání..."
	},
	search: {
		caption: "Vyhledávám...",
		Find: "Hledat",
		Reset: "Reset",
		odata: [
			{ oper: "eq", text: "rovno" },
			{ oper: "ne", text: "nerovno" },
			{ oper: "lt", text: "menší" },
			{ oper: "le", text: "menší nebo rovno" },
			{ oper: "gt", text: "větší" },
			{ oper: "ge", text: "větší nebo rovno" },
			{ oper: "bw", text: "začíná s" },
			{ oper: "bn", text: "nezačíná s" },
			{ oper: "in", text: "je v" },
			{ oper: "ni", text: "není v" },
			{ oper: "ew", text: "končí s" },
			{ oper: "en", text: "nekončí s" },
			{ oper: "cn", text: "obsahuje" },
			{ oper: "nc", text: "neobsahuje" },
			{ oper: "nu", text: "is null" },
			{ oper: "nn", text: "is not null" }
		],
		groupOps: [
			{ op: "AND", text: "všech" },
			{ op: "OR",  text: "některého z" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Click to select search operation.",
		resetTitle: "Reset Search Value"
	},
	edit: {
		addCaption: "Přidat záznam",
		editCaption: "Editace záznamu",
		bSubmit: "Uložit",
		bCancel: "Storno",
		bClose: "Zavřít",
		saveData: "Data byla změněna! Uložit změny?",
		bYes: "Ano",
		bNo: "Ne",
		bExit: "Zrušit",
		msg: {
			required: "Pole je vyžadováno",
			number: "Prosím, vložte validní číslo",
			minValue: "hodnota musí být větší než nebo rovná ",
			maxValue: "hodnota musí být menší než nebo rovná ",
			email: "není validní e-mail",
			integer: "Prosím, vložte celé číslo",
			date: "Prosím, vložte validní datum",
			url: "není platnou URL. Vyžadován prefix ('http://' or 'https://')",
			nodefined: " není definován!",
			novalue: " je vyžadována návratová hodnota!",
			customarray: "Custom function mělá vrátit pole!",
			customfcheck: "Custom function by měla být přítomna v případě custom checking!"
		}
	},
	view: {
		caption: "Zobrazit záznam",
		bClose: "Zavřít"
	},
	del: {
		caption: "Smazat",
		msg: "Smazat vybraný(é) záznam(y)?",
		bSubmit: "Smazat",
		bCancel: "Storno"
	},
	nav: {
		edittext: "",
		edittitle: "Editovat vybraný řádek",
		addtext: " ",
		addtitle: "Přidat nový řádek",
		deltext: "",
		deltitle: "Smazat vybraný záznam ",
		searchtext: "",
		searchtitle: "Najít záznamy",
		refreshtext: "",
		refreshtitle: "Obnovit tabulku",
		alertcap: "Varování",
		alerttext: "Prosím, vyberte řádek",
		viewtext: "",
		viewtitle: "Zobrazit vybraný řádek",
		savetext: "",
		savetitle: "Save row",
		canceltext: "",
		canceltitle: "Cancel row editing"
	},
	col: {
		caption: "Zobrazit/Skrýt sloupce",
		bSubmit: "Uložit",
		bCancel: "Storno"
	},
	errors: {
		errcap: "Chyba",
		nourl: "Není nastavena url",
		norecords: "Žádné záznamy ke zpracování",
		model: "Délka colNames <> colModel!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0.00" },
		currency: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0.00" },
		date: {
			dayNames:   [
				"Ne", "Po", "Út", "St", "Čt", "Pá", "So",
				"Neděle", "Pondělí", "Úterý", "Středa", "Čtvrtek", "Pátek", "Sobota"
			],
			monthNames: [
				"Led", "Úno", "Bře", "Dub", "Kvě", "Čer", "Čvc", "Srp", "Zář", "Říj", "Lis", "Pro",
				"Leden", "Únor", "Březen", "Duben", "Květen", "Červen", "Červenec", "Srpen", "Září", "Říjen", "Listopad", "Prosinec"
			],
			AmPm: ["do", "od", "DO", "OD"],
			S: function (j) {return j < 11 || j > 13 ? ["st", "nd", "rd", "th"][Math.min((j - 1) % 10, 3)] : "th";},
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
		locale: "cs"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		cs: $.extend({}, locInfo, { name: "čeština", nameEnglish: "Czech" }),
		"cs-CZ": $.extend({}, locInfo, { name: "čeština (Česká republika)", nameEnglish: "Czech (Czech Republic)" })
	}
});
}));
