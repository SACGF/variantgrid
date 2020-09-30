/**
 * @license jqGrid Polish Translation
 * Łukasz Schab lukasz@freetree.pl
 * http://FreeTree.pl
 *
 * Updated names, abbreviations, currency and date/time formats for Polish norms (also corresponding with CLDR v21.0.1 --> http://cldr.unicode.org/index)
 * Tomasz Pęczek tpeczek@gmail.com
 * http://tpeczek.blogspot.com; http://tpeczek.codeplex.com
 *
 * 22-01-2015: Updated locale strings
 * Robert 'Wooya' Gaudyn <wogaew@wp.pl>
 *
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
		recordtext: "Pokaż {0} - {1} z {2}",
		emptyrecords: "Brak rekordów do pokazania",
		loadtext: "Ładowanie...",
		pgtext: "Strona {0} z {1}",
		pgfirst: "Pierwsza strona",
		pglast: "Ostatnia strona",
		pgnext: "Następna strona",
		pgprev: "Poprzednia strona",
		pgrecs: "Rekordów na stronę",
		showhide: "Przełącz zwinięcie/rozwinięcie tabeli",
		savetext: "Trwa zapisywanie..."
	},
	search: {
		caption: "Wyszukiwanie...",
		Find: "Szukaj",
		Reset: "Czyść",
		odata: [
			{ oper: "eq", text: "dokładnie" },
			{ oper: "ne", text: "różne od" },
			{ oper: "lt", text: "mniejsze od" },
			{ oper: "le", text: "mniejsze lub równe" },
			{ oper: "gt", text: "większe od" },
			{ oper: "ge", text: "większe lub równe" },
			{ oper: "bw", text: "zaczyna się od" },
			{ oper: "bn", text: "nie zaczyna się od" },
			{ oper: "in", text: "jest w" },
			{ oper: "ni", text: "nie jest w" },
			{ oper: "ew", text: "kończy się na" },
			{ oper: "en", text: "nie kończy się na" },
			{ oper: "cn", text: "zawiera" },
			{ oper: "nc", text: "nie zawiera" },
			{ oper: "nu", text: "jest null" },
			{ oper: "nn", text: "nie jest null" }
		],
		groupOps: [
			{ op: "AND", text: "oraz" },
			{ op: "OR",  text: "lub" }
		],
		addGroupTitle: "Dodaj podgrupę",
		deleteGroupTitle: "Usuń grupę",
		addRuleTitle: "Dodaj regułę",
		deleteRuleTitle: "Usuń regułę",
		operandTitle: "Kliknij, aby wybrać operację wyszukiwania.",
		resetTitle: "Zresetuj wartość wyszukiwania"
	},
	edit: {
		addCaption: "Dodaj rekord",
		editCaption: "Edytuj rekord",
		bSubmit: "Zapisz",
		bCancel: "Anuluj",
		bClose: "Zamknij",
		saveData: "Dane zostały zmienione! Zapisać zmiany?",
		bYes: "Tak",
		bNo: "Nie",
		bExit: "Anuluj",
		msg: {
			required: "Pole jest wymagane",
			number: "Proszę wpisać poprawną liczbę",
			minValue: "Wartość musi być większa lub równa od",
			maxValue: "Wartość musi być mniejsza lub równa od",
			email: "nie jest poprawnym adresem e-mail",
			integer: "Proszę wpisać poprawną liczbę",
			date: "Proszę podaj poprawną datę",
			url: "jest niewłaściwym adresem URL. Pamiętaj o prefiksie ('http://' lub 'https://')",
			nodefined: " niezdefiniowane!",
			novalue: " wymagana jest wartość zwracana!",
			customarray: "Funkcja niestandardowa powinna zwracać tablicę!",
			customfcheck: "Funkcja niestandardowa powinna być obecna w przypadku niestandardowego sprawdzania!"
		}
	},
	view: {
		caption: "Pokaż rekord",
		bClose: "Zamknij"
	},
	del: {
		caption: "Usuń",
		msg: "Czy usunąć wybrany rekord(y)?",
		bSubmit: "Usuń",
		bCancel: "Anuluj"
	},
	nav: {
		edittext: "",
		edittitle: "Edytuj wybrany wiersz",
		addtext: "",
		addtitle: "Dodaj nowy wiersz",
		deltext: "",
		deltitle: "Usuń wybrany wiersz",
		searchtext: "",
		searchtitle: "Wyszukaj rekord",
		refreshtext: "",
		refreshtitle: "Przeładuj",
		alertcap: "Uwaga",
		alerttext: "Proszę wybrać wiersz",
		viewtext: "",
		viewtitle: "Pokaż wybrany wiersz",
		savetext: "",
		savetitle: "Zapisz wiersz",
		canceltext: "",
		canceltitle: "Anuluj edycję wiersza"
	},
	col: {
		caption: "Pokaż/Ukryj kolumny",
		bSubmit: "Zatwierdź",
		bCancel: "Anuluj"
	},
	errors: {
		errcap: "Błąd",
		nourl: "Brak adresu url",
		norecords: "Brak danych",
		model: "Długość colNames <> colModel!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ",", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0,00" },
		currency: { decimalSeparator: ",", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: " zł", defaultValue: "0,00" },
		date: {
			dayNames:   [
				"niedz.", "pon.", "wt.", "śr.", "czw.", "pt.", "sob.",
				"niedziela", "poniedziałek", "wtorek", "środa", "czwartek", "piątek", "sobota"
			],
			monthNames: [
				"sty", "lut", "mar", "kwi", "maj", "cze", "lip", "sie", "wrz", "paź", "lis", "gru",
				"styczeń", "luty", "marzec", "kwiecień", "maj", "czerwiec", "lipiec", "sierpień", "wrzesień", "październik", "listopad", "grudzień"
				],
			AmPm: ["", "", "", ""],
			S: function () { return ""; },
			srcformat: "Y-m-d",
			newformat: "d.m.Y",
			masks: {
				ShortDate: "d.m.y",
				LongDate: "l, j F Y",
				FullDateTime: "l, j F Y H:i:s",
				MonthDay: "j F",
				ShortTime: "H:i",
				LongTime: "H:i:s",
				YearMonth: "F Y"
			}
		}
	}
};
$.jgrid = $.jgrid || {};
$.extend(true, $.jgrid, {
	defaults: {
		locale: "pl"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		pl: $.extend({}, locInfo, { name: "polski", nameEnglish: "Polish" }),
		"pl-PL": $.extend({}, locInfo, { name: "polski (Polska)", nameEnglish: "Polish (Poland)" })
	}
});
}));
