/**
 * @license jqGrid Bosnian Translation
 * by Nimesin info@nimesin.com
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
		recordtext: "Pregled {0} - {1} od {2}",
		emptyrecords: "Nema zapisa",
		loadtext: "Učitavam...",
		pgtext: "Stranica {0} od {1}",
		pgfirst: "Prva stranica",
		pglast: "Zadnja stranica",
		pgnext: "Slijedeća stranica",
		pgprev: "Prethodna stranica",
		pgrecs: "zapisa po stranici",
		showhide: "Promijeni širenje/skupljanje grida",
		savetext: "Spremanje..."
	},
	search: {
		caption: "Traži...",
		Find: "Pretraživanje",
		Reset: "Poništi",
		odata: [
			{ oper: "eq", text: "jednak" },
			{ oper: "ne", text: "nije jednak" },
			{ oper: "lt", text: "manje" },
			{ oper: "le", text: "manje ili jednako" },
			{ oper: "gt", text: "veće" },
			{ oper: "ge", text: "veće ili jednako" },
			{ oper: "bw", text: "počinje sa" },
			{ oper: "bn", text: "ne počinje sa " },
			{ oper: "in", text: "je u" },
			{ oper: "ni", text: "nije u" },
			{ oper: "ew", text: "završava sa" },
			{ oper: "en", text: "ne završava sa" },
			{ oper: "cn", text: "sadrži" },
			{ oper: "nc", text: "ne sadrži" },
			{ oper: "nu", text: "je prazno" },
			{ oper: "nn", text: "nije prazno" }
		],
		groupOps: [
			{ op: "AND", text: "sve" },
			{ op: "OR",  text: "bilo koji" }
		],
		addGroupTitle: "Dodaj podgrupu",
		deleteGroupTitle: "Obriši grupu",
		addRuleTitle: "Dodaj pravilo",
		deleteRuleTitle: "Obriši pravilo",
		operandTitle: "Pritisnite za promjenu operacije traženja.",
		resetTitle: "Vrati vrijednosti na zadane"
	},
	edit: {
		addCaption: "Dodaj zapis",
		editCaption: "Promijeni zapis",
		bSubmit: "Preuzmi",
		bCancel: "Odustani",
		bClose: "Zatvori",
		saveData: "Podaci su promijenjeni! Preuzmi promjene?",
		bYes: "Da",
		bNo: "Ne",
		bExit: "Odustani",
		msg: {
			required: "Polje je obavezno",
			number: "Molim, unesite ispravan broj",
			minValue: "Vrijednost mora biti veća ili identična ",
			maxValue: "Vrijednost mora biti manja ili identična",
			email: "neispravan e-mail",
			integer: "Molim, unijeti ispravan cijeli broj (integer)",
			date: "Molim, unijeti ispravan datum ",
			url: "neispravan URL. Prefiks je obavezan ('http://' or 'https://')",
			nodefined: " nije definiran!",
			novalue: " zahtjevan podatak je obavezan!",
			customarray: "Opcionalna funkcija trebala bi bili polje (array)!",
			customfcheck: "Kod korisničke provjere, korisnička funkcija treba biti prisutna!"

		}
	},
	view: {
		caption: "Otvori zapis",
		bClose: "Zatvori"
	},
	del: {
		caption: "Obriši",
		msg: "Obriši označen zapis ili više njih?",
		bSubmit: "Obriši",
		bCancel: "Odustani"
	},
	nav: {
		edittext: "",
		edittitle: "Promijeni obilježeni red",
		addtext: "",
		addtitle: "Dodaj novi red",
		deltext: "",
		deltitle: "Obriši obilježeni red",
		searchtext: "",
		searchtitle: "Potraži zapise",
		refreshtext: "",
		refreshtitle: "Ponovo preuzmi podatke",
		alertcap: "Upozorenje",
		alerttext: "Molim, odaberi red",
		viewtext: "",
		viewtitle: "Pregled obilježenog reda",
		savetext: "",
		savetitle: "Spremi red",
		canceltext: "",
		canceltitle: "Prekini promjenu reda"
	},
	col: {
		caption: "Obilježi kolonu",
		bSubmit: "U redu",
		bCancel: "Odustani"
	},
	errors: {
		errcap: "Greška",
		nourl: "Nedostaje URL",
		norecords: "Bez zapisa za obradu",
		model: "colNames i colModel imaju različitu duljinu!"
	},
	formatter: {
		integer: { thousandsSeparator: ".", defaultValue: "0" },
		number: { decimalSeparator: ",", thousandsSeparator: ".", decimalPlaces: 2, defaultValue: "0,00" },
		currency: { decimalSeparator: ",", thousandsSeparator: ".", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0,00" },
		date: {
			dayNames:   [
				"Ned", "Pon", "Uto", "Sri", "Čet", "Pet", "Sub",
				"Nedjelja", "Ponedjeljak", "Utorak", "Srijeda", "Četvrtak", "Petak", "Subota"
			],
			monthNames: [
				"Sij", "Velj", "Ožu", "Tra", "Svi", "Lip", "Srp", "Kol", "Ruj", "Lis", "Stu", "Pro",
				"Siječanj", "Veljača", "Ožujak", "Travanj", "Svibanj", "Lipanj", "Srpanj", "Kolovoz", "Rujan", "Listopad", "Studeni", "Prosinac"
			],
			AmPm: ["am", "pm", "AM", "PM"],
			S: function () { return ""; },
			srcformat: "Y-m-d",
			newformat: "d.m.Y.",
			masks: {
				// see http://php.net/manual/en/function.date.php for PHP format used in jqGrid
				// and see http://docs.jquery.com/UI/Datepicker/formatDate
				// and https://github.com/jquery/globalize#dates for alternative formats used frequently
				// short date:
				//    d - Day of the month, 2 digits with leading zeros
				//    m - Numeric representation of a month, with leading zeros
				//    Y - A full numeric representation of a year, 4 digits
				ShortDate: "d.m.Y.",	// in jQuery UI Datepicker: "dd.mm.yy."
				// long date:
				//    l - A full textual representation of the day of the week
				//    j - Day of the month without leading zeros
				//    F - A full textual representation of a month
				//    Y - A full numeric representation of a year, 4 digits
				LongDate: "l, j. F Y", // in jQuery UI Datepicker: "dddd, d. MMMM yyyy"
				// long date with long time:
				//    l - A full textual representation of the day of the week
				//    j - Day of the month without leading zeros
				//    F - A full textual representation of a month
				//    Y - A full numeric representation of a year, 4 digits
				//    H - 24-hour format of an hour with leading zeros
				//    i - Minutes with leading zeros
				//    s - Seconds, with leading zeros
				FullDateTime: "l, j. F Y H:i:s", // in jQuery UI Datepicker: "dddd, d. MMMM yyyy HH:mm:ss"
				// month day:
				//    d - Day of the month, 2 digits with leading zeros
				//    F - A full textual representation of a month
				MonthDay: "d F", // in jQuery UI Datepicker: "dd MMMM"
				// short time (without seconds)
				//    H - 24-hour format of an hour with leading zeros
				//    i - Minutes with leading zeros
				ShortTime: "H:i", // in jQuery UI Datepicker: "HH:mm"
				// long time (with seconds)
				//    H - 24-hour format of an hour with leading zeros
				//    i - Minutes with leading zeros
				//    s - Seconds, with leading zeros
				LongTime: "H:i:s", // in jQuery UI Datepicker: "HH:mm:ss"
				// month with year
				//    F - A full textual representation of a month
				//    Y - A full numeric representation of a year, 4 digits
				YearMonth: "F Y" // in jQuery UI Datepicker: "MMMM yyyy"
			}
		}
	}
};
$.jgrid = $.jgrid || {};
$.extend(true, $.jgrid, {
	defaults: {
		locale: "bs"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		bs: $.extend({}, locInfo, { name: "bosanski", nameEnglish: "Bosnian" }),
		"bs-BA": $.extend({}, locInfo, { name: "bosanski (Bosna)", nameEnglish: "Bosnian (Bosnia)" })
	}
});
}));
