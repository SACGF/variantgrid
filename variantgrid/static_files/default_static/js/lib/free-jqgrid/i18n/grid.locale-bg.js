/**
 * @license jqGrid Bulgarian Translation
 * Tony Tomov tony@trirand.com
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
		recordtext: "{0} - {1} от {2}",
		emptyrecords: "Няма запис(и)",
		loadtext: "Зареждам...",
		pgtext: "Стр. {0} от {1}",
		pgfirst: "Първа Стр.",
		pglast: "Последна Стр.",
		pgnext: "Следваща Стр.",
		pgprev: "Предишна Стр.",
		pgrecs: "Брой записи на Стр.",
		showhide: "Свиване/Разтягане на таблицата",
		savetext: "Записване..."
	},
	search: {
		caption: "Търсене...",
		Find: "Намери",
		Reset: "Изчисти",
		odata: [
			{ oper: "eq", text: "равно" },
			{ oper: "ne", text: "различно" },
			{ oper: "lt", text: "по-малко" },
			{ oper: "le", text: "по-малко или=" },
			{ oper: "gt", text: "по-голямо" },
			{ oper: "ge", text: "по-голямо или =" },
			{ oper: "bw", text: "започва с" },
			{ oper: "bn", text: "не започва с" },
			{ oper: "in", text: "се намира в" },
			{ oper: "ni", text: "не се намира в" },
			{ oper: "ew", text: "завършва с" },
			{ oper: "en", text: "не завършава с" },
			{ oper: "cn", text: "съдържа" },
			{ oper: "nc", text: "не съдържа" },
			{ oper: "nu", text: "е NULL" },
			{ oper: "nn", text: "не е NULL" }
		],
		groupOps: [
			{ op: "AND", text: "&nbsp;И " },
			{ op: "OR", text: "ИЛИ" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Натисни за избор на операнд.",
		resetTitle: "Изчисти стойността"
	},
	edit: {
		addCaption: "Нов Запис",
		editCaption: "Редакция Запис",
		bSubmit: "Запиши",
		bCancel: "Изход",
		bClose: "Затвори",
		saveData: "Данните са променени! Да съхраня ли промените?",
		bYes: "Да",
		bNo: "Не",
		bExit: "Отказ",
		msg: {
			required: "Полето е задължително",
			number: "Въведете валидно число!",
			minValue: "стойността трябва да е по-голяма или равна от",
			maxValue: "стойността трябва да е по-малка или равна от",
			email: "не е валиден ел. адрес",
			integer: "Въведете валидно цяло число",
			date: "Въведете валидна дата",
			url: "e невалиден URL. Изискава се префикс('http://' или 'https://')",
			nodefined: " е недефинирана!",
			novalue: " изисква връщане на стойност!",
			customarray: "Потреб. Функция трябва да върне масив!",
			customfcheck: "Потребителска функция е задължителна при този тип елемент!"
		}
	},
	view: {
		caption: "Преглед запис",
		bClose: "Затвори"
	},
	del: {
		caption: "Изтриване",
		msg: "Да изтрия ли избраният запис?",
		bSubmit: "Изтрий",
		bCancel: "Отказ"
	},
	nav: {
		edittext: "",
		edittitle: "Редакция избран запис",
		addtext: "",
		addtitle: "Добавяне нов запис",
		deltext: "",
		deltitle: "Изтриване избран запис",
		searchtext: "",
		searchtitle: "Търсене запис(и)",
		refreshtext: "",
		refreshtitle: "Обнови таблица",
		alertcap: "Предупреждение",
		alerttext: "Моля, изберете запис",
		viewtext: "",
		viewtitle: "Преглед избран запис",
		savetext: "",
		savetitle: "Съхрани запис",
		canceltext: "",
		canceltitle: "Отказ редакция"
	},
	col: {
		caption: "Избери колони",
		bSubmit: "Ок",
		bCancel: "Изход"
	},
	errors: {
		errcap: "Грешка",
		nourl: "Няма посочен url адрес",
		norecords: "Няма запис за обработка",
		model: "Модела не съответства на имената!"
	},
	formatter: {
		integer: { thousandsSeparator: " ", defaultValue: "0" },
		number: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, defaultValue: "0.00" },
		currency: { decimalSeparator: ".", thousandsSeparator: " ", decimalPlaces: 2, prefix: "", suffix: " лв.", defaultValue: "0.00" },
		date: {
			dayNames:   [
				"Нед", "Пон", "Вт", "Ср", "Чет", "Пет", "Съб",
				"Неделя", "Понеделник", "Вторник", "Сряда", "Четвъртък", "Петък", "Събота"
			],
			monthNames: [
				"Яну", "Фев", "Мар", "Апр", "Май", "Юни", "Юли", "Авг", "Сеп", "Окт", "Нов", "Дек",
				"Януари", "Февруари", "Март", "Април", "Май", "Юни", "Юли", "Август", "Септември", "Октомври", "Ноември", "Декември"
			],
			AmPm: ["", "", "", ""],
			S: function (j) {
				if (j === 7 || j === 8 || j === 27 || j === 28) {
					return "ми";
				}
				return ["ви", "ри", "ти"][Math.min((j - 1) % 10, 2)];
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
		locale: "bg"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		bg: $.extend({}, locInfo, { name: "български", nameEnglish: "Bulgarian" }),
		"bg-BG": $.extend({}, locInfo, { name: "български (България)", nameEnglish: "Bulgarian (Bulgaria)" })
	}
});
}));
