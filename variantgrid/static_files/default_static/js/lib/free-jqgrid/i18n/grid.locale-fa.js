/**
 * @license jqGrid Persian Translation
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
			recordtext: "نمابش {0} - {1} از {2}",
			emptyrecords: "رکوردی یافت نشد",
			loadtext: "بارگزاري...",
			pgtext: "صفحه {0} از {1}",
			pgfirst: "صفحه اول",
			pglast: "صفحه آخر",
			pgnext: "صفحه بعد",
			pgprev: "صفحه قبل",
			pgrecs: "رکورد در صفحه",
			showhide: "تغییر وضعیت نمایش جدول",
			savetext: "در حال ذخیره..."
		},
		search: {
			caption: "جستجو...",
			Find: "يافته ها",
			Reset: "از نو",
			odata: [
				{ oper: "eq", text: "برابر" },
				{ oper: "ne", text: "نا برابر" },
				{ oper: "lt", text: "کوچکتر" },
				{ oper: "le", text: "کوچکتر يا مساوي" },
				{ oper: "gt", text: "بزرگتر" },
				{ oper: "ge", text: "بزرگتر يا مساوي" },
				{ oper: "bw", text: "شروع با" },
				{ oper: "bn", text: "شروع نشود با" },
				{ oper: "in", text: "عضو باشد" },
				{ oper: "ni", text: "عضو این نباشد" },
				{ oper: "ew", text: "اتمام با" },
				{ oper: "en", text: "تمام نشود با" },
				{ oper: "cn", text: "حاوی" },
				{ oper: "nc", text: "نباشد حاوی" },
				{ oper: "nu", text: "تهی" },
				{ oper: "nn", text: "تهي نباشد" }
			],
			groupOps: [
				{ op: "AND", text: "کل" },
				{ op: "OR",  text: "مجموع" }
			],
			addGroupTitle: "اضافه کردن عنوان گروه",
			deleteGroupTitle: "حذف عنوان گروه",
			addRuleTitle: "اضافه کردن  عنوانقانون",
			deleteRuleTitle: "حذف عنوان قانون",
			operandTitle: "برای انتخاب گزینه های جستجو کلیک کنید",
			resetTitle: "تنظیم مجدد مقدار جستجو"
		},
		edit: {
			addCaption: "اضافه کردن رکورد",
			editCaption: "ويرايش رکورد",
			bSubmit: "ثبت",
			bCancel: "انصراف",
			bClose: "بستن",
			saveData: "دیتا تغییر کرد! ذخیره شود؟",
			bYes: "بله",
			bNo: "خیر",
			bExit: "انصراف",
			msg: {
				required: "فيلدها بايد حتما پر شوند",
				number: "لطفا عدد معتبر وارد کنيد",
				minValue: "مقدار وارد شده بايد بزرگتر يا مساوي با",
				maxValue: "مقدار وارد شده بايد کوچکتر يا مساوي",
				email: "پست الکترونيک وارد شده معتبر نيست",
				integer: "لطفا يک عدد صحيح وارد کنيد",
				date: "لطفا يک تاريخ معتبر وارد کنيد",
				url: "این آدرس صحیح نمی باشد. پیشوند نیاز است ('http://' یا 'https://')",
				nodefined: " تعریف نشده!",
				novalue: " مقدار برگشتی اجباری است!",
				customarray: "تابع شما باید مقدار آرایه داشته باشد!",
				customfcheck: "برای داشتن متد دلخواه شما باید ستون با چکینگ دلخواه داشته باشید!"
			}
		},
		view: {
			caption: "نمایش رکورد",
			bClose: "بستن"
		},
		del: {
			caption: "حذف",
			msg: "از حذف گزينه هاي انتخاب شده مطمئن هستيد؟",
			bSubmit: "حذف",
			bCancel: "ابطال"
		},
		nav: {
			edittext: "",
			edittitle: "ويرايش رديف هاي انتخاب شده",
			addtext: "",
			addtitle: "افزودن رديف جديد",
			deltext: "",
			deltitle: "حذف ردیف هاي انتخاب شده",
			searchtext: "",
			searchtitle: "جستجوي رديف",
			refreshtext: "",
			refreshtitle: "بازيابي مجدد صفحه",
			alertcap: "اخطار",
			alerttext: "لطفا يک رديف انتخاب کنيد",
			viewtext: "",
			viewtitle: "نمایش رکورد های انتخاب شده",
			savetext: "",
			savetitle: "ذخیره سطر",
			canceltext: "",
			canceltitle: "انصراف از تغییر سطر"
		},
		col: {
			caption: "نمايش/عدم نمايش ستون",
			bSubmit: "ثبت",
			bCancel: "انصراف"
		},
		errors: {
			errcap: "خطا",
			nourl: "هيچ آدرسي تنظيم نشده است",
			norecords: "هيچ رکوردي براي پردازش موجود نيست",
			model: "طول نام ستون ها مخالف ستون هاي مدل مي باشد!"
		},
		formatter: {
			integer: {
				thousandsSeparator: " ",
				defaultValue: "0"
			},
			number: {
				decimalSeparator: ".",
				thousandsSeparator: " ",
				decimalPlaces: 2,
				defaultValue: "0.00"
			},
			currency: {
				decimalSeparator: ".",
				thousandsSeparator: " ",
				decimalPlaces: 2,
				prefix: "",
				suffix: "",
				defaultValue: "0"
			},
			date: {
				dayNames: [
					"يک", "دو", "سه", "چهار", "پنج", "جمع", "شنب",
					"یکشنبه", "دوشنبه", "سه‌شنبه", "چهارشنبه", "پنجشنبه", "جمعه", "شنبه"
				],
				monthNames: [
					"ژانویه", "فوریه", "مارس", "آوریل", "مه", "ژوئن", "ژوئیه", "اوت", "سپتامبر", "اکتبر", "نوامبر", "دسامبر",
					"ژانویه", "فوریه", "مارس", "آوریل", "مه", "ژوئن", "ژوئیه", "اوت", "سپتامبر", "اکتبر", "نوامبر", "دسامبر"
				],
				AmPm: ["ب.ظ", "ب.ظ", "ق.ظ", "ق.ظ"],
				S: function (b) {
					return b < 11 || b > 13 ? ["st", "nd", "rd", "th"][Math.min((b - 1) % 10, 3)] : "th";
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
			locale: "fa"
		},
		locales: {
			// In general the property name is free, but it's recommended to use the names based on
			// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
			// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
			// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
			// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
			fa: $.extend({}, locInfo, { name: "فارسى", nameEnglish: "Persian" }),
			"fa-IR": $.extend({}, locInfo, { name: "فارسى (ایران)", nameEnglish: "Persian (Iran)" })
		}
	});
}));
