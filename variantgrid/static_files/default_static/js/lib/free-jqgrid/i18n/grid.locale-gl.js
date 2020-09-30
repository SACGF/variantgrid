/**
 * @license jqGrid Galician Translation
 * Translated by Jorge Barreiro <yortx.barry@gmail.com>
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
		recordtext: "Amosando {0} - {1} de {2}",
		emptyrecords: "Sen rexistros que amosar",
		loadtext: "Cargando...",
		pgtext: "Páxina {0} de {1}",
		pgfirst: "First Page",
		pglast: "Last Page",
		pgnext: "Next Page",
		pgprev: "Previous Page",
		pgrecs: "Records per Page",
		showhide: "Toggle Expand Collapse Grid",
		savetext: "Gardando..."
	},
	search: {
		caption: "Búsqueda...",
		Find: "Buscar",
		Reset: "Limpar",
		odata: [
			{ oper: "eq", text: "igual " },
			{ oper: "ne", text: "diferente a" },
			{ oper: "lt", text: "menor que" },
			{ oper: "le", text: "menor ou igual que" },
			{ oper: "gt", text: "maior que" },
			{ oper: "ge", text: "maior ou igual a" },
			{ oper: "bw", text: "empece por" },
			{ oper: "bn", text: "non empece por" },
			{ oper: "in", text: "está en" },
			{ oper: "ni", text: "non está en" },
			{ oper: "ew", text: "termina por" },
			{ oper: "en", text: "non termina por" },
			{ oper: "cn", text: "contén" },
			{ oper: "nc", text: "non contén" },
			{ oper: "nu", text: "is null" },
			{ oper: "nn", text: "is not null" }
		],
		groupOps: [
			{ op: "AND", text: "todo" },
			{ op: "OR",  text: "calquera" }
		],
		addGroupTitle: "Add subgroup",
		deleteGroupTitle: "Delete group",
		addRuleTitle: "Add rule",
		deleteRuleTitle: "Delete rule",
		operandTitle: "Click to select search operation.",
		resetTitle: "Reset Search Value"
	},
	edit: {
		addCaption: "Engadir rexistro",
		editCaption: "Modificar rexistro",
		bSubmit: "Gardar",
		bCancel: "Cancelar",
		bClose: "Pechar",
		saveData: "Modificáronse os datos, quere gardar os cambios?",
		bYes: "Si",
		bNo: "Non",
		bExit: "Cancelar",
		msg: {
			required: "Campo obrigatorio",
			number: "Introduza un número",
			minValue: "O valor debe ser maior ou igual a ",
			maxValue: "O valor debe ser menor ou igual a ",
			email: "non é un enderezo de correo válido",
			integer: "Introduza un valor enteiro",
			date: "Introduza unha data correcta ",
			url: "non é unha URL válida. Prefixo requerido ('http://' ou 'https://')",
			nodefined: " non está definido.",
			novalue: " o valor de retorno é obrigatorio.",
			customarray: "A función persoalizada debe devolver un array.",
			customfcheck: "A función persoalizada debe estar presente no caso de ter validación persoalizada."
		}
	},
	view: {
		caption: "Consultar rexistro",
		bClose: "Pechar"
	},
	del: {
		caption: "Eliminar",
		msg: "Desexa eliminar os rexistros seleccionados?",
		bSubmit: "Eliminar",
		bCancel: "Cancelar"
	},
	nav: {
		edittext: "",
		edittitle: "Modificar a fila seleccionada",
		addtext: "",
		addtitle: "Engadir unha nova fila",
		deltext: "",
		deltitle: "Eliminar a fila seleccionada",
		searchtext: "",
		searchtitle: "Buscar información",
		refreshtext: "",
		refreshtitle: "Recargar datos",
		alertcap: "Aviso",
		alerttext: "Seleccione unha fila",
		viewtext: "",
		viewtitle: "Ver fila seleccionada",
		savetext: "",
		savetitle: "Save row",
		canceltext: "",
		canceltitle: "Cancel row editing"
	},
	col: {
		caption: "Mostrar/ocultar columnas",
		bSubmit: "Enviar",
		bCancel: "Cancelar"
	},
	errors: {
		errcap: "Erro",
		nourl: "Non especificou unha URL",
		norecords: "Non hai datos para procesar",
		model: "As columnas de nomes son diferentes das columnas de modelo"
	},
	formatter: {
		integer: { thousandsSeparator: ".", defaultValue: "0" },
		number: { decimalSeparator: ",", thousandsSeparator: ".", decimalPlaces: 2, defaultValue: "0,00" },
		currency: { decimalSeparator: ",", thousandsSeparator: ".", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0,00" },
		date: {
			dayNames:   [
				"Do", "Lu", "Ma", "Me", "Xo", "Ve", "Sa",
				"Domingo", "Luns", "Martes", "Mércoles", "Xoves", "Vernes", "Sábado"
			],
			monthNames: [
				"Xan", "Feb", "Mar", "Abr", "Mai", "Xuñ", "Xul", "Ago", "Set", "Out", "Nov", "Dec",
				"Xaneiro", "Febreiro", "Marzo", "Abril", "Maio", "Xuño", "Xullo", "Agosto", "Setembro", "Outubro", "Novembro", "Decembro"
			],
			AmPm: ["am", "pm", "AM", "PM"],
			S: function (j) {
				return j < 11 || j > 13 ? ["st", "nd", "rd", "th"][Math.min((j - 1) % 10, 3)] : "th";
			},
			srcformat: "Y-m-d",
			newformat: "d-m-Y",
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
		locale: "gl"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		gl: $.extend({}, locInfo, { name: "galego", nameEnglish: "Galician" }),
		"gl-ES": $.extend({}, locInfo, { name: "galego (galego)", nameEnglish: "Galician (Galician)" })
	}
});
}));
