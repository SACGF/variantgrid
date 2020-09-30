/**
 * @license jqGrid Spanish Translation
 * Traduccion jqGrid en Español por Yamil Bracho
 *
 * Traduccion corregida y ampliada por Faserline, S.L.
 * http://www.faserline.com
 *
 * Traducción corregida y ampliada por Marc lobato
 * http://www.navigatecms.com
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
		recordtext: "Mostrando {0} - {1} de {2}",
		emptyrecords: "Sin registros que mostrar",
		loadtext: "Cargando...",
		pgtext: "Página {0} de {1}",
		pgfirst: "Primera página",
		pglast: "Última página",
		pgnext: "Página siguiente",
		pgprev: "Página anterior",
		pgrecs: "Registros por página",
		showhide: "Alternar expandir/colapsar la tabla",
		savetext: "Guardando..."
	},
	search: {
		caption: "Búsqueda...",
		Find: "Buscar",
		Reset: "Limpiar",
		odata: [
			{ oper: "eq", text: "igual" },
			{ oper: "ne", text: "no igual a" },
			{ oper: "lt", text: "menor que" },
			{ oper: "le", text: "menor o igual que" },
			{ oper: "gt", text: "mayor que" },
			{ oper: "ge", text: "mayor o igual a" },
			{ oper: "bw", text: "empieza por" },
			{ oper: "bn", text: "no empieza por" },
			{ oper: "in", text: "está en" },
			{ oper: "ni", text: "no está en" },
			{ oper: "ew", text: "termina por" },
			{ oper: "en", text: "no termina por" },
			{ oper: "cn", text: "contiene" },
			{ oper: "nc", text: "no contiene" },
			{ oper: "nu", text: "es nulo" },
			{ oper: "nn", text: "no es nulo" }
		],
		groupOps: [
			{ op: "AND", text: "todo" },
			{ op: "OR",  text: "cualquier" }
		],
		addGroupTitle: "Crear subgrupo",
		deleteGroupTitle: "Eliminar grupo",
		addRuleTitle: "Crear regla",
		deleteRuleTitle: "Eliminar regla",
		operandTitle: "Clic para seleccionar la operación de búsqueda.",
		resetTitle: "Reiniciar valor de búsqueda"
	},
	edit: {
		addCaption: "Agregar registro",
		editCaption: "Modificar registro",
		bSubmit: "Guardar",
		bCancel: "Cancelar",
		bClose: "Cerrar",
		saveData: "Se han modificado los datos, ¿guardar cambios?",
		bYes: "Si",
		bNo: "No",
		bExit: "Cancelar",
		msg: {
			required: "Campo obligatorio",
			number: "Introduzca un número",
			minValue: "El valor debe ser mayor o igual a ",
			maxValue: "El valor debe ser menor o igual a ",
			email: "no es una dirección de correo válida",
			integer: "Introduzca un valor entero",
			date: "Introduza una fecha correcta ",
			url: "no es una URL válida. Prefijo requerido ('http://' or 'https://')",
			nodefined: " no está definido.",
			novalue: " valor de retorno es requerido.",
			customarray: "La función personalizada debe devolver un array.",
			customfcheck: "La función personalizada debe estar presente en el caso de validación personalizada."
		}
	},
	view: {
		caption: "Consultar registro",
		bClose: "Cerrar"
	},
	del: {
		caption: "Eliminar",
		msg: "¿Desea eliminar los registros seleccionados?",
		bSubmit: "Eliminar",
		bCancel: "Cancelar"
	},
	nav: {
		edittext: "",
		edittitle: "Modificar fila seleccionada",
		addtext: "",
		addtitle: "Agregar nueva fila",
		deltext: "",
		deltitle: "Eliminar fila seleccionada",
		searchtext: "",
		searchtitle: "Buscar información",
		refreshtext: "",
		refreshtitle: "Recargar datos",
		alertcap: "Aviso",
		alerttext: "Seleccione una fila",
		viewtext: "",
		viewtitle: "Ver fila seleccionada",
		savetext: "",
		savetitle: "Guardar fila",
		canceltext: "",
		canceltitle: "Cancelar edición de fila"
	},
	col: {
		caption: "Mostrar/ocultar columnas",
		bSubmit: "Enviar",
		bCancel: "Cancelar"
	},
	errors: {
		errcap: "Error",
		nourl: "No se ha especificado ninguna URL",
		norecords: "No hay datos para procesar",
		model: "Las columnas de nombres son diferentes de las columnas de modelo"
	},
	formatter: {
		integer: { thousandsSeparator: ".", defaultValue: "0" },
		number: { decimalSeparator: ",", thousandsSeparator: ".", decimalPlaces: 2, defaultValue: "0,00" },
		currency: { decimalSeparator: ",", thousandsSeparator: ".", decimalPlaces: 2, prefix: "", suffix: "", defaultValue: "0,00" },
		date: {
			dayNames:   [
				"Do", "Lu", "Ma", "Mi", "Ju", "Vi", "Sa",
				"Domingo", "Lunes", "Martes", "Miercoles", "Jueves", "Viernes", "Sabado"
			],
			monthNames: [
				"Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic",
				"Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"
			],
			AmPm: ["am", "pm", "AM", "PM"],
			S: function (j) {return j < 11 || j > 13 ? ["st", "nd", "rd", "th"][Math.min((j - 1) % 10, 3)] : "th";},
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
		locale: "es"
	},
	locales: {
		// In general the property name is free, but it's recommended to use the names based on
		// http://www.iana.org/assignments/language-subtag-registry/language-subtag-registry
		// http://rishida.net/utils/subtags/ and RFC 5646. See Appendix A of RFC 5646 for examples.
		// One can use the lang attribute to specify language tags in HTML, and the xml:lang attribute for XML
		// if it exists. See http://www.w3.org/International/articles/language-tags/#extlang
		es: $.extend({}, locInfo, { name: "español", nameEnglish: "Spanish" }),
		"es-ES": $.extend({}, locInfo, { name: "Español (España)", nameEnglish: "Spanish (Spain)" })
	}
});
}));
