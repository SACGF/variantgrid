/*
 * VGLoaders - a tiny registry of "loading animation" implementations.
 *
 * Each loader registers itself with a uniform interface so they are
 * interchangeable. The eventual goal is for the analysis node-loading overlay
 * (see showLoadingOverlay() in analysis.js) to pick one - either the user's
 * chosen animation, or a random one - and run it the same way regardless of
 * which animation it is.
 *
 * Loader interface:
 *   VGLoaders.register(id, {
 *     label: "Human readable name",
 *     start: function(container) {
 *        // container: a DOM element (or jQuery object) to draw inside.
 *        // The loader builds whatever it needs (usually a <canvas>) inside it.
 *        // Returns a stop() function that halts animation + cleans up.
 *        return function stop() { ... };
 *     }
 *   });
 */
(function (global) {
    "use strict";

    const registry = {};
    const order = [];

    function toEl(container) {
        // Accept a raw DOM node or a jQuery object.
        if (container && container.jquery) {
            return container[0];
        }
        return container;
    }

    const VGLoaders = {
        register: function (id, loader) {
            if (registry[id]) {
                console.warn("VGLoaders: overwriting loader '" + id + "'");
            } else {
                order.push(id);
            }
            registry[id] = Object.assign({ id: id }, loader);
            return VGLoaders;
        },

        get: function (id) {
            return registry[id];
        },

        ids: function () {
            return order.slice();
        },

        list: function () {
            return order.map(id => registry[id]);
        },

        random: function () {
            const id = order[Math.floor(Math.random() * order.length)];
            return registry[id];
        },

        /*
         * Start a loader by id inside container. Returns a stop() function.
         * This is the single entry point the VG overlay would call.
         * opts is passed through to the loader, e.g. { theme: "light" | "dark" }.
         */
        start: function (id, container, opts) {
            const loader = registry[id];
            if (!loader) {
                console.error("VGLoaders: no loader '" + id + "'");
                return function () {};
            }
            return loader.start(toEl(container), opts || {}) || function () {};
        }
    };

    global.VGLoaders = VGLoaders;
})(window);
