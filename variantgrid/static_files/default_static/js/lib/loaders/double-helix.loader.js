/*
 * Double Helix loader - wraps the existing jQuery DoubleHelix plugin
 * (js/lib/double-helix.jquery.js) in the VGLoaders interface so the current
 * animation sits alongside the new ones for comparison.
 *
 * Requires: jQuery + double-helix.jquery.js loaded first.
 */
(function () {
    "use strict";

    VGLoaders.register("double-helix", {
        label: "Double Helix (current)",
        start: function (container, opts) {
            const fgColor = opts && opts.theme === "light" ? "40,48,58" : "210,228,232";
            const $container = jQuery(container);
            const canvas = jQuery("<canvas />", { class: "node-load-animation" });
            canvas.attr({ width: 50, height: 180 });
            canvas.css("opacity", 0.35);
            canvas.DoubleHelix({ fps: 20, spinSpeed: 4, fgColor: fgColor });
            $container.append(canvas);

            return function stop() {
                canvas.each(function () { this.active = false; });
                canvas.remove();
            };
        }
    });
})();
