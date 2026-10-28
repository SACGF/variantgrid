// @ts-check
// uicore/templates/uicore/page/base.html, base_external.html
/* Stands in for Rollbar where settings.ROLLBAR['enabled'] is off, so pages don't download rollbar.min.js and
   callers of Rollbar.* still work. A callback (last argument) gets an error, as an unsent report would. */
(function () {
    /** @param {...any} args */
    function notSent(...args) {
        console.log("Rollbar disabled, not sent:", ...args);
        const callback = args[args.length - 1];
        if (typeof callback === "function") {
            callback(new Error("Rollbar is disabled"));
        }
    }

    // @ts-ignore
    window.Rollbar = {
        configure: () => {},
        log: notSent,
        debug: notSent,
        info: notSent,
        warn: notSent,
        warning: notSent,
        error: notSent,
        critical: notSent,
    };
})();
