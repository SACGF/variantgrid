// Re-initialises Autocomplete light select2 widgets if loaded via Ajax load.
// @see https://github.com/yourlabs/django-autocomplete-light/issues/1221
//
// Only runs when window.__dal__initialize is already set, meaning DAL's real window.load
// event has already fired and we are in an Ajax-loaded context. Calling window.dispatchEvent
// (new Event('load')) unconditionally causes autocomplete_light.js to re-dispatch
// 'dal-init-function', which triggers a "DAL function select2 already registered" error.
// DAL 3.9+ also uses a MutationObserver to catch dynamically added elements automatically.
$(document).ready(function() {
    if (window.__dal__initialize) {
        $('[data-autocomplete-light-function]').each(window.__dal__initialize);
    }
});
