# static_files — agent notes

Source JS/CSS/images live here, under `<site>_static/` (`default_static` unless site specific) - always edit here.
`variantgrid/sitestatic/` is collectstatic output: gitignored and overwritten.

`global.css` and friends are compiled from `.scss` by a PyCharm file watcher - do not run `sassc`/`sass` yourself (its
formatting creates huge diffs). Edit the `.scss`, then hand-apply the same minimal change to the generated `.css`
matching its formatting, so it works before the next recompile. Leave `.css.map` files alone. Compiling to a scratch
directory to check a change is fine: `sassc -t expanded` output matches the committed `.css` bar blank lines
(`diff -B -w`), but it cannot read `@use`, so a `global.scss` change needs a dart-sass binary.

`global.scss` is only a list of `@use` lines; the rules are in `css/global/_<concern>.scss`, and the `@use` order is the
cascade. Moving a rule between partials moves it in the cascade, so compile before and after and diff. Nothing styles
`<fieldset>`: a group is a Bootstrap `.card` (the legend as its `card-header`), forms are
`{% crispy form form_helper.horizontal %}` and label/value rows `{% labelled %}`. A whole-form `{{ form }}` wraps each
multiple-choice field in a bare `<fieldset>`, so render those fields one at a time; tables are Bootstrap `.table`.

`scripts/vg css unused` lists the class / id selectors in the `.scss` that nothing in the templates, JS or Python
names, so run it before deleting or moving a rule (`--dynamic` adds the names built up at runtime, `cs-{{ status }}`
style, which need a reader's eye). The edit hook (`.claude/hooks/post_edit.py`) runs it after every `.html`, `.scss`
or first-party `.js` edit and warns when anything is unused; it is deliberately not a CI check, since a class built
from data or added by JS can be live while nothing in the tree names it. It reads `library/vg/css.py`. `vg css unused --rendered` (`library/vg/css_rendered.py`)
crawls pages as `claude_agent` and reports which dynamic names reach rendered HTML - it can't see classes JS adds, so
most grid and flag names still need reading the code that builds them. A class whose name is a database value
(`SequencerModel.css_class`) goes in `DATA_DERIVED` there, or the scan calls it unused.

Bootstrap 4: `data-toggle` / `data-target`, not the Bootstrap 5 `data-bs-*` forms.

Third-party libraries live under `default_static/js/lib/` as our own copies, never loaded from a CDN (#1741), and the
version goes in the path (`variantgrid/static_files/default_static/js/lib/jquery-3.7.1/jquery-3.7.1.min.js`,
`variantgrid/static_files/default_static/js/lib/moment-2.30.1/moment-2.30.1.min.js`) so an upgrade
is a rename that shows in every template, and a scanner can read the version without opening the file (#1955).
Copy the file the project publishes (npm `min/` for moment) rather than a CDN's re-minified build, so the
`.map` matches. `double-helix.jquery.js` is the one unversioned file: a loading-animation snippet with no upstream
release.

Page JS goes in a static file by default, written there from the start rather than inline and moved later: static JS
is cached by the browser, linted, and can't change with the context it renders in. It lives in
`default_static/js/<app>/<template_name>.js` (#2021). The
static file knows nothing about the page's context: it defines functions whose parameters are the page values (one
`options` object with named keys once there are several, an `init<Name>(options)` entry point for a file called once per
page). The template keeps a short inline `<script>` after the `<script src>` that calls it with
`{{ value|jsonify }}` (`uicore/templatetags/js_tags.py`, a JS literal, so NaN passes through), e.g. `analysis/templates/analysis/node_data/node_data_grid.html`.
An inclusion tag shown several times on a page, or a fragment loaded by ajax, makes the same call from its own block.
The only JS left inline is the short init call and blocks that are mostly template logic (`{% if %}` around a line or
two) or callbacks that are nothing but page values (`classification/templates/classification/classification_dashboard.html`). URLs come from
`Urls.<name>(...)` (django_js_reverse). New files start with `// @ts-check` and must pass `npx eslint`. A new static
file 500s every page that references it until `collectstatic` has added it to the manifest
(`ManifestStaticFilesStorage`), so it ships via `scripts/upgrade.sh`, never a bare checkout on a served tree.
