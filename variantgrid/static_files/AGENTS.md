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
