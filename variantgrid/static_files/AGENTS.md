# static_files — agent notes

Source JS/CSS/images live here, under `<site>_static/` (`default_static` unless site specific) - always edit here.
`variantgrid/sitestatic/` is collectstatic output: gitignored and overwritten.

`global.css` and friends are compiled from `.scss` by a PyCharm file watcher - do not run `sassc`/`sass` yourself (its
formatting creates huge diffs). Edit the `.scss`, then hand-apply the same minimal change to the generated `.css`
matching its formatting, so it works before the next recompile. Leave `.css.map` files alone.

`scripts/vg css unused` lists the class / id selectors in the `.scss` that nothing in the templates, JS or Python
names, so run it before deleting or moving a rule (`--dynamic` adds the names built up at runtime, `cs-{{ status }}`
style, which need a reader's eye). It reads `library/vg/css.py`.

Bootstrap 4: `data-toggle` / `data-target`, not the Bootstrap 5 `data-bs-*` forms.
