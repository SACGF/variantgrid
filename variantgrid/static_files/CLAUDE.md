# static_files — agent notes

Source JS/CSS/images live here, under `<site>_static/` (`default_static` unless site specific) - always edit here.
`variantgrid/sitestatic/` is collectstatic output: gitignored and overwritten.

`global.css` and friends are compiled from `.scss` by a PyCharm file watcher - do not run `sassc`/`sass` yourself (its
formatting creates huge diffs). Edit the `.scss`, then hand-apply the same minimal change to the generated `.css`
matching its formatting, so it works before the next recompile. Leave `.css.map` files alone.

Bootstrap 4: `data-toggle` / `data-target`, not the Bootstrap 5 `data-bs-*` forms.
