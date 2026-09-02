# #1458 — Switch out UniSA for University of Adelaide

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-02

## Background

The Centre for Cancer Biology (CCB) moved from a UniSA alliance to a University of Adelaide one.
UniSA branding needs replacing with University of Adelaide branding, and the RUNX1db contact
email moved from `CCB-RUNX1DB@unisa.edu.au` to `CCB-RUNX1DB@adelaide.edu.au`.

## Already done (this repo, commits 424df8ef0 and 7717a32ea, March 2026)

- `variantgrid/templates/runx1_templates/map_matrix.html` — contact email now `CCB-RUNX1DB@adelaide.edu.au`.
- `variantgrid/static_files/runx1_static/images/unisa-1.png` deleted; `au-logo-dark-blue-horizontal.svg`
  (University of Adelaide logo, 800×191 viewBox) added.
- `variantgrid/static_files/runx1_static/help/about.html` — `.unisa-logo` now points at the Adelaide SVG.

## Remaining work

### 1. `variantgrid_sapath` repo (the item in the issue's last comment)

Sibling checkout at `../variantgrid_sapath`. Only one place references UniSA:

- `variantgrid/templates/sapathology_templates/registration/login.html`
  - lines 90–95: `.unisa-logo` CSS rule (193×150, `unisa-1.png`)
  - line 247: `<div class='logo unisa-logo'>` in `#bottom-logos`
- `variantgrid/static_files/sapathology_static/images/unisa-1.png` — the old logo.

Steps:
1. Copy `au-logo-dark-blue-horizontal.svg` from this repo's `runx1_static/images/` into
   `variantgrid_sapath/.../sapathology_static/images/`.
2. Delete `sapathology_static/images/unisa-1.png`.
3. In `login.html`, rename `.unisa-logo` → `.adelaide-logo` (CSS rule and the div), point it at the SVG,
   and size it to match the row. The other bottom logos are ~193–200px wide and 150–195px tall in a
   flex row; the Adelaide logo is wide (800:191, ≈4.2:1), so use something like `width: 320px; height: 76px;
   background-size: contain; background-repeat: no-repeat` rather than forcing 193×150 (which would squash it).
4. Check `sapathology_testvm_static`, `sapathology_training_static` and `sapathology_vg3upgrade_static` —
   they hold no `images/` dirs or UniSA references today, so nothing to do unless the login page is
   overridden there (it is not; only one `login.html` exists).
5. Run `collectstatic` on a sapath settings module and eyeball the login page.

### 2. RUNX1 about page — tidy-up in this repo

`variantgrid/static_files/runx1_static/help/about.html`:

- Line 51 text still says CCB is "an alliance between SA Pathology and the **University of South Australia**"
  — change to **University of Adelaide**.
- The `.unisa-logo` rule was retargeted but left at 628×150 next to 193px and 320px logos in an 800px-wide
  `#logos-container` (193 + 320 + 628 > 800, so the floats wrap). Rename the class to `.adelaide-logo` and
  shrink it to a height that matches the row (e.g. `width: 300px; height: 72px; background-size: contain`),
  or drop the fixed container width.

### 3. Verify nothing else

`grep -rniI unisa` across both repos (excluding `.venv`, `sitestatic`, `.git`) should return nothing once
the above is done. Also grep `university of south australia` — the only hit today is about.html line 51.

## Out of scope

- `paper/2018_paper.txt` affiliations are historical and stay as-is.
- `settings/env/vgaws.py` names "Centre for Cancer Biology" only — no UniSA reference.

## Testing

No unit tests apply; this is static content. Manual check: render the SA Path login page and the RUNX1
`/help/about.html` page after `collectstatic` and confirm the logo row lays out on one line.
