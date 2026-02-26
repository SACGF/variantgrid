# VariantGrid UICore App — Reference Document

## Purpose and Overview

The `uicore` app provides reusable UI utilities, template tags, form helpers, custom widgets, AJAX views, and JSON handling for the entire VariantGrid platform. It abstracts common UI patterns to reduce duplication.

**Key Design:** Utility-only — no models defined. Heavy template tag library. Custom form widgets with "other" option support. AJAX form rendering with multiple modes.

---

## Views

### AjaxFormView (Generic[T])
Base class for AJAX-enabled form views with multiple render modes.
- Method: `lazy_render(obj: T, context=None)` → LazyRender

### LazyRender (Generic[T], dataclass)
Deferred rendering wrapper for template content.
- Attributes: core_object, core_object_name, template_name, static_context, dynamic_context (callable)
- Methods:
  - `embed(request, mode=EMBEDDED_CARD, **kwargs)` → SafeString — returns safe HTML for template embedding
  - `render(request, saved=False)` → HttpResponse

**Render Modes:**
- **INLINE**: Inline div with ID for replacement on save
- **MODAL**: Opens in Bootstrap modal, closes on save
- **CARD**: Standard card div wrapper
- **EMBEDDED_CARD**: Card wrapper replaceable on edit/save

---

## Form Helpers

### FormHelperHelper
Builder for Crispy Forms FormHelper instances.
- Methods: `no_csrf()`, properties: `horizontal`, `horizontal_nested`, `horizontal_no_labels`, `fields_only`
- `form_helper_horizontal()` — convenience function returning horizontal FormHelper

---

## Custom Widgets

### RadioOtherWidget / CheckboxOtherWidget
Custom radio/checkbox widgets with "other" option.
- Templates: `uicore/widgets/radio_other_widget.html`, `uicore/widgets/checkbox_other_widget.html`
- **ValuesMissingOther** — wrapper for "other" selected but no value provided (triggers validation error)
- **OptionData** (frozen dataclass) — represents a single option: name, label, value, selected, index

### ChoiceFieldWithOther / MultiChoiceFieldWithOther
Form fields using the above widgets.
- Validation: raises ValidationError if ValuesMissingOther
- valid_value(): always returns True; to_python(): returns value as-is

---

## JSON Utilities (validated_json.py)

### JsonMessage (frozen dataclass)
- Fields: severity ("error", "warning", "info"), text
- Properties: is_error, bs (Bootstrap class)
- Methods: to_json(), deserialize()

### JsonMessages (frozen dataclass)
- Methods: errors(), warnings(), error()/warning()/info() (static constructors), `__add__()`, `__bool__()`, `__iter__()`

### ValidatedJson
Wraps JSON data with associated validation messages. Used for ClinVar exports and other serialization that needs validation results alongside data.

---

## Template Tags

### ui_utils.py — Filters
| Filter | Purpose |
|--------|---------|
| `jsonify(json_me)` | Python → JS-safe JSON |
| `jsonify_pretty(json_me)` | Pretty-printed JSON |
| `jsstring(text)` | Escape for JS strings |
| `js_symbol(value)` | Python value → valid JS identifier |
| `format_value(val, limit=0)` | Type-aware HTML formatting (None, empty str, dates, dicts, floats, ints) |
| `format_computer_text(val)` | Replaces & and _ with spaces |
| `format_unit_as_percent(val)` | Multiplies by 100, adds % |
| `dash_if_empty(val)` | Shows `-` for empty/None |
| `severity_icon(severity, title)` | FontAwesome icon for severity levels (C/E/W/I/D/S/CREATE/PROCESS/etc.) |
| `severity_bs(severity)` | Bootstrap class for LogLevel |
| `danger_badge(count)` | Red badge if >0, green if 0 |
| `boolean(test)` | Checkmark or X icon |
| `secret(value, length=-4)` | Masks with **** |
| `duration(td)` | Formats timedelta as "1d 2h 3m 4s" |
| `format_preference(value)` | Boolean/TextChoices display |
| `get_item(dictionary, key)` | Dict access workaround for templates |
| `enrich(text)` | Double-quoted text → highlighted spans |
| `emojify(text)` | Emoji codes → unicode |
| `bool_emojify(val)` | True → ✅, False → ❌ |
| `pretty_label(code)` | code_style → Pretty Label |

### ui_utils.py — Inclusion Tags
| Tag | Purpose |
|-----|---------|
| `update_django_messages(context)` | JS to update page messages from AJAX context |
| `code_json(data, css_class, dash_if_empty)` | JSON in code block |
| `code_xml(data, css_class)` | XML code block |
| `code_regex(data)` | Regex display with validation |
| `code_shell(data)` | Shell command display |
| `timedelta(time, show_micro)` | Timedelta breakdown display |
| `timestamp(timestamp, ...)` | Flexible timestamp display (Unix/datetime/date; time-ago mode; micro/second precision) |
| `diff_text(a, b)` | Side-by-side text diff |
| `admin_link(context, object)` | Django admin edit link (superuser only) |
| `value_with_icon(value, help, icon)` | Value with help popover and icon |
| `preview(obj)` | Renders PreviewModelMixin preview |

### ui_utils.py — Block Tags
| Tag | Purpose |
|-----|---------|
| `install-instructions(label, installed)` | Collapsible install instructions (superuser only) |
| `field_help` | Form-text help styling wrapper |
| `labelled(label, hint, ...)` | Labeled value row (hints: "tiny", "chunky", "inline", "large-label"); help popovers; admin-only; zero values styled distinctly |
| `modal(id, label, size, admin_only, button)` | Bootstrap modal with toggle link (sizes: xl/lg/md/sm) |
| `if_user_can_edit(object_expression)` | Renders content if user can write to object |

### english_tags.py
- `count(items, singular, plural=None)` — Pluralization: "3 variants"
- `plural(items, singular, plural)` — Returns singular or plural suffix
- `code_to_english(text)` — Filter: code_style → Pretty Labels

### ui_help.py
- `page_help(page_id, title, show_title, header_tag)` — Loads help text from `static/page_help/{page_id}.html`
- `page_help_embedded(title, help_url)` — Block tag; dynamic help content rendered from template block

### ui_menus.py
- `menu_top(context, url_name, app_name, title, method)` — Top-level menu item; checks visibility via get_visible_url_names(); auto-detects active state
- `menu_item(context, url_name, css_class, arg1, arg2, badge_count, title, icon, href, admin_only, external, method, other_urls)` — Side menu item; auto-titlecases url_name; badge display; active state detection
- `absolute_url(name, *args, **kwargs)` — URL reversal with full domain

### ui_login.py
- `login_form(context, next, form)` — Renders login form or OIDC link based on settings; maintenance mode awareness

### ui_menu_bars.py
Menu bar inclusion tags: `menu_bar_annotations`, `menu_bar_analysis`, `menu_bar_main`, `menu_bar_data`, `menu_bar_genes`, `menu_bar_patients`, `menu_bar_tests`, `menu_bar_sequencing` (conditional on SEQAUTO_ENABLED), `menu_bar_settings`, `menu_bar_classifications`, `menu_bar_variants`.

### ui_tabs_builder.py

**TabBuilderTab** (dataclass): tab_number, tab_id, label, badge, badge_status, tab_status, admin_only, url, resolved_url, param, content. Property: `active`.

**TabBuilder**: Manages set of tabs; `tabs_required` property.

**check_active_tab(tab_set, tab_id, request, context)** — Determines active tab by checking activeTab GET param, referer header, context.

---

## Utility Functions

**parse_tag(token, parser)** — Generic template tag argument parser; returns (tag_name, args_list, kwargs_dict). Used by custom block tags.

**TagUtils** — Static methods for template tag value resolution: value(), value_str(), value_int(), value_bool().

**jsonify_for_js(json_me, pretty=False)** — Converts Python objects to JS-safe JSON; handles strings/booleans/numbers/complex objects; escapes `</script>` as `<\/script>`.

---

## No URL Patterns
uicore does not define its own URL patterns — URLs registered by including apps.

---

## Integration Points

| App | Integration |
|-----|-------------|
| library | log_utils (log_level_to_bootstrap, report_message), health_check signals, preview_request (PreviewModelMixin), django_utils (get_url_from_view_path, is_ajax, require_superuser), utils (html_id_safe, emoji_to_unicode, diff_text, pretty_label) |
| snpdb | admin_utils (get_admin_url), user_settings_manager (timezone), perm_path (get_visible_url_names) |
| Django | Template system, forms, FormHelper, Http |
