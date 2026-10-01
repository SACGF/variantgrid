# review — research notes

Verified against a96540a68 on 2026-09-25

`review` records a multi-lab discussion about a shared object - in practice an Overlap (and, read-only now, a legacy
DiscordanceReport) in `classification` - and hands off to the owning app for the follow-up action. It owns no
classification logic: the form captures who met, how, which labs took part and a set of per-question "difference"
answers, all stored as JSON on `review/models/review_models.py:Review`. Models and URLs are in
[models](../maps/models.md#review) and [urls](../maps/urls.md#review); the one signal is listed in
[signals](../maps/signals.md). There are no tasks or commands.

## Flows

### Starting and saving a review

The owning app makes its model reviewable with `review/models/review_models.py:ReviewableModelMixin`, which adds a
nullable `reviews` FK to a `review/models/review_models.py:ReviewedObject` (created lazily by `reviews_safe`). Its
view then redirects to `start_review` with the ReviewedObject pk and a topic key - both callers hard-code
`"discordance_report"` (`classification/views/overlaps_view.py:overlap_report_review`,
`classification/views/discordance_report_views.py:discordance_report_review`).

`review/views/review_views.py:new_review` builds an unsaved Review, checks view and write permission and renders
`review/views/review_views.py:ReviewForm`. The form's lab choices come from the source object's `reviewing_labs`
(`review/widgets/multi_lab_selector.py:MultiChoiceLabField`, pre-selecting the user's default lab if it is one of them),
and one `DescribeDifferenceField` per enabled `ReviewQuestion` of the topic. `clean` requires at least one question
answered. `ReviewForm.save` writes `meeting_meta = {"participants": {"review_method", "review_participants"},
"answers": {question_key: DescribeDifference JSON}}`, sets the M2M labs (saving first if new), logs via
`log_admin_change`, and the view redirects to `Review.next_step_url()` - i.e. the source object's `post_review_url`.

### Completing: the action page lives in the owner

The post-review step is not in this app. `classification/views/overlaps_view.py:action_overlap_review` (and the
discordance equivalent) either postpones or applies per-lab value changes, then calls
`Review.complete_with_data_and_save(post_review_data)`. Once `is_complete`, `review/views/review_views.py:edit_review`
serves the read-only detail instead of the form. To render `post_review_data`, `Review.post_review_data_formatted`
sends `review_detail_signal` with the source object's class as sender; receivers are
`classification/signals/overlap_review_formatting.py` and `classification/signals/discordance_report_review_detail.py`,
falling back to raw JSON.

### Permissions

View: `Review.can_view` defers to `source_object.can_view` if it exists, else allows - Overlap has none, so overlap
reviews are visible to every logged-in user (matching the overlap pages). Write: `ReviewableModelMixin.can_review` -
the user (admin-checked) belongs to one of `reviewing_labs`. `edit_review` falls back to the detail view for
non-writers, complete reviews, or `is_review_locked` sources (DiscordanceReport always returns True: reviews are
being moved to Overlaps). Covered by `classification/tests/views/test_overlap_review_permissions.py`.

## Why it is shaped this way

- `ReviewedObject` is an indirection instead of a GenericFK or per-model M2M: each reviewable model carries one FK to
  it, and `ReviewedObject.source_object` finds the owner by scanning its reverse `*_set` managers
  (hardened for Django 6 in ab4de89de). Adding a new reviewable model needs no change here.
- Topics and questions are data (edited in the admin, `review/admin/review_admin.py:ReviewTopicAdmin` with a question
  inline), so wording can change per deployment. The only question type is `QuestionValueType.Disagreement`.
- `ReviewMedium` / `ReviewParticipants` are fixed choices with a free-text "other"; stored values that are not a choice
  are read back as `ValueOther(key="other")`. `review_method` was single-select originally, and `Review.review_method`
  still accepts a string.

## History

Built May 2023 for discordance reports (0c4f089ae; "only one review per discordance report" 7c7ac005e), with the UI
wording later changed from "Review" to "Discussion" (ebdbfb9f1) while code names stayed. The Overlaps rework (September 2026,
25ed4f32f) made `Overlap` reviewable and locked DiscordanceReport reviews. Write permission by reviewing-lab
membership came in e26588083 / e54821fba (variantgrid_private#3827).

## Traps

- `review/migrations/0005_discordance_report_topic.py` seeds the `discordance_report` topic both callers hard-code,
  with a default question set; a deployment that had already created it in the admin keeps its own name and questions.
  Without the row the "review" buttons 500 on `ReviewTopic.objects.get`.
- Overlaps start a new Review each time (resume is commented out), DiscordanceReports resume the first existing one.
- `Review.user` means "last actor", not author: `ReviewForm.save` sets it to the editor (and logs as them), and the
  action views overwrite it with whoever completed the review.
- Editing a review keeps its stored `review_date`; only a new review defaults to today.
- `Review.answers` skips answers whose question was deleted in the admin; disabling (`enabled=False`) is still the
  way to retire a question, since a deleted one silently drops out of historical reviews.
- Question `key` is the table-wide primary key, not per topic: prefix it with the topic when adding questions.
