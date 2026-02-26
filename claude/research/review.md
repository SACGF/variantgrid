# VariantGrid Review App — Reference Document

## Purpose and Overview

The `review` app implements a collaborative review workflow system for multi-lab discussions on shared objects (particularly variant classifications and discordance reports). It allows teams to document review meetings, record decisions with detailed reasoning, and track post-review actions.

**Key Design:** Structured meeting documentation with templates; multi-lab participant tracking; signal-based domain-specific processing.

---

## Models

### ReviewTopic
Template for a type of review (e.g., "Discordance Discussion").
- Fields: key (PK, text), name, heading
- Methods: questions property → ordered list of enabled ReviewQuestions

### ReviewQuestion
A specific question prompted in a review form.
- Fields: topic FK, key (PK), label (question text), help (nullable), heading (groups related questions), order (int), value_type (QuestionValueType, currently only 'D'=Disagreement), enabled (bool)

### ReviewedObject
Container representing an object being reviewed.
- Fields: label (display name/description)
- Methods:
  - `new_review(topic, user)` → Review: Create new Review for this object
  - `source_object` (cached_property): Lazy-loads actual related object (Allele, Classification, etc.) via generic relation

### Review
Actual documentation of a review meeting.
- Fields:
  - reviewing FK (ReviewedObject), topic FK (ReviewTopic), user FK (PROTECT)
  - review_date (DateField)
  - reviewing_labs (M2M → Lab)
  - meeting_meta (JSONField) — core review data (see structure below)
  - is_complete (bool) — marks ready for post-review action
  - post_review_data (JSONField) — populated after completion
- Methods:
  - `as_json()` → dict
  - `can_view(user_or_group)` → bool (delegates to source_object)
  - `check_can_view(user)` — raises PermissionDenied
  - `next_step_url()` — calls source_object.post_review_url()
  - `review_method` property → list of ValueOther
  - `participants` property → list of ValueOther
  - `answers` property → list of ReviewAnswer objects parsed from meeting_meta
  - `post_review_data_formatted` (cached_property) — formatted via review_detail_signal
  - `complete_with_data_and_save(data)` — marks complete

### ReviewableModelMixin (Abstract)
Mixin for models that can be reviewed.
- Fields: reviews FK (ReviewedObject, nullable)
- Methods:
  - `reviewing_labs` property — returns {self.lab} by default; subclass for multi-lab
  - `reviews_safe` property — lazy-creates ReviewedObject
  - `post_review_url(review)` → str — URL after review (default: get_absolute_url())
  - `reviews_all()` → QuerySet[Review] — in reverse date order
  - `is_review_locked` property — default False
- Known implementations: **DiscordanceReport**

---

## Enums and Data Classes

### ReviewMedium
- `email = 'email'`
- `phone = 'phone'`
- `video = 'mtm'` (Multi-disciplinary Team Meeting)

### ReviewParticipants
- `curation = 'curation'` (Curation Scientists)
- `clinicians = 'clinicians'`
- `external_experts = 'external_experts'`

### ValueOther
Represents a selected value with "other" option support.
- Fields: key (identifier or "other"), label
- Static method: `from_str(value, text_choices_class)` — returns as "other" if unrecognized
- Sortable (other values go to end)

### ReviewAnswer
Represents one answer in a review.
- Fields: question (ReviewQuestion), details (text), resolution (DifferenceResolution enum value)

---

## Views and Forms

### ReviewForm
- Fields: review_date (DateField), review_method (MultiChoiceFieldWithOther), review_participants (MultiChoiceFieldWithOther), reviewing_labs (MultiChoiceLabField), dynamic question fields (DescribeDifferenceField)
- Validation: at least one question must have an answer
- Save logic: serializes to meeting_meta JSON, associates reviewing_labs M2M, logs change
- Init: populates from existing review.meeting_meta if editing; defaults to user's default lab

### View Functions

- `new_review(request, reviewed_object_id, topic_id)` — creates new Review
- `edit_review(request, review_id)` — edits existing Review; checks can_view() first
- `view_discussion_detail(request, review_id)` — AJAX endpoint for viewing review detail template
  - Query params: edit, show_source_object, show_outcome
- `_handle_review(request, review, reviewing=None)` — internal helper for GET/POST

---

## URL Patterns
```python
path('new/<int:reviewed_object_id>/<str:topic_id>/', views.new_review, name='start_review'),
path('detail/<int:review_id>', views.view_discussion_detail, name='review_detail'),
path('<int:review_id>/', views.edit_review, name='edit_review'),
```

---

## Admin Interface

- **ReviewTopicAdmin**: id, key, name, heading; ReviewQuestionInline (editable, ordered by order field)
- **ReviewAdmin**: default ModelAdminBasics display

---

## Signal

### review_detail_signal
- Sender: ReviewedObject model class (source_object's class)
- Args: instance (Review)
- Purpose: Allow source objects to provide formatted representation of post_review_data
- Usage: Classification/discordance app implements handler to format domain-specific outcomes

---

## Data Storage Structure

### Review.meeting_meta JSON
```json
{
  "participants": {
    "review_method": ["email", "phone"],
    "review_participants": ["curation", "clinicians", "other_value"]
  },
  "answers": {
    "question_key": {
      "details": "Discussion details here",
      "resolution": "AGREED"
    }
  }
}
```

---

## Workflow

1. **Identify Object** → get ReviewedObject (or use reviews_safe on ReviewableModelMixin object)
2. **Start Review** → `/review/new/<object_id>/<topic_id>/` — displays ReviewForm with topic questions
3. **Document Meeting** → fill form: date, method, participants, lab(s), answers to questions
4. **Save Review** → meeting data serialized to meeting_meta JSON
5. **Post-Review Action** → review.next_step_url() redirects to source_object's post-review page

---

## Integration Points

| App | Integration |
|-----|-------------|
| Classification/Discordance | DiscordanceReport uses ReviewableModelMixin; reviews document how labs resolved disagreements; review_detail_signal hooks into classification formatting |
| Lab system | M2M reviewing_labs; lab permission system controls access |
| User settings | UserSettings.get_for_user() for default lab pre-population |
