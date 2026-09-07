# Classify-queue tags by property, not by name

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-07
Status: landed

`RequiresClassification` and `SomaticReportable` started life as magic tag names. `Tag` now carries the two
properties that made them special (`snpdb/models/models.py:Tag`), and the queue, the New Classification buttons and
tag resolution all already ask `Tag.requires_classification`. This plan finishes the job: every remaining place that
names a tag asks the tag's properties instead, and the somatic wizard in
[`somatic_curation_reuse_issue_1419_plan.md`](somatic_curation_reuse_issue_1419_plan.md) is specified on the same
terms so it lands without reintroducing a named tag.

## The data

Unchanged. The two fields are the whole design; a lab configures them per tag on the tag settings page
(`snpdb/templates/snpdb/settings/tag_settings.html`, `snpdb/tag_operations.py`).

```python
class Tag(models.Model):
    id = models.CharField(max_length=50, primary_key=True)
    retired = models.DateTimeField(null=True, blank=True)
    merged_into = models.ForeignKey('self', null=True, blank=True, on_delete=SET_NULL)
    # Both / Germline only / Somatic only - which classifications the tag is about
    allele_origin_bucket = models.CharField(max_length=1, choices=TAG_ALLELE_ORIGIN_CHOICES,
                                            default=AlleleOriginBucket.UNKNOWN)
    # A to-do: tagging asks for a classification, and classifying the variant resolves the tagging
    requires_classification = models.BooleanField(default=False)
```

A **classify-queue tag** is a live tag with `requires_classification=True`. A classify-queue tag **for a bucket** is one
whose `allele_origin_bucket` is that bucket or Both, which is the same rule
`snpdb/models/models_enums.py:AlleleOriginFilterDefault` already applies to its `buckets`. Both queries get a home on
`Tag` as classmethods next to `live_qs`, so callers share one definition.

### "SomaticReportable survives the classification" needs no third field

sapath#246 says the somatic lab keeps the tag on the variant after classifying. Resolution already keeps it:
`analysis/variant_tag_operations.py:resolve_variant_tag` stamps `resolved` and `resolved_classification` rather than
deleting, the row still counts on the tags page and grids, and the tags node shows it again with
`TagNode.include_resolved`. The wizard follows the same rule. If a lab later wants a queue tag that stays *visible* in
the work lists after classification, that is one more boolean on `Tag`, set from the same settings page, and this plan
leaves room for it by keeping every lookup on `Tag` fields.

## Part A - runtime code stops reading `settings.TAG_REQUIRES_CLASSIFICATION`

The setting (`variantgrid/settings/components/default_settings.py:766`) stays, because three pushed migrations read it
to seed the tag (`snpdb/migrations/0002_initial_data.py`, `snpdb/migrations/0089_one_off_user_tag_colors_to_collection.py`,
`snpdb/migrations/0251_one_off_tags_requiring_classification.py`). Its comment says so: it is the name of the tag a
fresh install is seeded with, and behaviour comes from `Tag.requires_classification`.

Every non-migration reader today is a test (`analysis/tests/test_variant_tags.py`,
`variantopedia/tests/test_tagged_variant_grid.py`, plus the 27 test references to the two names). They move to a fixture
builder, `create_classify_queue_tag(tag_id="ToDo", bucket=AlleleOriginBucket.UNKNOWN)`, in the tags section of
`claude/guides/testing.md`, so a test that needs a queue tag makes one with the flag rather than assuming a name has it.
Tests that check the somatic path make one with `bucket=AlleleOriginBucket.SOMATIC`. Tag ids in tests become
descriptive of the property under test (`ToDo`, `SomaticToDo`) so a reader sees the property is what matters.

## Part B - the help text names whatever tags are configured

Three classification list templates tell the user to tag a variant "as RequiresClassification"
(`classification/templates/classification/classifications.html:183`,
`classification/templates/classification/classification_groupings.html:120`,
`classification/templates/classification/classifications_legacy.html:197`). A simple template tag in
`snpdb/templatetags/` renders the live classify-queue tag names, so the sentence reads "by tagging it as
RequiresClassification or SomaticReportable" from the database, and the tagging clause is left out entirely when no
queue tag exists. The three templates use it.

## Part C - the 1419 plan is specified on properties

Edits to `somatic_curation_reuse_issue_1419_plan.md`:

- **Part C, stage 1**: the wizard's rows are the unresolved taggings (`VariantTag.unresolved_q`) in scope whose tag is a
  classify-queue tag for the somatic bucket. Dropping a row deletes the tagging as now. Classifying resolves the
  tagging exactly as the germline queue does, which keeps the row.
- **Part D**: the launch points count and scope by the same query. The `TAG_SOMATIC_REPORTABLE` setting and the
  `formatVariantTagFirstColumn` sentence go, replaced by one line saying the launch points use the `Tag` classmethods.
  A deployment turns the feature on by flagging a tag as a somatic queue tag on the settings page, which SA Path already
  has from migration 0251.

## Part D - docs

- `claude/domain.md` gains a **Classify-queue tag** entry next to **VariantTag**, defining the property and the bucket rule.
- `snpdb/CLAUDE.md` gets one line: the setting is seed data only, behaviour is on `Tag`.
- `analysis/CLAUDE.md:61` already says the buttons ask the property; it drops the parenthetical naming `SomaticReportable`
  as the example once the wizard plan no longer needs it.

## Tests

`scripts/vg tests --explain` after each part. The tests that earn their keep are the existing queue and resolution
tests rewritten on the fixture, one test for the bucket rule on the `Tag` classmethod (Both matches either bucket,
Germline only does not match somatic), and one for the template tag's empty and populated renderings.

## Done when

- `grep -rn "TAG_REQUIRES_CLASSIFICATION"` outside `migrations/` and the settings file returns nothing.
- `grep -rn "RequiresClassification\|SomaticReportable"` outside migrations returns only the seed name in the settings
  file and prose in plans.
- `scripts/vg docs check` passes for this plan and the 1419 plan; `scripts/vg map --check` passes.
