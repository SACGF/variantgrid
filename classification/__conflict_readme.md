conflict_populate : Management command to populate all Conflict objects that should exist
classification_grouping.py AlleleOriginGrouping.update calls conflict_services.calculate_and_apply_conflicts_for

ConflictLab links the ClassificationGrouping to the Conflict

TODO:
Need to introduce cross contexts into conflicts
Two approaches:
    1) Make cross conflicts their own kind of Conflict (very annoying as a user to have to go through cross-conflict as a separate kind of conflict to the one they're interested in)
    2) Allow ClassificationGroupings to belong to multiple Conflicts in multiple contexts

Assuming we're going for option 2:
Have to work out how cross conflict Pathogenic, changes the overall value of a VUS vs LO - which difference takes priority?
Could even have 2 different values - 1 for differences purely within and 1 for cross context - alternatively just have to work out priority.

Need to make sure cross conflict concerns don't create a conflict - or at least if it does that conflict is marked as invalid.
e.g.

Somatic Non-Cancer Oncogenic would be a primary part of the Somatic Non-Cancer Onc/Path conflict, but a cross context for
Germline Onc/Path conflict. Yet if there are no Germline classifications it seems noisy to create one.

Probably easier to create the conflict but be able to mark the Conflict as invalid unless there are non cross conflict contributions.
So would have several DataRows that would also indicate that they're cross conflict

It does mean the results of a cross conflict could appear in 2 or multiple conflicts, but that's probably a good thing