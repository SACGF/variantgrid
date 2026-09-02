""" Migration helper for composite grid columns - takes `apps`, never live models.

    @see snpdb.models.models_columns.CompositeColumnMember """


def collapse_into_composite(apps, composite_id: str) -> None:
    """ In every collection: replace the run of this composite's members with the composite, at the
        position of the first member (or the composite itself, if already present). A collection with
        none of them is untouched. Bumps version_id on each collection changed (historical models skip
        CustomColumn.save(), which is what normally bumps the node grid definition cache key) """
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    members = set(CompositeColumnMember.objects.filter(composite_id=composite_id)
                  .values_list("column_id", flat=True))
    if not members:
        raise ValueError(f"Composite column '{composite_id}' has no members")

    for ccc in CustomColumnsCollection.objects.all():
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        ordered = list(cc_qs.order_by("sort_order").values_list("column_id", flat=True))
        group = [c for c in ordered if c == composite_id or c in members]
        if not group:
            continue

        position = ordered.index(group[0])
        new_order = [c for c in ordered if c not in group]
        new_order.insert(position, composite_id)
        if new_order == ordered:
            continue

        cc_qs.exclude(column_id__in=new_order).delete()
        CustomColumn.objects.get_or_create(custom_columns_collection=ccc, column_id=composite_id,
                                           defaults={"sort_order": position})
        for sort_order, column_id in enumerate(new_order):
            cc_qs.filter(column_id=column_id).update(sort_order=sort_order)

        ccc.version_id += 1
        ccc.save()
