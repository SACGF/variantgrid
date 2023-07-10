# Discordances

Discordances are only currently enabled in the Shariant environment. The definition of a discordance follows the one used by ClinVar - by dividing clinical significances into 3 buckets:

Bucket 1: Benign, Likely Benign
Bucket 2: VUS, VUS-A, VUS-B, VUS-C
Bucket 3: Likely Pathogenic, Pathogenic

If for any given allele with multiple classifications, if they fall into two or more buckets, a discordance is raised.

## Clinical Groupings
_formerly clinical contexts_

While not heavily exposed in the interface currently, published classifications are assigned to a clinical grouping, and the clinical groupings are assigned to an allele.

This will allow us to sub-divide classifications within the same allele by condition, if the conditions were significantly different.

## Discordance Report

A DiscordanceReport marks the start, and optionally the end, of a discordant state for a clinical grouping.
The general hope is, when the relevant labs are presented with a discordance report, they will discuss between themselves how to solve it.

## Review

VariantGrid can support many kinds of Reviews, though as of this writing (July-2023) the only review type supported is for DiscordanceReports. This gives the users a template where they can record why they believe a discordance occurred in the first place and how it was (or wasn't solved).

## Future Development

A DiscordanceReport is opened when a discordance is first detected, and closed when it's resolved - or when users agree that it can't be solved. Determining if a DiscordanceReport is at "pending concordance" is a live calculation rather than stored in the database.

These can result in some messy edge cases - we would prefer it if a DiscordnaceReport is marked as "can't be solved", it technically doesn't close but is marked as a certain status.