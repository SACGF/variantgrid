# MatchMaker Exchange — onboarding / admin checklist

Hand-off list for the lab head / research assistant. This is the **non-code** onboarding;
the software integration is tracked in `claude/matchmaker_exchange_plan.md`.

## First: decide the scope

The technical connection to MatchMaker Exchange (MME) is **bilateral** — you exchange an
authentication token directly with each node you want to talk to. How much formal
onboarding is needed depends on what we want:

- **Path A — connect with a specific partner node** (lightest). E.g. a collaborator
  running seqr / PatientMatcher. Mostly a direct token exchange with that node; no
  central Steering Committee process (the partner still sets their own terms).
- **Path B — become a recognised, discoverable MME Service.** Needed to be listed in the
  network and to serve inbound reciprocal matches to the established nodes (GeneMatcher,
  DECIPHER, etc.). This is the layer with central oversight and agreements.

The big well-known nodes generally expect Path B before they'll issue a token, so if the
goal is broad querying, plan for Path B.

## Path A — bilateral partner connection

- Identify the partner node and their technical contact.
- Agree terms directly with them (their data-sharing / use conditions).
- **Exchange auth tokens**: they issue us a token to query them; we issue them a token to
  query us. Hand both tokens to our technical contact.
- Agree public contact details to expose (see "Common tasks" below).

## Path B — become a recognised MME Service

- **Make first contact**: email **api@matchmakerexchange.org** stating we (SA Pathology /
  VariantGrid) want to connect a database as a new MME node. Ask for the current
  onboarding pack.
- **Meet and sign the three governance documents**:
  - **Informed Consent Policy**
  - **End User Agreement**
  - **Service Requirements** (2026 revision)
  - Signing needs someone with authority to sign for the organisation (lab head or above).
- **Steering Committee review**: new sites are reviewed/approved by the MME Steering
  Committee. Find out their timeline and whether they expect an application or short
  presentation.
- **Serve inbound matches**: being a Service means other labs can query us — confirm we
  are willing and cleared to do this (see ethics below).

## Common tasks (both paths)

- **Nominate contacts**:
  - a **technical contact** for the API integration, and
  - a **scientific/clinical contact** — the person other labs email when a patient
    matches. This becomes the public `contact` on every submission.
- **Decide public contact details** to publish on outbound submissions: a **name**, an
  **institution**, and a **contact URL or email**. No patient-identifying information.
- **Local ethics / privacy sign-off**:
  - Sharing curated variant + condition data (and HPO phenotype terms where present) to
    external third-party databases.
  - Our rule that only classifications shared at **"3rd Party Databases" level** (the same
    bar as ClinVar) are eligible to go to MME.
  - Serving inbound queries from other labs, if we take Path B.
- **Per-node token exchange**: for each node we connect to, complete their connection
  agreement and swap tokens. Give the tokens to the technical contact — these are the only
  credentials the software needs.
- **Test then go live**: connect in **test mode** against each node's test instance first;
  switch to production only after sign-off.

## Notes for whoever picks this up

- **This is the critical path.** The software can be ready well before the agreements,
  token exchange, and any Steering Committee approval clear — so start the first contact
  (Path A partner, or api@matchmakerexchange.org for Path B) early.
- **Timelines are set by MME and the other nodes**, not us. The intake email is what
  unblocks the rest.
