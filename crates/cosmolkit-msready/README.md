# `cosmolkit-msready`

Detached implementation of `cosmolkit_msready_v1`, version/schema `1.0.0`,
behavior commitment **Native**. The public Rust entrypoint is `cosmolkit`:
`Molecule::with_msready()` or `with_msready_with_params(&MsReadyParams)`, under
opt-in `cap-msready`. Defaults, `full`, `core` and existing bundles do **not**
include this leaf. Use `--features cap-msready` or `full,cap-msready`; select
`cap-smiles` separately when constructing molecules from SMILES.

This implements a finite published MS-Ready concept, **not** historical
v2.5.8/KNIME/Indigo or EPA production equivalence. Static source assets are not
an execution oracle. Python snake_case and JS `withMsready` /
`withMsreadyWithParams` projections are registered but actual MSReady language
bindings are not delivered.

## Support catalogue and source order

Five distinct directed group/tautomer rules occur in six ordered appearances:

| Order / rule | Source | Edit |
|---|---|---|
| 1 `nitro_initial` | `N mesos (#3406)/Table Creator (#2133)`, row 3 | Neutral pentavalent nitro to `[O-][N+](R)=O`; normalize charged O roles |
| 2 `enamine_to_imine` | `Metanode (#3383)/MarvinSketch (#3381)` | C=C–N–H to C–C=N; H goes to opposite alkene C |
| 3 `enol_to_carbonyl` | `Metanode (#3412)/MarvinSketch (#2253)` | C=C–O–H to C–C=O |
| 4 `ynol_to_ketene` | `Metanode (#3412)/MarvinSketch (#2136)` | Concrete C≡C–O–H to C=C=O |
| 5 `azide_mesomer` | `Metanode (#3412)/MarvinSketch (#2256)` | R–N−–N+≡N to R–N=N+=N− |
| 6 `nitro_final` | `N mesos (#3430)/Table Creator (#2158)`, local row 3 | Same nitro chemistry, distinct staged appearance |

All paths are under `Tauto_Neut (#3652)` in the pinned inventory. Guards check
complete local bond/H/charge/aromatic/radical state and validated product
valence. Unmatched branches and R attachments remain unchanged. Proton-moving
centers are neutral/nonaromatic/nonradical; nitro/azide centers are also
nonaromatic/nonradical. Two guarded Native proton classes follow the group
rules: single-bond aliphatic N+ with removable ordinary H becomes a neutral
amine; eligible O− singly bonded to C/P/S receives H. Balanced internal
transfer precedes residual reservoir exchange. Permanent/aromatic onium,
nitro/azide-required charges and blocked isotope-only donors remain charged.

Native adaptations preserve R correspondence (a slot may be filled by valid
implicit/bracket/explicit H), ignore drawing E/Z orientation, transfer D/T as
labelled state, and rank nitro O roles by isotope/chemical invariants, not
input IDs. Matches are reacquired after edits; overlapping candidates use the
lexicographically smallest valid isotope-preserving canonical graph. The six
appearances plus proton classes repeat to an exact fixed point, unlike
historical one-pass processing. `max_edits` defaults to 1,024 committed edits
per original input, accumulated across disconnection, all components, group
rules and proton edits; new sweeps never reset this edit count.

`max_matches_per_sweep` defaults to 4,096 candidate attempts, **not** a
total-input candidate cap. The initial disconnection scan is a separate sweep.
For each retained component, a group sweep pools candidates across all six
ordered appearances and their repeated rematching; it resets only when a new
complete six-appearance sweep starts, not between rules or individual edits.
A proton-engine invocation starts a separate sweep pooling internal transfer,
reservoir deprotonation and reservoir protonation, including rematching within
all three phases. Subsequent component/group/proton sweeps start fresh candidate
counts. Thus fixed-point processing may perform multiple candidate sweeps while
still consuming one input-wide edit allowance. Enumeration never silently
truncates. Limits and repeated distinct states fail with typed limit/cycle
errors, not partial success.

The classified 115 connected source appearances comprise 6 adopted / 105
unadopted normalization appearances and 3 adopted / 1 unadopted pre-component
appearances. Omitted rules are outside v1, conflicting radical repair, or
unmodeled query semantics. The sulfur/iodine rewrite and separate 2024 QSAR
supplement are not adopted. See the immutable
[source classification](../../testdata/msready/rule_classification.json) for
full paths/hashes, and `src/data/rules.json` and `src/data_exclusion_variants.json`
for packaged rule/finite-template data.

## Exactly selected metals and disconnections

Three selected Li/Na/K–O disconnections require a nonaromatic/nonradical metal
with degree 1 and no H, and an ordinary single bond to nonaromatic/nonradical O.
O has either one other heavy neighbor C/P/S on a single bond and no H, or no
other heavy neighbor and exactly one H (hydroxide). No other neighbors except
permitted explicit H are allowed. Charge pairs are exactly `(0,0)→(+1,−1)`
and `(+1,−1)→(+1,−1)`. Cutting preserves atoms/H/isotopes/net charge and validates
both resulting fragments. Only the recognized charged bonded representation
may defer pre-cut strict valence. Other charge pairs, coordination/multiple
bonds, metal-bound H and broader metal disconnection are unsupported.
Isolated q+1 degree-0 H-free Li/Na/K spectators are accepted; unlabelled ions
are excluded by explicit templates, not a generic salt filter.

Retained covalent metal tuples `(element, charge, degree, total H)` are exactly
`(Hg,+1,1,0)`, `(Hg,0,2,0)`, `(Sn,0,4,0)`. Every neighbor is a definite neutral,
nonaromatic/nonradical carbon on an ordinary nonaromatic single bond, with
complete validated carbon/component valence. No metal query/radical/aromatic
state, H, coordination or unspecified bond is accepted. Supported isotope
labels and C–metal connectivity remain in product topology even if InChI
collapses distinctions. Other valid metal/complex states are structured
unsupported chemistry; invalid valence and owner/identifier failure on an
accepted state remain failures, never reclassified unsupported for convenience.

## Exact exclusions, including nitrate

All 32 pinned source rows apply to standalone inputs and mixtures; both
carbonate IDs remain. Whole-component comparison checks cardinality, complete
bijection and exact bond/charge/H/radical/isotope state. Only stereo/maps/drawing
metadata are ignored, with owner-backed ordinary explicit-H folding on a
comparison copy. Analysis graphs are not rewritten to force exclusions.

The finite catalogue includes explicit water/ammonia/hydrogen-halide proton
states; specified singleton-metal ionic variants; reviewed terminal-O states
of phosphate/sulfate/carbonate/hypophosphite retaining P-bound H; raw PF6/PF6−
and source BF4−. Other rows match their stored state. This is not general
uncharging, solvent/salt removal or a formula/substructure filter. Labels are
never erased. `[S-2]` does not imply HS−/H2S; atomic halogen radicals do not
match hydrogen-halide templates.

Row 83 has **exactly three** comparison-only templates, each with one N, three
terminal O, one N=O and two N–O bonds, no isotope/radical/query/aromatic state:

| Template ID | Graph | Formula / charge |
|---|---|---|
| `83_raw_h2no3` | `ON(=O)O` | H2NO3 / 0; verbatim raw depiction, exclusion-only strict-N-valence exception |
| `83_nitric_acid` | `[OH][N+](=O)[O-]` | HNO3 / 0; validated adopted comparison template |
| `83_nitrate` | `[O-][N+](=O)[O-]` | NO3 / −1; validated adopted comparison template |

`nitrate_comparison_charge_adaptation_v1` records source→acid N +1, one O −1,
H −1 (net charge delta 0); source→nitrate N +1, both single O −1, H −2 (net
charge delta −1). These are catalogue metadata, **not chemical edits or proton
reservoir events**. No generated nitrate proton ladder, labelled template,
organic nitrate, nitrite or nitro exclusion is introduced. Unmatched invalid
neutral-N variants fail valence. Raw nitrate/PF6 exceptions are exclusion-only;
no invalid unmatched input survives as a retained product.

## Values, preservation and failure

The graph boundary borrows topology/properties; the owner also provides a text
boundary reusing the SMILES parser. Original graphs/properties and input
linkage remain immutable. Text parsing preserves actual text separately from
optional archival metadata. Non-exclusion text components initialize radicals
through the core owner before guards; explicit radicals are preserved. Exact
pending charged-alkali repairs remove only recognized bonds on a transient
**observation copy** before inference, with stable row correspondence. This
Native adaptation never substitutes the copy for analysis/original topology or
logs fictitious disconnections. Graph callers' supplied radical counts remain
unchanged. Radicals are retained and flagged, never closed-shell repaired.

Derived products physically clear atom/bond/enhanced/axial stereo, directions,
maps and perception carriers and use no original coordinates for identity.
Canonical isomeric SMILES preserves isotope labels, with stereo already absent.
D/T stays labelled; ordinary explicit H folds with composed dispositions.
Selected owner valence/aromaticity/H operations are reused, without generic
cleanup, organometallic cleanup or an InChI graph roundtrip.

Every product requires real Standard InChI with default options, a full
27-character key, engine version/diagnostics, exact normalized canonical SMILES,
isotope-aware formula/composition, formal charge and owner-backed **neutral-atom
formula exact mass**. Mass is not electron-corrected ionic mass, m/z or an adduct
calculation. Identity groups by `(full key, exact SMILES)` and sorts by that
pair; no cross-toolkit SMILES-byte equivalence is promised. All eligible
components, including organic counterions, survive: no largest-fragment,
implicit-carbon, solvent, mass or size filter. Dedup retains every origin and
multiplicity with graph-aligned product→input mappings.

`MsReadyStatus` distinguishes products from genuine `NoEligibleComponent`.
Every original atom has a retained/folded/reservoir-deleted/excluded disposition.
Exclusions retain source ID/hash/selected field/template/adaptation/comparison
metadata, matched formula/charge and original mappings. Events report rule,
stage, original rows, H/charge deltas and carrier. Structured failures retain
original input, stage/rule/component/atom context and completed evidence, never
a partial successful product collection. Invalid input/valence, unsupported
chemistry/identifier, owner, invariant, limit and cycle failures remain distinct.

## Provenance and runtime independence

Source: `kmansouri/MS-ready` v2.5.8 (2022-02-02), McEachran et al. (2018), DOI
`10.1186/s13321-018-0299-2`; archive SHA-256
`042516c36088ffa8b9e9ae8e6e2e9154188cdb790c8ac3c510643bb7b01acc5c`.
The supplied archive README attests MIT repository terms and original article
data terms including CC BY 4.0 / applicable CC0; this is retained attribution,
not independent online legal verification. Full source/selection hashes, URLs
and original terms are in [fixture provenance](../../testdata/msready/README.md)
and its manifest. Run information exposes packaged catalogue hashes.

Builds/runtime/examples use package-owned immutable data, never external
absolute fixture paths. Development-only preparation is local-path-only,
validates immutable hashes and fails on missing/corrupt sources. The
[approved plan](../../dev/plans/msready_pipeline_plan.md) records gates;
profile tests and static evidence do not establish historical parity.
