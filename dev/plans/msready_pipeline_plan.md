# Implementation Plan

## Goal
Deliver an audited, isotope-preserving, stereo-free published-concept MSReady pipeline in detached `cosmolkit-msready`, an opt-in thin Rust `Molecule` API, a validated local commit, and local PR notes without remote actions.

## Revision disposition and evidence

This replaces the original `00edfd4c-7a95-45f3-a1b1-f2d9ee8ef332/plan.md`; the context-builder handoff in that directory remains authoritative. Retain the original five-rule Native profile, early real-InChI proof, descriptor-owned neutral-atom mass helper, provenance model, source discipline and mandatory validation. Resolve the oracle's corrections as follows:

1. `cap-msready` is **opt-in only**. Do not add it to `full`, `core`, any existing bundle or `default`; do not modify normative API/architecture documents. Explicitly enable it in its validation commands/example.
2. Nitrate row 83 has **three separate, exact comparison templates**, specified below. The raw neutral `H2NO3` depiction is not described as conventional nitrate. Adopted nitric-acid/nitrate templates are documented comparison-only adaptations, not O-only proton variants or runtime repairs.
3. Freeze exact local graph-state predicates for retained Hg/Sn and neutral/charged alkali disconnection below, with accepted/rejected cases, owner validation and blocking real-identifier proof.

No broad repository exploration was repeated. Targeted reads were necessary because the handoff/original plan omitted these precise contracts: `dev/public_api_design.md` §7 confirms unchanged nine-bundle `full`; archive `exclusion_list_v2_5_8.json` confirms row 83 selects `CANONICAL_SMILES=ON(=O)O`, formula `H2NO3`; `crates/cosmolkit-core/src/periodic_table.rs:607,637` supplies Sn valences 2/4 and unconstrained Hg valence -1; `src/sanitize.rs:578-650` confirms borrowed selected sanitization and separable generic/organometallic cleanup. Hg's unrestricted owner valence entry is **not** a substitute for the narrower v1 predicate. No build, test, Git or network command was run by this planner. During Task 1 source verification, the inventory-attested `final_nitro` source was corrected from Table Creator (#2133) to `Tauto_Neut (#3652)/N mesos (#3430)/Table Creator (#2158)/settings.xml` (manifest hash verified); its classified source appearance remains local order 3 with the same raw reaction as initial nitro. This correction changes source attribution only, not chemistry or stage order.

## Authority and execution contract

- Worktree: `/data/dev/COSMolKit-msready`; branch: `feat/msready-pipeline`; local-main base: `d39640716a6c9cf007350288349a0604fde0936a`.
- Reference: `/data/dev/COSMolKit/MSREADY_SMILES_implementation_reference.md`; extracted source root: `/data/dev/COSMolKit-msready-reference/MSREADY_reference_data/`. ZIP SHA-256: `042516c36088ffa8b9e9ae8e6e2e9154188cdb790c8ac3c510643bb7b01acc5c`.
- User authorizes the new crate/branch, local commit, and narrowly scoped **uncommitted** `AGENTS.md` exception to use a dedicated task plan. All other repository rules hold. Neither checkout's exception may be staged or committed. Preserve unrelated edits and existing architecture/standards/completion ledger.
- The parent/oracle APPROVE verdict authorizes this plan. This persisted copy is the sequential execution ledger; Steps 1–7 (Tasks 1–3) are completed at this handoff and all later steps remain unchecked. Completion of this source/fixture seam is not a chemistry implementation or runtime test pass.
- Once approved, Task 1 persists this dedicated plan as the task's sole sequential ledger. Each unchecked step executes in order, continues until complete/blocked/interrupted, and is marked only by `[ ]` → `[x]` on completion. Do not skip or summarize required Read Steps, rely on earlier reading, or implement stubs/placeholders/partial supported behavior. After compaction reread source reproduction rules.
- Every real action has an independent immediately preceding Read Step for `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`; each task below additionally lists its applicable focused reading. After writing tests, the next real task runs their most specific validation. Final checks do not replace these immediate gates.
- Offline environment for every command: `CARGO_NET_OFFLINE=true`, `UV_OFFLINE=1`, `PIP_NO_INDEX=1`. No fetching, dependency downloads, submodule update, remote lookup, push, publication or remote PR. Missing cached tools/dependencies/fixtures are explicit blockers, never authorization to use the network.

## Frozen v1 contract

### Profile, rules and limits

Profile ID `cosmolkit_msready_v1`, version/schema `1.0.0`, behavior commitment `Native`: a project implementation of the published concept, **not** v2.5.8, Indigo or EPA production equivalence.

Retain all eligible connected components, including organic counterions. Apply all 32 pinned exclusion rows to standalone inputs too; retain both duplicate carbonate IDs. No largest-fragment, generic solvent/salt, mass/size or implicit carbon-presence filter.

Adopt **five distinct group/tautomer rules in six ordered source appearances**:

| Order | Project rule | Pinned source | Frozen chemistry |
|---|---|---|---|
| 1 | `nitro_initial` | Inventory `initial_nitro`, Table Creator 2133, row 3, raw `[*]N(=O)=O>>[O-][N+]([*])=O` | Repair neutral pentavalent nitro and normalize its valid charge-separated target; N has one preserved substituent and two terminal O atoms. |
| 2 | `enamine_to_imine` | `mrv/Tauto_Neut (#3652)__Metanode (#3383)__MarvinSketch (#3381).mrv` | C=C–N–H → C–C=N, transferring N-bound H to the opposite alkene C while preserving each R attachment. |
| 3 | `enol_to_carbonyl` | `mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2253).mrv` | C=C–O–H → C–C=O, transferring O-bound H to the opposite C. |
| 4 | `ynol_to_ketene` | `mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2136).mrv` | C≡C–O–H → C=C=O, transferring H to the other C; concrete C triple-bond endpoint required. |
| 5 | `azide_mesomer` | `mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2256).mrv` | R–N−–N+≡N → R–N=N+=N−, preserving branch, formula and net charge. |
| 6 | `nitro_final` | Inventory `final_nitro`, `N mesos (#3430)/Table Creator (#2158)`, local order 3 | Same raw reaction chemistry as order 1, retained as a separate staged appearance. |

Proton-moving centers are neutral/nonaromatic/nonradical, use supported bond orders, possess the required H and yield valid product valence. Nitro/azide edit centers are also nonaromatic/nonradical. Preserve all unmatched branches, atom identities, formula/isotope composition and net charge except explicitly logged proton exchanges.

Classify all **115 connected appearances**: 6 adopted/105 unadopted normalization appearances and 3 adopted/1 unadopted pre-component appearances. Record every unadopted appearance as `outside_v1_rule_set`, conflicting radical repair or unmodeled query semantics. These are omitted profile rules, not supported-input runtime errors relabelled unsupported. Do not import the separate 2024 QSAR supplement.

Explicit, versioned Native adaptations:

- Preserve MRV R-group attachment correspondence; do not turn schematic R groups into independent wildcard atoms. A concrete attachment slot can represent available implicit/bracket/ordinary explicit H when valid, including acetone/enol convergence; record this as an adaptation, not proven Indigo behavior.
- Ignore drawing E/Z orientation for the stereo-free profile; guard aromatic centers. Transfer D/T as actual labelled atom/bond state, never erase it.
- Normalize nitro O roles using isotope/chemical invariants, including already-charged/isotope-asymmetric nitro. Do not choose by input ID.
- Re-match after edits. For overlapping matches choose the lexicographically smallest valid next stereo/map-free isotope-preserving canonical graph representation; input IDs may break ties only between chemically identical results. Candidate enumeration may not silently truncate.
- Repeat the six-rule sweep then guarded proton classes until unchanged exact normalized identity; detect repeated distinct states. Defaults: 1,024 edits/input and 4,096 candidate matches/rule sweep. Limit/cycle errors are typed failures, not partial success. Fixed-point behavior is explicitly different from historical one-pass processing.

### Exact metal-domain predicates

Definitions: `degree` counts all graph neighbors; `heavy_degree` excludes H atoms; `H_total` includes bracket H, explicit H neighbors and owner-inferred implicit H, with isotopic H counted separately. All inspected states include isotope/charge/radical/bond attributes, not just element symbols. Every accepted product requires model topology validation plus selected core valence/sanitization; no generic cleanup/organometallic rewrite may widen this domain.

#### Retained covalent Hg/Sn

At each retained Hg/Sn atom:

- No query, aromaticity, radical electrons, explicit/bracket/implicit H, H neighbors, coordination/dative/zero/unspecified bond, or unsupported atom annotations. Isotope labels on a supported state are preserved; owner isotope/mass/identifier errors are not bypassed.
- Every neighbor is a definite **nonaromatic, neutral, nonradical carbon**, joined by a nonaromatic ordinary single bond. Carbon's complete local valence and all remaining component atoms must validate. Substituents on those carbons remain unchanged.
- Accepted `(element, formal charge, degree, H_total)` tuples are exactly **`(Hg,+1,1,0)`, `(Hg,0,2,0)`, `(Sn,0,4,0)`**. No implicit completion of missing ligands or metal H is permitted. Hg is not accepted merely because its owner periodic-table valence is unconstrained.
- Preserve accepted metal–C bonds, charge, isotopes and radical state elsewhere through every subsequent stage and identity conversion. InChI is not allowed to round-trip/replace the graph.

Required positive tests: `C[Hg+].[Cl-]` → retained `C[Hg+]`; `C[Hg]C` → unchanged dimethylmercury graph; `C[Sn](C)(C)C` → unchanged tetramethyltin graph; atom permutations and isotope-labelled accepted Hg/Sn state. Required negative tests: neutral singly bonded `C[Hg]`, `C[Hg+2]`, `C[Hg]Cl`, undercoordinated `C[Sn](C)C`, Sn with a metal-bound H, metal radicals, aromatic-carbon metal attachment and coordination/multiple/query bonds. Construct graph-level fixtures where parser syntax cannot represent the tested annotation. Valid-but-unmodeled tuples return contextual `UnsupportedChemistry`; an independently invalid atom valence returns `InvalidValence`. A predicate-accepted case whose owner/identifier implementation fails is a **supported-input failure/blocker**, not a newly unsupported tuple.

#### Selected Li/Na/K–O disconnection

Only these three alkali elements are disconnected. The metal is nonaromatic/nonradical, has **degree 1, H_total 0**, and one ordinary single bond to a nonaromatic/nonradical O. Preserve isotope labels and all other atoms. Oxygen has no neighbors except the metal, optionally explicit H atoms, and one of these two ligand forms:

1. **Organic/oxo-acid ligand:** exactly one other heavy neighbor C/P/S, connected by an ordinary single bond, **H_total 0**. The C/P/S neighbor is definite/nonradical and its remaining bonds/valence are unchanged; an aromatic C neighbor is allowed because only the ordinary exocyclic O bond is edited.
2. **Hydroxide ligand:** no other heavy neighbor and **H_total exactly 1** (ordinary H or isotope-bearing H, preserved).

Allowed input charge pairs and edits are exactly:

| `(q_metal,q_O)` before | After severing the M–O bond | Accounting |
|---|---|---|
| `(0,0)` | `(+1,-1)` | Δq_M=+1, Δq_O=−1; total charge, atoms/H/isotopes unchanged. |
| `(+1,-1)` | `(+1,-1)` | Bond removal only; no charge/H edit. |

No other input pair is silently redistributed. In particular, `(+1,0)` and `(0,-1)` are outside this disconnection capability; do not invent M2+, H or electron/radical changes to preserve the total. Model topology/syntax always validates first. The **exact** `(+1,-1)` bonded representation may defer its pre-disconnection O/alkali strict valence check, since the bond depiction itself is the recognized representation repair; after removal, both fragments must pass owner validation with recorded unchanged H counts. This bounded repair is not general permissive retention. The `(0,0)` representation must satisfy its ordinary pre-edit ligand H/valence predicates. Both cases strictly validate post-edit with formula/net-charge equality before any proton edit.

Disconnected Li/Na/K ions with q=+1, degree 0, H_total 0 and no radical/query/aromatic state are accepted spectator states. Unlabelled forms match the explicit exclusions; an isotope-labelled accepted ion must not match an unlabelled template by label removal and requires real identity/mass support if retained. Neutral isolated source metals and other archived singleton states are exclusion-only unless another frozen retained-domain rule covers them.

Required positives, **for each of Li/Na/K**: neutral bonded acetate `CC(=O)O[M]`, already-charged bonded acetate `CC(=O)[O-][M+]`, ionic dotted acetate `CC(=O)[O-].[M+]`, neutral alkoxide `CO[M]`, raw source hydroxide `[M]O`, and charged hydroxide `[OH-][M+]`; labelled ligand H and metal labels test preservation/non-exclusion. Required negatives: `CC(=O)O[M+]`, `CC(=O)[O-][M]`, metal degree>1, M–C/M–N/multiple/dative bonds, O with extra heavy attachment/H, local radicals and metal-bound H. Verify exact charge/H deltas and metal origin even when subsequently excluded.

Do not adopt the sulfur/iodine source rewrite. Retained metal/complex states outside these explicit predicates are unmodeled independent capabilities, not candidates for generic disconnection. Archive whole-component exclusion matches are handled before retained-domain rejection so all 32 source entries remain testable; this does not authorize retaining their unusual valence states.

### Guarded proton/H and radical policy

Use two explicit Native proton classes from reference §7:

1. N+ deprotonation: nonaromatic/nonradical N with only single bonds and removable ordinary H, yielding validated neutral amine valence. Preserve quaternary/aromatic onium, nitro/azide N and structurally necessary charges.
2. O− protonation: nonradical O with one single bond to C/P/S and no H, yielding validated neutral O–H. Nitro O− does not qualify.

Apply to charged molecules and net-zero zwitterions. Prefer balanced internal ordinary-H transfer; log any remaining proton-reservoir exchange. D/T can transfer internally but are not removable ordinary reservoir protons: retain charge and emit a specific warning if only labelled H is available. Preserve radicals/electron counts and flag them; no radical-to-closed-shell repair. A valid center with no adopted edit remains charged, not incorrectly unsupported. Recompute owner H/valence after edits, compose maps and enforce H/charge deltas and isotope-aware formula invariants.

### Exact exclusion comparison and nitrate templates

Compare whole components on a temporary comparison graph, with equal atom/bond cardinality, a complete bijection, exact bond state and explicit charge/H/radical/isotope predicates using `cosmolkit-search`. Ordinary explicit-H folding delegates to core mappings; it never edits the analysis graph or erases isotope labels. Ignore only stereo/maps/drawing metadata for comparison. Formula equality or substructure matching is insufficient.

Preserve the **raw 32-row table** independently of reviewed production variants. Select `Structure(SMILES)` when present, otherwise the recorded `CANONICAL_SMILES` field (including rows 58 and 83); record the field choice. Never use this data-field selection as a chemical error fallback. Include both carbonate IDs when matched.

Finite non-nitrate variants remain: H2O/OH−; NH3/NH4+; hydrogen-halide/halide-ion F/Cl/Br/I; archived singleton metals and explicit Li+/Na+/K+/Ag+, Mg2+/Ca2+/Sr2+/Ba2+/Zn2+, Al3+/Bi3+ variants; terminal-O proton states of source phosphate/sulfate/carbonate/hypophosphite with validated valence and preserved P-bound H; raw PF6 and PF6−; source BF4−. Other source rows match their stored exact graph/H/charge state. Enumerate/review terminal-O variants once during preparation, not by broad runtime uncharging. `[S-2]` does not imply HS−/H2S, and atomic halogen radicals do not match hydrogen-halide/halide templates.

**Nitrate row 83: exactly these three templates; no generated proton ladder.** Each has one N plus three terminal O, exactly one N=O and two N–O bonds, no extra atoms/bonds, no isotope labels/radicals/aromatic/query state, no N-bound H, and explicit terminal-O H counts. O permutations are exact graph symmetries, not input-ID choices.

| Template ID | Exact representative graph | N charge; single-bond O `(q,H)`; double-bond O `(q,H)` | Formula/net charge | Status |
|---|---|---|---|---|
| `83_raw_h2no3` | `ON(=O)O` | N=0; `(0,1),(0,1)`; `(0,0)` | H2NO3 / 0 | Verbatim stored source depiction. Recognized **exclusion-only** strict-N-valence exception; never a retained valid nitrate product. |
| `83_nitric_acid` | `[OH][N+](=O)[O-]` | N=+1; `(0,1),(-1,0)`; `(0,0)` | HNO3 / 0 | Adopted comparison-only nitric-acid template; requires ordinary owner validation. |
| `83_nitrate` | `[O-][N+](=O)[O-]` | N=+1; `(-1,0),(-1,0)`; `(0,0)` | NO3 / −1 | Adopted comparison-only nitrate-ion template; requires ordinary owner validation. |

The derived templates are **not** asserted to result from O-only proton edits of the raw source. Record catalogue adaptation `nitrate_comparison_charge_adaptation_v1`: source→acid changes N charge by +1, one single-O charge by −1 and H count by −1, net Δcharge 0; source→nitrate changes N by +1, both single-O charges by −1 and H count by −2, net Δcharge −1. These are source-to-template metadata, **not edits/exchanges applied to a component**. Each row-83 exclusion stores `template_id`, source ID/hash/field, adaptation ID (none for raw), `comparison_only=true`, actual matched formula/charge and original atom mapping. Do not emit a false chemical transformation/proton-reservoir event, modify input formula or repair the component to force a match.

Exclude each of the three exact depictions standalone and in mixtures; they map to source row 83 with distinct template evidence. Do not add protonated `[OH][N+](=O)[OH]` or neutral-N anionic proton variants. Required negatives: labelled nitrate/nitric acid/raw source; nitrite `O=N[O-]`; organic nitrate `CO[N+](=O)[O-]`; nitro `C[N+](=O)[O-]`; unmatched neutral-N variants. Invalid unmatched neutral-N depictions fail `InvalidValence`, not exclusion/unsupported-by-convenience. Labelled variants cannot be erased to match; any later capability failure is reported honestly.

Source PF6/nitrate exact raw depictions and archived unusual metal entries can be comparison-only exclusions after syntax/model topology validation and before strict **retained** valence/domain validation. Recognized nitro/disconnection repairs likewise precede their bounded strict checks. No unmatched invalid graph can survive as a product. Standalone matcher tests cover every source row even if pipeline disconnection first splits a hydroxide into two excluded origins.

### Types, facade, identity and finalization

- Detached entrypoints: `standardize(&TopologyBlock, &MoleculeProperties, &MsReadyParams) -> Result<MsReadyResult, MsReadyFailure>` and `standardize_smiles(&str, &MsReadyParams) -> Result<MsReadyResult, MsReadyFailure>`. Reuse owner SMILES parser; no replacement parser/string edits.
- Parameters: optional record ID/original-SMILES archival metadata and the two explicit limits; chemistry frozen by profile. Text entrypoint preserves actual source text separately from optional metadata.
- Public vocabulary: `MsReadyParams`, `MsReadyResult`, `MsReadyRunInfo`, `MsReadyProduct`, `MsReadyOrigin`, `MsReadyExclusion`, `MsReadyEvent`, `MsReadyWarning`, `MsReadyStatus::{Products,NoEligibleComponent}`, `MsReadyFailure`, structured kind/stage/context.
- Preserve original topology/properties and input linkage. Provenance uses original input IDs, not renumbered fragment IDs: original/post-disconnection component indices, product→input mappings, H dispositions, all origins and multiplicities. Every input atom belongs to a retained/excluded/logged disposition.
- Failure retains original input, stage/rule/component/atom context, underlying category and completed exclusions/events. No partially successful product collection is a success. Distinguish InvalidInput, InvalidValence, UnsupportedChemistry, UnsupportedIdentifier, TransformationLimit, TransformationCycle, InvariantViolation and owner errors. Supported-input failures do not become newly unsupported inputs.
- Products own normalized detached topology, canonical SMILES, required Standard InChI/full 27-character key, isotope-aware composition/formula, formal charge, neutral-atom formula mass and origins. Standard defaults `options=None`; retain engine version/diagnostics. No silently absent key or SMILES fallback. Deduplicate by `(full Standard InChIKey, exact normalized canonical SMILES)` with pinned profile/writer, sort by that key and retain all origins/counts. Exact representation protects against identifier collapse of metal/mobile-H/radical distinctions; no cross-toolkit SMILES-byte promise.
- Early real identifier proof covers ordinary molecules, permanent ions, C/D/T, radicals, **all three retained Hg/Sn tuples**, labelled accepted metal/alkali-ion states. If a required independent identifier capability is absent, stop with evidence before broadening another owner's port or disabling identity. An ordinary supported-input bug must be fixed within scope rather than mislabeled unsupported.
- Under **opt-in** `cap-msready`, thin shared-borrow facade methods `Molecule::with_msready()` and `Molecule::with_msready_with_params(&MsReadyParams)` return detached results. Re-export all public vocabulary through `cosmolkit`. No live product `Molecule`, mutating runtime path or new operation registry.
- Canonical binding registration uses compiler-checked Native entries, shared receiver, `state: value_returning`, `operation: none`, owning feature `cap-msready`, registered nested types/errors; Python names remain snake_case, JS `withMsready`/`withMsreadyWithParams`. Actual Python/WASM bindings are outside this Rust implementation scope and must not be claimed delivered by declared projections. Never hand-edit stubs.
- Preserve existing `default`, `full`, `core` and all nine bundle memberships verbatim. Neither default nor `full` alone exposes MSReady; `full,cap-msready` is the explicit broad-build combination. The leaf enables required implementation dependencies directly, not `cap-smiles`, `cap-inchi`, `cap-descriptors` or unrelated public capabilities. Add the leaf to the existing feature CI loop without weakening existing gates.
- Physically clear atom/bond/enhanced/axial stereo, directions/groups and perception carriers on derived state only. Use empty derived coordinates for serialization/identity; retain original coordinates/properties. Writer options: isomeric=true (isotopes preserved), canonical=true, kekule=false, clean_stereo=false, ignore_atom_map_numbers=true, no root/dative output, ordinary implicit-H output. Test no `@`, `/`, `\\` or stereo identifier layers.
- Reuse selected owner H/valence/aromaticity operations, not generic cleanup/organometallic cleanup. Formula uses `molecular_formula_with_valence(topology, true, false, Some(assignment))`; retain labelled graphs/composition and distinguish raw exclusion metadata from validated retained-product formula.
- Add `cosmolkit-descriptors::neutral_atom_exact_mass_with_valence` as a narrow detached helper using its existing exact-mass kernel and owner-private `RDKIT_ELECTRON_MASS` to undo existing net-charge electron correction. Reuse prepared valence/checked charge; do not copy mass tables/constants into MSReady or substitute average molecular weight. Preserve existing ionic exact-mass semantics. No adduct/m/z feature is added.

## Tasks

1. **Freeze the contract and complete source inventory**
   - Read Step: `AGENTS.md`, `dev/README.md`, `dev/agent_plan_standard.md`, architecture/public API/fixture rules, full reference, archive README/manifest/inventory, selected seven MRV files and nitro settings.
   - Files: new `dev/plans/msready_pipeline_plan.md`, `crates/cosmolkit-msready/README.md`, `testdata/msready/rule_classification.json`.
   - Changes: Persist this approved sequential ledger, frozen predicates/nitrate templates/opt-in feature policy, source paths/indices/hashes, mapped edits/invariants/adaptations and omissions for all 115 appearances. Retain archive-attested licensing context.
   - Acceptance: 32 source exclusions; 3/1 adopted/unadopted disconnections; 6/105 normalization appearances; five group rules/two proton classes; no historical/QSAR claim or existing-bundle change.

2. **Prepare immutable fixtures and local preparation tests**
   - Read Step: archive manifest/exclusions, selected MRV/nitro XML, native sample `.smi`/discarded CSV, repository fixture rules and exact nitrate template contract.
   - Files: new `tools/testdata/prepare_msready.py`, `tools/testdata/test_prepare_msready.py`, `testdata/msready/{README.md,manifest.json,fixtures/}`.
   - Changes: Deterministic local-path-only preparation validates hashes and copies 32-row table/inventory/seven selected MRV/nitro settings/archive README; selects small native sample records with IDs/offsets/selection hashes. Preserve original bytes; no bulk workflow/XLS/SDF/QSAR copy. Enumerate exact reviewed variant data including the three nitrate records and adaptation metadata; never regenerate expected fixtures during tests. Test missing/corrupt source, variant cardinalities/adaptation provenance, deterministic manifest and fixture immutability.
   - Acceptance: Each fixture has version/source path/URL/hash/license/selection identity. No production/build dependence on external absolute paths. Raw row 83 remains unchanged H2NO3.

3. **Run preparation tests immediately**
   - Read Step: required policy/source rules and preparation test entrypoint.
   - Command: `python3 -m unittest discover -s tools/testdata -p test_prepare_msready.py -v`.
   - Acceptance: Record exit/count/log; missing/corrupt inputs fail rather than skip; nitrate table is exactly three templates; no network/fixture refresh.

4. **Implement package, vocabulary and complete input validation**
   - Read Step: model atom/bond/topology/properties; core selected valence/sanitize; detached SMILES parser/parameters; frozen repair/exclusion ordering and metal predicates.
   - Files: `Cargo.toml`, `Cargo.lock`; new `crates/cosmolkit-msready/Cargo.toml`, `src/{lib,types,profile,validate}.rs`.
   - Changes: Rust 2024 workspace member, existing version/license, downward model/core/search/SMILES/descriptors/InChI dependencies. Implement result/error/params/profile types/original snapshot/limit checks and finite graph capability validation. Reject empty/reaction syntax, wildcard/query/polymer/unspecified attachments, malformed topology and unsupported bonds. Separate syntax/model checks from retained-domain checks so exact exclusion-only states and exact repair templates can reach their bounded stages; no global strict-valence bypass. Encode the precise Hg/Sn/alkali predicates, selected-owner validation and invalid-versus-unmodeled categories. Add private `validation_tests`/`profile_tests`, including metal positive/negative graph predicates, repair deferral limits, source raw nitrate exception isolation and feature/profile metadata. Do not expose stub pipeline functions.
   - Acceptance: Originals immutable; errors contextual; predicate-accepted owner failures not relabeled unsupported; Hg unrestricted owner valence does not admit additional states; metadata contains hashes/versions and nitrate adaptation IDs.

5. **Run input/vocabulary checks immediately**
   - Read Step: policy/source rules and validation/profile tests.
   - Commands: `cargo check -p cosmolkit-msready`; `cargo test -p cosmolkit-msready --release --lib validation_tests`; `cargo test -p cosmolkit-msready --release --lib profile_tests`.
   - Acceptance: Nonzero tests/counts/logs; negative graph predicates and narrow valence exceptions execute without cfg-hidden coverage.

6. **Implement the real Standard InChI bridge and early proof**
   - Read Step: InChI `lib.rs:210-283`, `source/rdkit/inchi.rs:95-211,344-385,854-930`; core valence/kekulization/periodic-table APIs; accepted metal/ion tuples.
   - File: new `crates/cosmolkit-msready/src/identify.rs`.
   - Changes: Private transient TopologyBlock→InchiMolecule bridge preserves H/isotopes/radicals/charge/bonds without maps/stereo/coordinates. Implement all eight toolkit methods using owner algorithms: cache observation/update, kekulization, symbol/atomic weight, total H, implicit valence and degree; synchronize adapter changes, no no-op/hardcoded-zero required method. Generate once, derive full key using `inchi_to_inchi_key`, check Standard prefix and keep diagnostics. Add `identity_tests` for acetone diagnostic key, distinct unnormalized enol, amine/quaternary ammonium, C/D/T, methyl radical, `C[Hg+]`, `C[Hg]C`, `C[Sn](C)(C)C`, isotope-labelled accepted metal and q+1 singleton alkali states.
   - Acceptance: Actual engine, not scripted doubles; graph metal bonds remain in product topology even if identifier representation collapses them. All frozen metal-state identity capabilities are proved before further implementation.

7. **Run the blocking identity proof immediately**
   - Read Step: policy/source rules and identity test cases.
   - Command: `cargo test -p cosmolkit-msready --release --lib identity_tests`.
   - Acceptance: Every proof case executes with keys/version/diagnostics. Stop before Task 8 if required engine support is absent; do not silently drop identity, expand owner scope or narrow the accepted domain.

8. **Implement all exact exclusions and explicit nitrate adaptations**
   - Read Step: all 32 raw rows, variant/nitrate contract, core H-mapping and search query/target/whole-graph matching APIs.
   - Files: new `src/exclusions.rs`, `src/data/{exclusions,exclusion_variants}.json`.
   - Changes: Package immutable production-owned data, independent of testdata/archive paths. Exact matching checks cardinality/bijection/all graph attributes, preserves analysis graph and duplicate source IDs. Store exact three nitrate templates with no generated proton ladder; record template/source/adaptation metadata in exclusions rather than fake chemical events. Raw nitrate/PF6 exceptions remain comparison-only. Add `exclusion_tests` for every raw row/reviewed variant, absent source fields, labelled negatives, constitutional isomers/CCl/connected phosphate/sulfate, each nitrate standalone/mixed/permuted/explicit-H form and all stated nitrate negatives. Verify input charge/H/formula/atoms unchanged by comparison and no source-to-template delta mistaken for actual edits.
   - Acceptance: 32 source rows covered; raw row 83 and both valid adapted depictions explicitly distinguished; unsupported/invalid unmatched graphs do not receive false exclusions.

9. **Run exclusion tests immediately**
   - Read Step: policy/source rules and exclusion test entrypoint.
   - Command: `cargo test -p cosmolkit-msready --release --lib exclusion_tests`.
   - Acceptance: Report source-row/variant/negative counts, including exactly three nitrate templates; zero matches/ignored cases not passes.

10. **Implement guarded neutral/charged alkali disconnection and mapped components**
    - Read Step: three raw K/Li/Na MRV sources, frozen charge-pair/ligand predicates, core root-exported `connected_components`/`get_molecule_fragments`, model edits/mappings.
    - Files: new `src/{disconnect,components}.rs`.
    - Changes: Implement exactly `(0,0)→(+1,-1)` and `(+1,-1)→(+1,-1)` for the two permitted ligand forms. Require atom/H/isotope/formula/net-charge conservation and strict post-edit owner validation before proton edits. Narrowly defer bonded charged-case pre-valence, preserving its stated H count; no broad metal workaround. Reuse core splitting/materialization with `sanitize_fragments=false`, `copy_conformers=false`, graph-only coordinates; compose original/component/H mappings. Add `component_tests` for every positive/negative Li/Na/K case in the contract; accepted Hg/Sn preservation/negative tuples; mixtures/repeats/order/accounting. Direct exclusion matching still covers raw source hydroxide entries; pipeline splitting records both resulting exclusions and all origins.
    - Acceptance: Ionic/bonded acetate convergence for both permitted pairs; rejected pair never changes charge/H; no generic metal deletion; preserved metal bonds and later-excluded ions traceable.

11. **Run component/disconnection tests immediately**
    - Read Step: policy/source rules and component tests.
    - Command: `cargo test -p cosmolkit-msready --release --lib component_tests`.
    - Acceptance: All three metals, both charge pairs, both ligand forms, negative pair/coordination/H cases and isotope/accounting assertions execute.

12. **Implement all five directed rules and the bounded deterministic engine**
    - Read Step: selected MRV reactant/product/R/ligandOrder/bond data, both nitro raw rows, source protocol and model/search editing.
    - Files: new `src/{rules,rewrite}.rs`, `src/data/rules.json`.
    - Changes: Encode reviewed project queries/correspondences/guards/edits/invariants, not preview SMILES or a new generic MRV framework. Preserve branches/stable IDs; fresh matching, deterministic overlap selection, six appearances, bounded cycle/state control and multi-nitro repair before strict validation. Add `rule_tests` for every rule/guard, implicit/bracket/explicit/D/T transfer, adjacent/overlap/order-sensitive motifs, asymmetric branches/nitro isotopes, aromatic controls and atom permutations.
    - Acceptance: Complete five-rule implementation, no lost substituent; per-edit formula/net-charge/isotope/map checks; independent source fidelity/complexity assessment.

13. **Run directed-rule tests immediately**
    - Read Step: policy/source rules and rule tests.
    - Command: `cargo test -p cosmolkit-msready --release --lib rule_tests`.
    - Acceptance: Every rule/both nitro appearances and negative guards execute; invariant assertions active.

14. **Implement guarded charge/H bookkeeping and radical retention**
    - Read Step: reference §7, frozen two proton classes, core H/valence and model isotope/radical carriers.
    - Files: new `src/{charge,hydrogens}.rs`.
    - Changes: N deprotonation/O protonation, balanced internal zwitterion transfer and residual reservoir events; composed fold/deletion maps; preserve D/T/necessary charges/radical electrons with specific warnings. Validate valence and H/charge delta. Add `charge_tests` for primary/secondary/tertiary versus permanent/aromatic ammonium, carboxylate/alkoxide/glycine/multiple sites, nitro/azide, explicit H/D/T and radicals. Nitrate exclusion-template metadata must never enter this proton engine as an actual source-to-product rewrite.
    - Acceptance: Internal transfer preserves isotope composition/net charge; reservoir H and q change together; unedited valid charges are not unsupported.

15. **Run proton/H tests immediately**
    - Read Step: policy/source rules and charge tests.
    - Command: `cargo test -p cosmolkit-msready --release --lib charge_tests`.
    - Acceptance: Graph/H/composition/charge assertions, not string-only or total-charge-only tests.

16. **Implement finalization and owner-backed neutral-atom mass**
    - Read Step: selected core H/sanitize/aromaticity, model stereo/groups/properties, writer parameters, descriptor exact mass/formula/private electron constant.
    - Files: new `src/{finalize,composition}.rs`; `crates/cosmolkit-descriptors/src/lib.rs`.
    - Changes: Clear all derived atom/bond/enhanced/axial/perception stereo/maps/drawing carriers, use empty derived coordinates and pinned isotope-preserving writer; ordinary-H folding retains D/T. Add narrow descriptor-owned `neutral_atom_exact_mass_with_valence` using existing kernel/private constant/checked charge/prepared assignment; no duplicated tables or changes to existing ionic mass. Add MSReady `finalize_tests` and owner-private `neutral_atom_exact_mass` tests: neutral/cation/anion/multiple-charge/isotope/implicit-H masses with tolerance and unchanged ionic outputs; tetrahedral/EZ/enhanced/axial/coordinate-induced stereo, D/T, aromatic/Kekule forms, retained metal connectivity and identifier/writer no-stereo checks.
    - Acceptance: Correct labelled neutral-atom formula mass and isotope-aware formula; originals preserved; no isotope suppression to remove stereo.

17. **Run finalization and mass-owner regressions immediately**
    - Read Step: policy/source rules and both changed-owner test modules.
    - Commands: `cargo test -p cosmolkit-msready --release --lib finalize_tests`; `cargo test -p cosmolkit-descriptors --release neutral_atom_exact_mass`; `cargo test -p cosmolkit-descriptors --release`.
    - Acceptance: Existing ionic regressions retained; owner/molecule graph/isotope/stereo and identifier checks pass, with actual counts/logs.

18. **Implement complete pipeline, deduplication and acceptance/invariant tests**
    - Read Step: reference §§11,13,14; final stage interfaces and frozen provenance/error/identity/domain/nitrate contract.
    - Files: new `src/pipeline.rs`, complete `src/lib.rs`, new `tests/{acceptance,invariants,source_data,release_regressions}.rs`.
    - Changes: Compose structural validation→selected disconnection→components/exact exclusions→strict retained validation/recognized repair→group/proton fixed point→isotope/H/stereo/aromatic finalization→validation/final component filters/accounting→identity/composition→dedup. Honor exact bounded pre-validation repair/exclusion exceptions, not arbitrary permissive graphs. Re-split/re-filter if connectivity changes. NoEligibleComponent only for real exclusions; other failures retain context rather than partial success. Add all 20 §13 worked cases, accepted/rejected metal predicates, three nitrate templates/negative chemistry, idempotence/atom+component permutations/multiplicity/original immutability/limits/source hashes. Native release samples are small fixed regressions with documented profile differences, not historical parity; CSV InChI is not authoritative.
    - Acceptance: Organic counterions/all eligible components survive, duplicates retain every origin; nitrate adaptations are comparison-only; retained metals/isotopes/radicals have required identities and mappings; no silent omission or supported-input relabeling.

19. **Run isolated-crate completion immediately**
    - Read Step: policy/source rules and crate unit/integration test inventory.
    - Commands: `cargo test -p cosmolkit-msready --release`; `cargo check -p cosmolkit-msready`.
    - Acceptance: Nonzero complete stage/integration coverage, exact counts/limitations; Native profile completion, not historical parity.

20. **Integrate opt-in thin facade and canonical registration**
    - Read Step: public API §7, Cargo dependency/feature/binding patterns, return classifications and existing feature-selection tests.
    - Files: `crates/cosmolkit/{Cargo.toml,src/lib.rs,src/binding_contract/registry.rs,tests/feature_selection.rs}`, new `src/msready.rs`, `tests/msready.rs`, root lockfile as needed.
    - Changes: Optional downward dependency and independent **opt-in** `cap-msready`, preserving `default`/`full`/every existing bundle. Re-export public types and two shared-borrow forwarders; register actual Native nested types/methods before exposure. No runtime mutation/operation matrices/live-product construction. Facade tests focus on registration/signature/default/error/sharing/export/leaf availability and dependency isolation, not duplicated chemistry suites. Add explicit compiler availability tests: no-default missing leaf; default/full missing leaf; leaf alone present without unrelated APIs; full+leaf present. Extend existing feature-selection/compile-case mechanisms, not a production handwritten bundle registry.
    - Acceptance: `full` membership unchanged, default/full alone does not expose MSReady, leaf enables only required implementation dependencies, original blocks shared/unchanged on success/failure.

21. **Run facade, opt-in feature, strict and privacy gates immediately**
    - Read Step: policy/source rules, facade/feature/privacy entrypoints.
    - Commands: `cargo check -p cosmolkit --no-default-features`; `cargo check -p cosmolkit --no-default-features --features full`; `cargo check -p cosmolkit --no-default-features --features cap-msready`; `cargo check -p cosmolkit --no-default-features --features full,cap-msready`; `cargo test -p cosmolkit --release --no-default-features --features cap-msready,cap-smiles,op-contracts-strict --test msready`; `cargo test -p cosmolkit --release --no-default-features --features cap-msready,op-contracts-strict --test feature_selection`; `cargo test -p cosmolkit --release --no-default-features --features full,op-contracts-strict --test feature_selection`; `cargo test -p cosmolkit --release --test migration_run_privacy`.
    - Acceptance: Absence/presence compile tests genuinely execute; no hidden cfg coverage or weakened privacy cases. Live-product API changes require a new approved generated multi-output contract and revised plan.

22. **Document exact delivery and extend existing CI feature loop**
    - Read Step: `.github/workflows/features.yml`, feature/example/README patterns and final source/support catalogue.
    - Files: `.github/workflows/features.yml`, `crates/cosmolkit/README.md`, `crates/cosmolkit-msready/README.md`, new `crates/cosmolkit/examples/msready.rs`, its gated Cargo entry.
    - Changes: Add opt-in leaf to the existing loop without replacing gates. Document finite rules/order/variants (especially three nitrate templates/adaptation-only semantics)/metal predicates/limits/identity/mass/accounting; explicitly explain that `full` and defaults do not include MSReady and show `--features cap-msready` or `full,cap-msready`. Add facade-only counterion/multiplicity example requiring `cap-smiles,cap-msready`. Mark Python/JS projections declared, not delivered. Scratch/PR notes stay outside worktree.
    - Acceptance: No normative-doc/bundle edits; examples use only `cosmolkit`; production build/runtime independent of external fixtures; limitations visible.

23. **Run example, formatting and all mandatory final/pre-commit gates**
    - Read Step: current `dev/README.md`/AGENTS validation requirements and every `run` step of the updated features.yml feature-matrix.
    - Commands: `cargo run -p cosmolkit --no-default-features --features cap-smiles,cap-msready --example msready`; `cargo fmt --all`; `cargo fmt --all --check`; `cargo check -p cosmolkit --features cap-msready,op-contracts-strict`; `cargo test -p cosmolkit --release --features cap-msready,op-contracts-strict`; `cargo test --workspace --release --features cosmolkit/op-contracts-strict,cosmolkit/cap-msready`; all existing/new feature-matrix run steps; `UV_OFFLINE=1 uv sync --locked --group dev`; `cargo test --workspace --locked --profile dev-test --no-fail-fast --exclude cosmolkit-parity-tests --exclude cosmolkit-parity-tests-fixed --features cosmolkit/op-contracts-strict,cosmolkit/cap-msready`; `cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen`; `.venv/bin/maturin develop --profile dev-test --manifest-path python/Cargo.toml`; `.venv/bin/pytest python/tests`.
    - Acceptance: Exact state/head/exit/count/log in managed artifacts. Release workspace remains **unexcluded**; README exclusions apply only to its dev-test gate. New opt-in leaf explicitly enabled on changed-behavior gates, not assumed from full. Missing caches/fixtures, zero matches/ignored suites or unrun checks are not passes. No hand-edited stubs. If in-scope owner fixes touch core, additionally `cargo check -p cosmolkit-core` and `cargo test -p cosmolkit-core --release`; affected InChI owner tests if touched. Fix in-scope failures/repeat affected gates; genuine offline/scope blockers require a decision, not network/weakened checks.

24. **Obtain independent read-only acceptance review**
    - Read Step: final diff, source/catalogue hashes, completed ledger and exact validation evidence.
    - Artifact: managed reviewer output outside worktree.
    - Acceptance: Review source/query/edit constraints and complexity, positive/negative metal predicates, comparison-only nitrate metadata, identity/H/mass correctness, all-origin accounting, immutability, license evidence, unchanged bundle membership/opt-in feature, scope and AGENTS exclusion. Parent dispositions control fixes; revalidate affected gates. No commit with blockers or mandatory failures.

25. **Prepare local PR notes and explicit staging allowlist**
    - Read Step: delivered support/API docs and accepted review/validation evidence.
    - Artifact: managed local PR description outside Git.
    - Changes: Include base/branch, Native rule support/adaptations, nitrate template policy, exact metal predicates, opt-in API, commands/results, explicit unsupported capabilities, terms and lack of historical equivalence. Allowlist only intentional code/data/docs/tests/lockfile. Exclude AGENTS, original untracked reference/ZIP/unrelated changes, runtime/log/corpus artifacts.
    - Acceptance: Paste-ready honest notes without remote PR; no invented engine run, pass or delivered language binding.

26. **Create authorized local Conventional Commit**
    - Read Step: allowlist and final gate/review evidence.
    - Commands: stage individual allowlisted files (never `git add -A`, `git add .` or unrestricted directory); `git diff --cached --name-only`; `git diff --cached --check`; assert `git diff --cached -- AGENTS.md` empty; commit `feat(msready): add audited published-concept pipeline`.
    - Acceptance: Reviewed intentional implementation only; no AGENTS exception, normative-policy change, unrelated edits, bulk archive/corpus/logs/version bump or remote action.

27. **Verify local delivery and preserve both checkouts**
    - Read Step: local commit/branch/status and original-checkout status.
    - Commands: `git rev-parse HEAD`; `git status --short --branch`; `git diff d39640716a6c9cf007350288349a0604fde0936a..HEAD -- AGENTS.md`; `git merge-base HEAD d39640716a6c9cf007350288349a0604fde0936a`; local original-checkout status.
    - Acceptance: Correct base/ancestry; original changes intact; AGENTS exception unstaged/local and absent from committed diff. Report commit/worktree/branch/base/PR-notes path/results/risks; preserve worktree, no push/remote PR/status.

## Sequential execution ledger

Task descriptions are the acceptance map; this is the sole execution checklist to persist at Task 1 after approval. No step is already completed by reconnaissance or this planning revision.

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Freeze the corrected contract and source classification specified in Task 1.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Implement immutable fixture preparation and its tests specified in Task 2.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Run the immediate preparation gate specified in Task 3.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Implement package/vocabulary/input and exact-domain validation specified in Task 4.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Run the immediate validation/profile gates specified in Task 5.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Implement the real-engine identity bridge and expanded metal-state proof specified in Task 6.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run the blocking real-identity proof specified in Task 7.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Implemented the versioned finite variant catalogue and complete exclusion matcher in `src/exclusions.rs` and `src/data_exclusion_variants.json`. The matcher covers all 32 raw rows, finite variants and duplicate carbonate IDs, preserves source fields/mappings, delegates ordinary-H folding to core mappings, and records exact three-template nitrate/PF6 comparison metadata.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Exclusion gate passed: `CARGO_NET_OFFLINE=true UV_OFFLINE=1 PIP_NO_INDEX=1 cargo test --offline -p cosmolkit-msready --release --lib exclusion_tests` (5 passed, 0 failed, 10 filtered). Follow-up `cargo test --offline -p cosmolkit-msready --release` passed all 15 unit tests; `cargo check --offline -p cosmolkit-msready` passed. Evidence: `subagent-artifacts/outputs/c70ee8d9-8edf-4648-b49d-c2c25ccfaae3/exclusions-handoff.md`; prior parent acceptance was provisional pending independent review.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Implemented bounded neutral/charged Li/Na/K–oxygen disconnection, full frozen pair/degree/ligand/bond/H/radical checks, cloned input edits, formula/charge/per-atom H/isotope/radical and retained atom/bond invariants, strict post-edit component validation, and original atom/bond mapping composition through disconnection and core component splitting. Regression coverage includes both charge pairs, organic and hydroxide ligand classes, unsupported connectivity/coordination/H/radical/isotope cases, Hg/Sn domain retention, exclusion origins, repeated/permuted components and map accounting.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Task 11 acceptance gate passed offline in release: `cargo test --offline -p cosmolkit-msready --release --lib component_tests` (7 passed), `component_mapping_tests` (2 passed), `validation_tests` (5 passed), `exclusion_tests` (5 passed), `identity_tests` (1 passed), and the complete `cargo test --offline -p cosmolkit-msready --release` suite (21 passed, 0 failed; doc-tests 0). `cargo check --offline -p cosmolkit-msready`, `cargo fmt --all`, and `cargo fmt --all --check` passed.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Added the six-appearance versioned rule catalogue and five source-anchored match/edit rules in `src/{rules,rewrite}.rs`. Fresh matching, lexicographic isotope-preserving canonical selection, repeated fixed-point sweeps, bounded candidate/edit limits with typed errors, canonical cycle detection, original atom/bond map composition, explicit/bracket/implicit/tracked-D/T transfer, and per-edit formula/charge/isotope/H/radical/row-state invariants are implemented. Regression tests cover every appearance, all five edits, isotope-asymmetric nitro and atom reordering, multiple nitro centers, branch/row retention, aromatic negatives, hydrogen carriers, D/T, and both typed limits. Source execution behavior/complexity remains conservatively marked unresolved against static MRV assets; independent review remains required.
Step 26 [x]: Re-read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before the immediate directed-rule gate.
Step 27 [x]: Directed-rule gate passed offline: `cargo test --offline -p cosmolkit-msready --release --lib rule_tests -- --nocapture` (11 passed, 0 failed, 21 filtered); complete `cargo test --offline -p cosmolkit-msready --release` (32 passed, 0 failed, 0 doc-tests); complete debug library suite (32 passed, 0 failed); `cargo fmt --all -- --check` and `cargo check --offline -p cosmolkit-msready` passed. Existing upstream warnings remain. Independent source-fidelity/complexity review is pending.
Step 28 [x]: Re-read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before Task 14.
Step 29 [x]: Implemented guarded single-bond aliphatic N+ deprotonation, eligible terminal O− protonation, same-component balanced transfer before residual reservoir exchanges, explicit/bracket/inferred/tracked-isotope H carriers, source-row atom/bond mapping composition, per-input-H folded/retained/reservoir-deletion dispositions, post-edit owner valence recomputation, formula/charge/isotope/radical invariants, and stable-radical/blocked-D/T warnings. Nitro/azide, permanent/aromatic ammonium and nitrate comparison templates remain unedited. Added 11 `charge_tests`, including multiple centers, mapped explicit-H folding/deletion, explicit and tracked D/T, and formula/charge/center-state assertions. Source anchors use the pinned reference text; this is not a claim of historical Indigo parity.
Step 30 [x]: Re-read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before the immediate proton/H gate.
Step 31 [x]: Offline release gate passed: `cargo test --offline -p cosmolkit-msready --release --lib charge_tests -- --nocapture` (11 passed, 0 failed, 32 filtered); complete `cargo test --offline -p cosmolkit-msready --release` (43 passed, 0 failed, 0 doc-tests; all prior validation/identity/exclusion/component/rule regressions included); `cargo fmt --all -- --check` and `cargo check --offline -p cosmolkit-msready` passed. Existing upstream/dead-code warnings remain.
### Tasks 1–15 independent chemistry-review remediation (before Step 32)

The historical focused passes above did not establish the missing chemistry contracts identified by independent review F1–F11. Parent-authorized remediation now implements H-filled R slots, canonical already-charged nitro roles/no-op handling, transient tracked-D/T identity/ranking materialization, complete internal carrier enumeration, finite retained-metal rejection, shared exact raw-nitrate exception predicates, bounded repair-batch validation, cumulative group/proton edit and candidate budgets, and the distal ynol state guard. Sixteen additional private behavioral regressions retain all prior coverage; the set-mechanics test is separately named, and a real checked-graph cycle test exercises the actual private engine boundary. Disconnection now has verbatim pinned MRV/hash anchors with separate Native behavior/complexity assessments.

A live supervisor decision additionally confirmed the already-frozen neutral-nitro repair ordering: strict disconnection validation may validate a transient copy repaired ONLY by the adopted guarded neutral-nitro matcher/edit. The analysis graph remains neutral/unmodified until the later actual group stage logs its repair. Standalone/no-cut, mixed neutral/charged alkali, multiple nitro, exact-exclusion neighbors and unmatched-invalid controls are covered. This is not a generic valence exemption. Integration must share one `rewrite::TransformationBudget` across `_with_budget` group/proton calls for each original input; standalone wrappers intentionally create their own budget.

Final offline release gates: `cargo test --offline -p cosmolkit-msready --release --lib` and `cargo test --offline -p cosmolkit-msready --release` both exit 0, **59 passed, 0 failed/ignored/filtered** (full crate has 0 doc-tests). Focused final validation/component/rule/charge/identity gates pass **7/10/16/15/3** tests respectively; `cargo check --offline -p cosmolkit-msready`, `cargo fmt --all` and `cargo fmt --all --check` pass. Source audit verifies all 9 implementing-function source blocks against pinned fixture hashes/verbatim excerpts. All failed pre-fix regressions and intervening test-harness corrections remain in durable logs. No staging/network/owner-mass-helper/finalization/pipeline/facade changes were made. **READY_FOR_REVIEW** applies only to F1–F11 remediation; independent acceptance review and Step 32 onward remain pending. Evidence: `/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/ddffa8cd-32e2-47c6-a7ea-56e5ea81a371/chemistry-fixes.md` and adjacent `logs/`.

Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Implemented private finalization/composition stages with physical stereo/map/drawing/perception clearing, empty derived coordinates, mapped ordinary-H folding, validated tracked-D/T materialization, selected owner aromaticity/H/valence, pinned canonical isotope-preserving writer and real Standard identity. Reuses the completed descriptor neutral-atom helper unchanged. Twelve private finalization regressions plus one direct real aromatic identity regression pass (12/4 focused tests). Supervisor authorized the narrow shared `identify.rs` aromatic-bond flag transport fix after two real-identifier failures; no lower-owner rewrite, domain narrowing or identity roundtrip. Initial failure and corrected focused logs are retained in the finalization handoff artifact.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Offline release gates passed with `CARGO_NET_OFFLINE=true UV_OFFLINE=1 PIP_NO_INDEX=1`: `cargo test --offline -p cosmolkit-msready --release --lib finalize_tests` (12 passed, 60 filtered); `cargo test --offline -p cosmolkit-descriptors --release neutral_atom_exact_mass` (2 passed); full descriptors release (191 unit + 57 integration + 10 doc tests, 258 total passed); full MSReady release (72 unit passed, 0 doc tests); focused identity (4 passed). All exit 0, no failures/ignored tests. `cargo check --offline -p cosmolkit-msready`, `cargo fmt --all`, `cargo fmt --all -- --check`, verbatim Native-anchor audit and read-only `git diff --check` pass; index inspection empty. Initial two aromatic identifier failures and a 120-second descriptor-build tool timeout are retained, not counted as passes; full successful rerun retained. Finalization handoff and durable logs/counts/exits: `/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/45f0925a-d45a-4f8d-94b5-76139a7c55e5/finalization-handoff.md` and adjacent `logs/`. Steps 36 onward and independent stage review remain pending.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Implemented detached `standardize`/`standardize_smiles`/`run_info`, complete result/failure/origin/event/H/atom-accounting vocabulary, raw snapshot and actual/archival source linkage, selected disconnection and exact filters, strict retained checks with existing bounded nitro deferral, whole-profile group/proton fixed point, one input-wide budget (including full disconnection candidate attempts/committed cuts), private isotope/stereo finalization, real Standard identity/composition, sorted full-key/exact-SMILES deduplication and graph-aligned all-origin mappings. Added public-boundary `acceptance` (6 tests/all 20 §13 cases), `invariants` (10), `source_data` (4/all 32 raw rows and 26 finite variants/immutable hashes), and `release_regressions` (2/all 8 selected Native records); each focused offline release gate passes. One additional private whole-profile real-canonical-state cycle/no-op boundary regression passes. Supervisor approved owner radical initialization on retained text components before guards, narrowly observing only exact pending charged-alkali repairs on a bond-deleted copy; raw/original graphs and explicit radicals remain unchanged. This adaptation is documented in-function and in the managed pipeline handoff, not historical parity or a chemistry-owner rewrite. Initial failures and evidence-based malformed-fixture/expected-edit-count corrections remain in logs; no facade/CI/owner-data/fixture/policy/AGENTS edits or staging.
Step 38 [x]: Re-read `dev/policy_invariants.md` and the complete `dev/source_reproduction_protocol.md`, Task 19's commands/acceptance and the crate unit/integration inventory before the isolated-crate completion gate.
Step 39 [x]: Final Task 19 gates passed offline with `CARGO_NET_OFFLINE=true UV_OFFLINE=1 PIP_NO_INDEX=1`: unfiltered `cargo test --offline -p cosmolkit-msready --release` exit 0, 73 unit + 22 integration = 95 passed, 0 failed/ignored/filtered (0 doc tests); `cargo check --offline -p cosmolkit-msready`, `cargo fmt --all` and `cargo fmt --all -- --check` all exit 0. Full run includes every earlier private stage regression, new actual whole-profile cycle/no-op boundary, and acceptance/invariants/release_regressions/source_data suites (6/10/2/4). Expanded audit verifies 18 pinned/verbatim MRV/XML/reference blocks and the descriptor owner's unchanged Task 17 hash. It exposed one inherited paraphrased charge anchor, corrected to the exact pinned sentence without behavior changes; focused charge tests (15 passed, 58 filtered) and the complete final gates were rerun successfully. Read-only tracked-diff/index checks pass with no staged files. Durable command/exit/count/hash/source-audit/index evidence and the authoritative detached handoff are under `/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/45f0925a-d45a-4f8d-94b5-76139a7c55e5/pipeline-handoff.md` and adjacent `pipeline-logs/`. Native finite-profile detached completion only: no historical KNIME/Indigo parity, facade/bindings/CI or workspace acceptance is claimed. Independent parent-managed acceptance review remains required; Step 40 onward is untouched.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Implement opt-in facade/registrations/features/tests specified in Task 20.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Task21 offline gates passed: all four minimal/full/leaf/full+leaf checks; focused facade (3 passed), feature_selection leaf/full/full+leaf (6/8/8 passed), four-block sharing unit (1 passed), native rustdoc minimal/default/full/leaf absence-isolation probes (59/14/14/55 passed), and unchanged real-module migration_run_privacy (1 passed; eight compiler invocations preserve default/strict compile-pass/fail coverage). All exit 0; exact commands/counts/durations/logs in integration-handoff adjacent gates.jsonl/logs.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Implement exact support/docs/example/CI-loop updates specified in Task 22. Parent follow-up clarified candidate budgets: README and public facade parameter documentation now distinguish input-wide committed edits from separate disconnection/group/proton candidate sweeps, pooling all six group appearances/rematches or all internal/reservoir proton phases/rematches without per-edit resets. Read-only source/doc audit passes; no detached code change or barrier-dependent gate restart.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Task23 task-specific scoped validation disposition is accepted under the user's explicit local-delivery decision recorded below; full workspace/Python/stub validation has NOT passed. Historical failed gate evidence follows unchanged: Task23 was BLOCKED, not complete: example, fmt/write/check and full+leaf strict check pass, but exact `cargo test -p cosmolkit --release --features cap-msready,op-contracts-strict` exits 101 (428 passed / 3 failed / 0 ignored before early termination). All three `conformer_original_fixed19` failures require missing `testdata/conformer/expected/rdkit/smiles_small/conformer_generation.jsonl` at line 21. Test source matches pinned HEAD; the fixture is absent in both isolated and original checkouts. Supervisor authorized inventorying remaining exact mandatory gates without exclusions, fixture generation, network or out-of-scope fixes. Durable commands/exits/counts/logs and final disposition: integration-handoff.md / gates.jsonl adjacent to output e73e97e0-f64d-47ef-ae90-4e395686af1a. This failure is not MSReady readiness or a full-suite pass. During the unexcluded workspace-release attempt, the parent started an authorized nonoverlapping detached F12/F13 fix writer. Per live decision, already-running workspace Cargo completed safely as mixed-state/non-acceptance evidence (exit101, 559.73s, 428 passed / 3 failed / 0 ignored at the same missing fixed19 fixture); The parent lifted the writer barrier after the stable F12/F13 handoff (`outputs/1683071d-8064-4f8e-8666-ec2accfc2497/pipeline-fixes.md`, reported 98 detached release passes; final independent `final-code-review.md` verdict Code OK closes F1–F13, and raw descriptor logs correct the owner total to 191 unit + 57 integration + 10 doctests = 258, not the earlier 253); its six recorded source hashes match and final `cargo fmt --all` changes no Rust bytes. Revival inventory executed all 23 prescribed/focused commands on that stable seam: final focused facade/leaf/sharing/full+leaf/privacy gates pass (3/6/1/8/1), fmt/check/example/strict check and all 17 CI leaf/composed checks pass. Exact full facade and unexcluded workspace release each exit101 with 428 passed / 3 missing-fixed19 failures / 0 ignored before early termination. README-only-excluded dev-test workspace exit101: 9,714 passed / 37 failed / 82 ignored across 338 executed targets; includes THREE current-diff registry prefix/list/Native-status failures, missing-reference prerequisites and one pre-existing hardcoded WASM output-path permission failure. CI release tooling exit1 (10 tests, 2 current-diff release dependency-catalogue failures); CI coverage/minimal/core gates pass. Offline `uv sync --locked --group dev` succeeds (41 cached packages); official stub generator exit101 because `deranged v0.5.8` is uncached, no stub output/delta; fresh prescribed maturin build succeeds; pytest exit2 after collecting 2,322 items with 1 missing-fixed19 collection error (no tests executed). No implementation/checker/expectation/fixture edits or network workaround. Both owned gate controllers completed and disappeared; no SIGSTOPped chain remains. Step47 remains UNCHECKED. Exact command/exit/count/hash/failure-classification/index evidence and required parent/user dispositions: `/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/5a487d48-a082-4fc5-b62d-d1714ad2c535/final-gate-inventory.md` and adjacent logs/JSON.
Task23 authorized integration-correction follow-up: **INTEGRATION_FIXES_COMPLETE for the five current-diff failures only; Step47 remains UNCHECKED.** Added exactly eight approved downward production edges to the canonical `dev/tools/bump_rust_version.py` catalogue without changing versions/dependencies or any legacy checker. Moved the unchanged 26 Native MSReady declarations to the registry suffix, preserving every legacy registry byte/order and shared/compiler-checked contracts. Extended `migration_run_schema` with the exact opt-in list/Native commitment and vocabulary/suffix/absence regression; unchanged `migration_q_parser_api` Element-prefix check passes. Initial focused reproduction: release tooling 8 passed/2 failed (exit1), registry 26 passed/3 failed (exit101). Patched immediate release tooling: all 11 tests pass; final release strict default+leaf registry/feature tests 38 pass, default/no-leaf 38 pass, minimal/no-leaf and minimal/leaf 15 each pass. Focused facade/sharing/macro-binding/privacy pass 3/1/31/1 (sharing retains 108 filtered; privacy executes eight actual compiler invocations); default/minimal no-leaf rustdoc probes each pass four compile-fail cases. Formatting/check, strict/minimal checks and diff check pass. Exact README no-fail-fast dev-test workspace rerun exits101: **9,718 top-level passed / 34 failed / 82 ignored / 0 filtered**, 338 top-level result groups; one nested conformer child result is separately logged (raw totals 9,719/34/82/203 across339 groups), not double-counted as an additional suite pass. Failure-name multiset equals the prior inventory minus precisely the three corrected registry failures: 33 missing-prerequisite and one unchanged WASM-path failure across nine targets. Full strict/release and Python/offline stub blockers from the prior inventory remain unresolved; those known failing gates were not restarted. No chemistry/facade/owner/manifest/lock/AGENTS/standards/checker/WASM changes, fixture preparation, downloads, staging or commits. Parent-managed re-review remains required; no mandatory acceptance/commit-readiness claim. Exact commands/exits/counts/logs, final hashes and preservation evidence: `/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/336fa0be-baea-43dc-8c9b-4e2800620b0a/integration-fixes.md` and adjacent JSON/logs.

### Explicit user scoped-acceptance disposition for Tasks 23–27

The user expressly authorized LOCAL commit/PR preparation after accepting only failures relevant to this change as blocking. Independently established pre-existing missing fixtures, the unchanged WASM output-path permission failure and the offline existing stub dependency cache are **nonblocking for this change**, not passes or supported skips. All current-diff failures (F1–F13 and release-catalogue/registry-prefix/exact-inventory/Native-commitment corrections) are closed by the independent final reviews. The exact failed full gates, ignored/unexecuted coverage and raw logs above remain authoritative; no checker, fixture, normative standard or AGENTS rule was changed or weakened. No network, push or remote PR is authorized. This is a task-specific acceptance disposition, not full-workspace acceptance.

Step 48 [x]: Reread the complete policy/source protocol, final diff and source/catalogue hashes, completed ledger and exact validation evidence before Task 24; all 3,112 final reviewed snapshot hashes match the delivery entry state.
Step 49 [x]: Independent Task 24 read-only code reviews accepted: `outputs/5a487d48-a082-4fc5-b62d-d1714ad2c535/final-code-review.md` closes F1–F13 with Code OK; `outputs/5b76ae69-b3fd-4ba9-bd0d-0c2ab04ea047/integration-final-review.md` accepts the final five-file integration correction with no issues. Their baseline/resource delivery blockers are disposed as nonblocking solely by the explicit user decision above. No new implementation changes are authorized or made; reviewed source hashes match.
Step 50 [x]: Reread the complete policy/source protocol, delivered owner/facade support documentation and accepted review/validation/fixture-handoff evidence before Task 25.
Step 51 [x]: Prepared honest paste-ready local PR notes at `/data/dev/COSMolKit-msready-PR.md` and an explicit 63-file allowlist at `outputs/e53f3483-be26-40a1-bfd3-1f5efed48715/staging-allowlist.txt`. No broad directory staging. Source/owner review snapshot matches; separately reviewed fixture/tooling seam confirmed from prior fixture handoff and fresh preparation (8 passed) / release source_data (4 passed) gates, each exit 0. All 17 selected immutable outputs and manifest compare byte-for-byte to pinned sources/selection in memory without fixture writes. Five immutable originals/selected originals have CRLF-only default Git whitespace diagnostics: native discarded CSV, selected SMI, exclusion CSV and initial/final nitro XML. Exact source hashes and zero real trailing spaces/tabs are verified. Explicit supervisor disposition permits exact-path per-command `cr-at-eol` recognition for those five only, while all other paths retain default checks; no fixture/manifest/checker/config/gitattributes change. Both AGENTS exceptions, generated artifacts and original unrelated inputs remain excluded. Final commit identity/status will be appended to external notes/receipt only after actual execution.
Step 52 [x]: Immediately before Task 26, reread the complete policy/source protocol, exact 63-file allowlist, accepted final review/gate/fixture evidence and updated local PR notes; final staging will reassert branch/base/main, empty index, reviewed source identities, intended status and original/source-byte preservation. Default all-path whitespace diagnostics on the five exact-original CRLF fixtures will be retained; ordinary paths and each of those five exact-path CRLF-aware checks must pass before commit. No commit or committed-head verification is asserted in this source snapshot.
Step 53 [ ]: Create the authorized local Conventional Commit specified in Task 26. Pending at this pre-commit source snapshot; actual outcome is recorded in the external terminal receipt below, never asserted before execution.
Step 54 [ ]: Reread policy/source protocol and actual local commit/branch/status plus original-checkout status immediately before Task 27; pending at this pre-commit source snapshot, with completed read evidence recorded externally.
Step 55 [ ]: Verify and report local delivery/preservation specified in Task 27. Pending at this pre-commit source snapshot. To avoid a post-commit ledger edit or circular commit-hash assertion, the actual completion of Steps 53–55 and verified commit/base/branch/index/worktree/original-preservation outcomes will be appended durably to `/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/e53f3483-be26-40a1-bfd3-1f5efed48715/local-delivery.md`; accompanying logs/JSON and external PR notes retain the exact receipt. Source checkboxes are deliberately not advanced ahead of the Git actions.

## Files to Modify

- `Cargo.toml`, `Cargo.lock` — workspace member/local dependency delta only, no upgrade/version bump.
- `crates/cosmolkit-descriptors/src/lib.rs` — narrow neutral-atom exact-mass helper and owner regressions, unchanged existing ionic mass behavior.
- `crates/cosmolkit/Cargo.toml` — optional downward dependency, **opt-in** `cap-msready`, gated example; unchanged bundle/default memberships.
- `crates/cosmolkit/src/lib.rs` — gated module and public vocabulary exports.
- `crates/cosmolkit/src/binding_contract/registry.rs` — compiler-checked Native type/method contracts with cap-msready gating.
- `crates/cosmolkit/tests/feature_selection.rs` — leaf isolation and absent-under-default/full availability checks via existing mechanisms.
- `crates/cosmolkit/README.md` — delivered Rust support/opt-in API details.
- `.github/workflows/features.yml` — new opt-in leaf in existing loop; retain all prior gates.
- `AGENTS.md` — **no further edits**; existing authorized local-only exception remains forbidden from staging/commit.

## New Files

- `dev/plans/msready_pipeline_plan.md` — approved dedicated sequential ledger, authorized only for isolated MSReady task.
- `crates/cosmolkit-msready/{Cargo.toml,README.md}` — detached package and complete finite Native contract.
- `crates/cosmolkit-msready/src/{lib,types,profile,validate,identify,exclusions,disconnect,components,rules,rewrite,charge,hydrogens,finalize,composition,pipeline}.rs` — narrowly owned stages/private stage tests.
- `crates/cosmolkit-msready/src/data/{exclusions,exclusion_variants,rules}.json` — package-owned production data including exact nitrate templates/adaptation IDs, no runtime fixture dependency.
- `crates/cosmolkit-msready/tests/{acceptance,invariants,source_data,release_regressions}.rs` — public detached-boundary/invariant/fixed-source regressions.
- `crates/cosmolkit/src/msready.rs`, `crates/cosmolkit/tests/msready.rs`, `crates/cosmolkit/examples/msready.rs` — thin facade, focused integration tests and explicit-feature example.
- `tools/testdata/{prepare_msready.py,test_prepare_msready.py}` — deterministic offline fixture preparation/checks.
- `testdata/msready/{README.md,manifest.json,rule_classification.json,fixtures/}` — immutable selected originals/selection hashes/licensing context.
- Managed validation/review/PR-note artifacts outside worktree — delivery evidence, not committed scratch.

## Dependencies and exclusive ownership seams

| Seam | Exclusive write scope | Gate and durable handoff |
|---|---|---|
| Source/fixture preparation | `tools/testdata/prepare_msready*`, `tools/testdata/test_prepare_msready.py`, `testdata/msready/**`, dedicated ledger/support contract | Tasks 1–3: immutable source hashes/classification/reviewed nitrate variants and preparation logs; expectations never refreshed to fit failures. |
| Detached chemistry | `crates/cosmolkit-msready/**`, initial workspace/lock delta | Tasks 4–19: early real-identity gate, complete stage tests, frozen rules/metal predicates/nitrate catalogue/provenance; handoff tested detached API/data hashes/logs. Rule/H/identity/finalization changes are tightly coupled and have one writer. |
| Existing mass owner | `crates/cosmolkit-descriptors/src/lib.rs` only | Tasks 16–17: narrow helper using existing kernel/private constant plus unchanged ionic regression output; reusable helper signature is the only integration contract. |
| Facade/integration | Named `crates/cosmolkit` files, feature CI, final lock delta | Tasks 20–23, after detached gate: opt-in methods/types/registration, unchanged bundles and full strict/offline validation. No chemistry reimplementation. |
| Read-only review/local delivery | Managed artifacts, then explicit intentional staging/commit | Tasks 24–27 after mandatory gates; AGENTS excluded at index and committed-diff checks. |

Execute the sequential ledger even for independently testable seams; parent assigns ownership explicitly. Root/lock ownership transfers sequentially. No overlapping simultaneous writers. Component completion must yield a durable reviewed contract/log handoff before integration-only work starts.

Logical dependencies: approved contract/fixtures precede package/exclusions; vocabulary/exact predicates precede identity; **Task 7 blocks all later chemistry**; connectivity/maps precede rule/H/accounting; charge/finalization/mass-owner gates precede complete pipeline; complete detached tests precede facade; documentation/final checks and read-only review precede staging/commit. Nitrate comparison metadata and input-domain validation share the frozen template boundary; changing either requires revalidation, not inferred broad repair.

## Source anchors and review strategy

- Native design is not a fictitious C++ port. Inside each adopted rule function include verbatim relevant raw MRV atom/bond excerpts or nitro reaction string, pinned path/hash and mapped edits/guards/invariants, with separate behavior and complexity review. Raw source files remain reviewable.
- Explicitly mark R/H expansion, isotope-aware O selection, rematching/fixed point, exact charged alkali guard and nitrate comparison templates as project adaptations; none proves original Indigo execution parity. Hashes/comments alone are not behavior evidence.
- For actual changed/ported lower-owner source functions follow verbatim pinned C/C++ in-function anchors and independent behavior/complexity markers. Reuse correct owner APIs rather than copying algorithms into the facade/standardizer.
- Review real local graph/H/charge/isotope/bond attributes and composed maps, not only strings/formula/key equality. Review bounded candidate enumeration/canonical ranking/cloning costs separately from behavioral reproduction; no unjustified complexity marker upgrades.

## Risks and stop conditions

1. Real InChI integration remains unproved, including all accepted Hg/Sn tuples and labelled alkali/metal states. Task 7 is blocking; unavailable independent support needs authorization, not identity omission/domain shrinking or mislabelled supported failures.
2. Five-rule Native scope is finite; no full historical 111-appearance graph/KNIME parity or EPA equivalence claim. Parent/final oracle approval freezes this corrected contract; later chemistry changes require an explicit revision.
3. Raw nitrate H2NO3 violates conventional strict N valence. It is a narrowly recognized exact exclusion-only input, never a validated retained analyte or ordinary proton variant. Nitric-acid/nitrate adaptations must remain catalogue/comparison evidence and never mutate the analysis graph.
4. Bonded `(+1,-1)` alkali–O input can fail pre-edit strict valence. Only the exact frozen repair defers that check; post-edit strict validation and formula/H/charge equality are mandatory. Rejected pairs cannot use charge/H guessing.
5. Owner Hg valence -1 does not define the supported metal domain; every accepted tuple/neighbor/H/radical predicate is enforced. Unknown valid complexes are structured unmodeled capabilities; invalid valence and failures on accepted inputs remain distinct.
6. Source MRV is static, not an execution oracle. H slots, isotope transfer, overlap/order, branch retention and nitro O roles need mechanistic tests; no fallback to non-executable previews/QSAR rules.
7. Renumbered fragment/H-fold IDs, isotope-labelled H or deduplication can lose accounting/multiplicity; complete input-atom disposition and all origins remain invariants.
8. Stereo suppression flags can erase isotopes; physically clear all derived stereo/perception carriers and avoid original-coordinate re-perception. Keep graph/radical/metal topology independent from identifier normalization.
9. Existing descriptor exact mass includes electron correction. Neutral-atom formula mass requires owner-backed correction, no copied constant/table or average-mass fallback; existing mass behavior/regressions remain unchanged.
10. Offline dependencies/corpora/pre-existing workspace failures may block mandatory pre-commit gates. Capture exact logs/head/state and stop for decision; never fetch, weaken CI/checkers, silently exclude members, ignore tests or commit a mandatory failure.
11. Archive licensing is attested by supplied provenance, not independently verified online. Retain MIT/CC BY/CC0 attribution/hash context; an explicit redistribution prohibition blocks copying that asset.
12. `cap-msready` must remain opt-in and all existing bundles/normative docs unchanged. Declared Python/JS names are not delivered bindings. AGENTS exceptions and unrelated original changes remain local; actual PR opening/push/remote actions are prohibited.

## Planning acceptance

The only write made by this planner is this managed revised-plan artifact. No project code/tests/plans/AGENTS edits, staging, Git operations, builds/tests, dependency setup or network actions occurred. Targeted read evidence resolved the oracle's feature/nitrate/metal contract gaps; runtime feasibility and mandatory validation remain future gates.

```acceptance-report
{
  "criteriaSatisfied": [
    {
      "id": "criterion-1",
      "status": "satisfied",
      "evidence": "Complete 27-task/55-step revision preserves original gates and freezes opt-in feature membership, exactly three nitrate comparison templates/adaptation metadata, and explicit Hg/Sn/neutral+charged alkali predicates with positive/negative tests."
    }
  ],
  "changedFiles": [
    "/home/analytit_admin/.pi/agent/sessions/--data-dev-COSMolKit--/subagent-artifacts/outputs/49b2d10f-600d-4cd5-8640-2299bd1c0d41/revised-plan.md"
  ],
  "testsAddedOrUpdated": [],
  "commandsRun": [],
  "validationOutput": [
    "Targeted reads confirmed exact existing full bundle policy, raw row-83 ON(=O)O/H2NO3 data, Sn 2/4 and Hg unconstrained owner valence, and separable selected sanitization/cleanup. No runtime validation claimed."
  ],
  "residualRisks": [
    "Real-engine InChI support for every accepted isotope/radical/Hg/Sn/alkali state remains a blocking proof gate.",
    "Static rule data cannot certify historical equivalence; archive licensing and mandatory offline/workspace gates remain unverified."
  ],
  "noStagedFiles": true,
  "diffSummary": "Managed revised plan only; no repository mutation or staging by this planner.",
  "reviewFindings": [
    "Oracle corrections are addressed in the proposed contract; final oracle and parent approval remain required before implementation."
  ],
  "manualNotes": "Use the context-builder report as primary handoff. Keep both local AGENTS exceptions outside the index/commit. The revised plan deliberately does not add cap-msready to full/default/core, does not chemically rewrite nitrate to match an exclusion, and does not relabel failures on accepted metal graph states as unsupported. No tests or commands were run during planning."
}
```
