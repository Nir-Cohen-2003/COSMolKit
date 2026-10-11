MS Ready SMILES reference data
Research date: 2026-10-08

This archive accompanies MSREADY_SMILES_implementation_reference.pdf/docx.
It is source and inspection material, not an executable MS-Ready implementation.

Baseline: kmansouri/MS-ready GitHub release v2.5.8, published 2022-02-02.
Asset inspected: MS-ready_2.5.8-knime_4.5.knwf.
Release: https://github.com/kmansouri/MS-ready/releases/tag/v2.5.8
Paper: https://doi.org/10.1186/s13321-018-0299-2
Repository: https://github.com/kmansouri/MS-ready
Repository declares an MIT license. See source URLs and retain attribution.
Supplementary data: authors' original files; source article licenses and data
terms are supplied on the source pages, including CC BY 4.0 and the applicable
CC0 data waiver. No claim of authorship of original reaction or exclusion data.

FILES
connected_rule_inventory.json: static connected definition inventory. Contains
four pre-component disconnection definitions and 111 normalization definition
appearances (including a repeated five-entry nitro table). Counts are not counts
of unique tautomer rules. Stage ordering follows the visible source graph.
For conditional branching, valence bypass, column routing and rule selection,
consult the preserved workflow_configuration XML. No runtime conformance claim.

mrv/: original connected ChemAxon MRV reaction drawings, decoded from the XML
string encoding. Formal charges, explicit hydrogens, R-group labels, radical
states, special valences and query bonds have not been flattened.

Any display_reaction_smiles_not_executable field is a convenience preview only.
It was generated with Indigo 1.46.0 from a molecule-mode CML import and can lose
query semantics. Some previews explicitly fail on query bonds. These previews
are not SMIRKS, not a query-preserving conversion, and not safe to execute.
Raw reaction_string entries are retained from the source Table Creator. Their
Indigo query import and R-group handling still require a compatible engine.

workflow_configuration/: original settings.xml and workflow.knime files only;
cached tables, binary state and complete bundled binaries are omitted.
Disconnected example nodes remain in this original configuration tree.

exclusion_list_v2_5_8.json/.csv: the 32 rows of Table Creator 2202, as stored,
including original names, structure strings and numeric columns. These are not
recanonicalized or corrected. Original mass fields are source data, not a
recommended modern mass calculator. Exact structure semantics require attention
to hydrogen and charge states and the matcher configuration.

exclusion_list_2018_original.sdf: unchanged Additional file 1 of the 2018 paper.
It is distinct from the extracted later table. Valence/implicit-H conventions
must be inspected before assuming identical serialized molecular graphs.

sample_input/ and sample_output/: unchanged v2.5.8 release sample assets.
The summary output has 654 data rows; 63 Canonical_MSr values contain slash or
backslash stereo symbols, and none contains @. This is provided sample behavior,
not a fresh execution. InChI strings in the summary CSV should not be assumed a
lossless authoritative identifier serialization; prefer native SDF/SMI and
regenerate/verify identifiers.

qsar_rules_2024_original_NOT_msready.txt: unchanged Additional file 2 to the
2024 QSAR-ready paper https://doi.org/10.1186/s13321-024-00814-3 . Related rules,
not the 2022 MS-ready release profile, and not current EPA production proof.

manifest.json: hashes, original URLs, extracted-file inventory and inspection
scope. The downloaded complete KNIME asset is not included because the useful
configuration/rule material is extracted here and its original URL is retained.

See the report for published intent, evidence limits, conservative engineering
choices, examples, implementation interfaces and conformance requirements.
