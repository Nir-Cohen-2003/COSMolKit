# MSReady Native profile fixtures

These files are a small, immutable source/reference inventory for the approved
`cosmolkit_msready_v1` Native profile. They do **not** represent an executed
v2.5.8 workflow, an executable rule port, or a historical parity oracle. The
workflow bundle was inspected statically only.

## Provenance and license context

- Baseline: `kmansouri/MS-ready` v2.5.8 (2022-02-02),
  <https://github.com/kmansouri/MS-ready/releases/tag/v2.5.8>.
- Downloaded KNIME asset URL:
  <https://github.com/kmansouri/MS-ready/releases/download/v2.5.8/MS-ready_2.5.8-knime_4.5.knwf>.
- Downloaded asset SHA-256:
  `042516c36088ffa8b9e9ae8e6e2e9154188cdb790c8ac3c510643bb7b01acc5c`.
- The archive README reports that the repository declares MIT and that
  included original article data is subject to source-page terms including CC
  BY 4.0 and the applicable CC0 waiver. Original rule/exclusion data are not
  claimed as COSMolKit-authored; retain attribution to the source repository
  and McEachran et al. (2018), DOI `10.1186/s13321-018-0299-2`.
- `fixtures/source/README.txt` is the archive's original terms and scope note.
  `manifest.json` binds every copied/selected byte sequence to its source path,
  source checksum, URL, selection identity, and resulting checksum.

## Contents and selection

- `fixtures/source/exclusion_list_v2_5_8.{csv,json}`: all 32 original rows,
  unchanged. The selected match field is `Structure(SMILES)` when present and
  nonempty, otherwise the recorded `CANONICAL_SMILES`; duplicate carbonate
  IDs remain distinct. Source row 83 remains the original `ON(=O)O` / H2NO3.
- `fixtures/source/connected_rule_inventory.json`: the 115 source appearances
  and visible source ordering. `rule_classification.json` separately labels
  each appearance with source path/hash, adoption status, and rationale; it
  records 6 adopted of 111 normalization appearances and 3 adopted of 4
  pre-component appearances. Unadopted definitions are outside the frozen
  profile or conflict with its radical/query semantics. No 2024 QSAR supplement
  is imported.
- `fixtures/source/mrv/` and selected `workflow_configuration/.../settings.xml`:
  original, unflattened selected MRV and nitro table settings. The initial
  nitro table is `N mesos (#3406)/Table Creator (#2133)`; final nitro is
  `N mesos (#3430)/Table Creator (#2158)` per the inventory and archive
  manifest. Their row-3 raw reaction strings agree but source paths/hashes stay
  distinct. Display reaction SMILES in the inventory are not executable
  transformations.
- `fixtures/native_sample_selection.smi` and
  `fixtures/native_sample_selection.json`: fixed native sample input lines
  1–8, selected by one-based line number; source record IDs and per-line SHA-256
  values are retained. This is an input selection, not expected profile output.
- `fixtures/native_discarded_sample.csv`: the native discarded-structure CSV
  is retained as source context only. Neither sample output nor XLS/SDF assets
  are used as authoritative expected results.
- `fixtures/nitrate_templates.json`: exactly three row-83 comparison templates:
  verbatim raw H2NO3, comparison-only nitric acid, and comparison-only nitrate.
  The versioned source-to-template charge/H deltas are metadata only, not
  component edits or proton-transfer events.
- `manifest.json`: schema and source/selection/output identities.

The selected source inventory deliberately excludes the full sample input,
XLS/SDF files, sample output structures, and `qsar_rules_2024_original_NOT_msready.txt`.

## Local preparation and tests

Preparation needs a local extracted source root as an explicit argument; there
is no default path, download, or network fallback:

```sh
python3 tools/testdata/prepare_msready.py \
  --source-root /path/to/MSREADY_reference_data \
  --output testdata/msready
```

The tool validates each input against the archive manifest and refuses to
replace existing fixture bytes with different contents. Tests validate
preparation determinism, hash identity, missing/corrupt-source failures,
fixture immutability, source appearance counts, and the exact three nitrate
templates. Normal tests do not regenerate committed fixture expectations.
