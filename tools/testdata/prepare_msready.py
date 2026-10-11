#!/usr/bin/env python3
"""Prepare the pinned, small MSReady source fixtures from a local archive.

This tool is development-only. It never downloads inputs and refuses source
files that do not match the archive's recorded SHA-256 digests.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1.0.0"
ARCHIVE_URL = "https://github.com/kmansouri/MS-ready/releases/download/v2.5.8/MS-ready_2.5.8-knime_4.5.knwf"
ARCHIVE_SHA256 = "042516c36088ffa8b9e9ae8e6e2e9154188cdb790c8ac3c510643bb7b01acc5c"
SAMPLE_INPUT = "sample_input/Sample_input.smi"
SAMPLE_DISCARDED = "sample_output/Sample_input_DiscardedStructures.csv"
SAMPLE_SELECTION_LINES = (1, 2, 3, 4, 5, 6, 7, 8)
NITRATE_TEMPLATES = [
    {
        "template_id": "83_raw_h2no3",
        "smiles": "ON(=O)O",
        "formula": "H2NO3",
        "net_charge": 0,
        "status": "source_verbatim_exclusion_only",
        "adaptation_id": None,
        "comparison_only": True,
    },
    {
        "template_id": "83_nitric_acid",
        "smiles": "[OH][N+](=O)[O-]",
        "formula": "HNO3",
        "net_charge": 0,
        "status": "adopted_comparison_template",
        "adaptation_id": "nitrate_comparison_charge_adaptation_v1",
        "comparison_only": True,
    },
    {
        "template_id": "83_nitrate",
        "smiles": "[O-][N+](=O)[O-]",
        "formula": "NO3",
        "net_charge": -1,
        "status": "adopted_comparison_template",
        "adaptation_id": "nitrate_comparison_charge_adaptation_v1",
        "comparison_only": True,
    },
]


class PreparationError(RuntimeError):
    """A required local source is absent, corrupt, or inconsistent."""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_source_manifest(root: Path) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    manifest_path = root / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PreparationError(f"cannot read archive manifest {manifest_path}: {exc}") from exc
    indexed = {entry["path"]: entry for entry in manifest.get("files", [])}
    if not indexed:
        raise PreparationError("archive manifest has no file inventory")
    return manifest, indexed


def verified_bytes(root: Path, indexed: dict[str, dict[str, Any]], relpath: str) -> bytes:
    expected = indexed.get(relpath)
    if expected is None:
        raise PreparationError(f"source is not pinned in archive manifest: {relpath}")
    source = root / relpath
    try:
        data = source.read_bytes()
    except OSError as exc:
        raise PreparationError(f"required source is missing: {source}: {exc}") from exc
    actual = sha256(data)
    if actual != expected["sha256"]:
        raise PreparationError(
            f"source checksum mismatch for {relpath}: expected {expected['sha256']}, got {actual}"
        )
    if len(data) != expected["bytes"]:
        raise PreparationError(f"source byte count mismatch for {relpath}")
    return data


def destination_files(root: Path, indexed: dict[str, dict[str, Any]]) -> dict[str, tuple[str, bytes]]:
    selected = [
        "README.txt",
        "connected_rule_inventory.json",
        "exclusion_list_v2_5_8.csv",
        "exclusion_list_v2_5_8.json",
        "mrv/Metanode (#1222)__MarvinSketch (#2131).mrv",
        "mrv/Metanode (#1222)__MarvinSketch (#2135).mrv",
        "mrv/Metanode (#1222)__MarvinSketch (#2136).mrv",
        "mrv/Tauto_Neut (#3652)__Metanode (#3383)__MarvinSketch (#3381).mrv",
        "mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2136).mrv",
        "mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2253).mrv",
        "mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2256).mrv",
        "sample_input/Sample_input.smi",
        "sample_output/Sample_input_DiscardedStructures.csv",
    ]
    # These settings contain the raw nitro entries at initial/final stages and
    # preserve the source row context rather than relying on preview strings.
    selected.extend(
        [
            "workflow_configuration/Tauto_Neut (#3652)/N mesos (#3406)/Table Creator (#2133)/settings.xml",
            "workflow_configuration/Tauto_Neut (#3652)/N mesos (#3430)/Table Creator (#2158)/settings.xml",
        ]
    )
    result: dict[str, tuple[str, bytes]] = {}
    for relpath in selected:
        if relpath in (SAMPLE_INPUT, SAMPLE_DISCARDED):
            continue
        result["fixtures/source/" + relpath] = (relpath, verified_bytes(root, indexed, relpath))

    sample_bytes = verified_bytes(root, indexed, SAMPLE_INPUT)
    lines = sample_bytes.splitlines(keepends=True)
    chosen = []
    for line_number in SAMPLE_SELECTION_LINES:
        if line_number > len(lines):
            raise PreparationError(f"sample input has no line {line_number}")
        raw = lines[line_number - 1]
        text = raw.decode("utf-8").rstrip("\r\n")
        fields = text.split("\t", 1)
        if len(fields) != 2:
            raise PreparationError(f"sample line {line_number} is not SMILES plus record ID")
        chosen.append(
            {
                "line_number": line_number,
                "record_id": fields[1],
                "smiles": fields[0],
                "line_sha256": sha256(raw),
            }
        )
    selected_lines = b"".join(lines[number - 1] for number in SAMPLE_SELECTION_LINES)
    result["fixtures/native_sample_selection.smi"] = (SAMPLE_INPUT, selected_lines)
    result["fixtures/native_discarded_sample.csv"] = (
        SAMPLE_DISCARDED,
        verified_bytes(root, indexed, SAMPLE_DISCARDED),
    )

    # Preserve the raw source table; row 83 is verified below without rewriting.
    exclusions = json.loads(verified_bytes(root, indexed, "exclusion_list_v2_5_8.json"))
    row83 = [row for row in exclusions if row["ID"] == "83"]
    if len(row83) != 1 or row83[0].get("CANONICAL_SMILES") != "ON(=O)O":
        raise PreparationError("pinned exclusion row 83 is absent or no longer the raw H2NO3 depiction")
    nitrate = {
        "schema_version": SCHEMA_VERSION,
        "source_id": "83",
        "source_field": "CANONICAL_SMILES",
        "source_smiles": row83[0]["CANONICAL_SMILES"],
        "source_formula": row83[0]["Formula"],
        "source_net_charge": 0,
        "adaptation_id": "nitrate_comparison_charge_adaptation_v1",
        "templates": NITRATE_TEMPLATES,
        "source_to_template_metadata": {
            "raw_to_acid": {"delta_n_charge": 1, "delta_single_o_charge": -1, "delta_h": -1, "delta_net_charge": 0},
            "raw_to_nitrate": {"delta_n_charge": 1, "delta_single_o_charge": -1, "delta_h": -2, "delta_net_charge": -1},
            "is_component_edit": False,
        },
    }
    result["fixtures/nitrate_templates.json"] = (
        "exclusion_list_v2_5_8.json#ID=83",
        (json.dumps(nitrate, indent=2, sort_keys=True) + "\n").encode(),
    )
    sample_selection = {
        "source_path": SAMPLE_INPUT,
        "source_sha256": indexed[SAMPLE_INPUT]["sha256"],
        "selection_method": "fixed one-based source line numbers; preserve complete original records",
        "selected_records": chosen,
    }
    result["fixtures/native_sample_selection.json"] = (
        SAMPLE_INPUT,
        (json.dumps(sample_selection, indent=2, sort_keys=True) + "\n").encode(),
    )
    return result


def build_manifest(
    files: dict[str, tuple[str, bytes]], indexed: dict[str, dict[str, Any]],
    archive_manifest: dict[str, Any],
) -> dict[str, Any]:
    file_entries = []
    for dest, (source, data) in sorted(files.items()):
        if source.endswith("#ID=83"):
            source_info = {"path": source, "sha256": indexed["exclusion_list_v2_5_8.json"]["sha256"]}
        else:
            original_sources = archive_manifest.get("original_sources", [])
            if source.startswith("sample_input/"):
                url = next((entry["url"] for entry in original_sources if entry.get("local_directory") == "sample_input"), ARCHIVE_URL)
            elif source.startswith("sample_output/"):
                url = next((entry["url"] for entry in original_sources if entry.get("local_directory") == "sample_output"), ARCHIVE_URL)
            else:
                url = archive_manifest.get("workflow_asset_url", ARCHIVE_URL)
            source_info = {
                "path": source,
                "sha256": indexed[source]["sha256"],
                "url": url,
                "release": archive_manifest["baseline"],
            }
        file_entries.append({"path": dest, "sha256": sha256(data), "bytes": len(data), "source": source_info})
    return {
        "schema_version": SCHEMA_VERSION,
        "profile": "cosmolkit_msready_v1",
        "profile_commitment": "Native published-concept implementation; not historical equivalence",
        "source_archive": {
            "release": archive_manifest["baseline"],
            "asset_url": archive_manifest.get("workflow_asset_url", ARCHIVE_URL),
            "asset_sha256": ARCHIVE_SHA256,
            "static_inspection_only": True,
        },
        "license_context": {
            "archive_repository": "MS-ready repository declares MIT",
            "original_data": "archive README attributes source article data under CC BY 4.0 and applicable CC0 waiver; retain source attribution and terms",
            "attribution": "kmansouri/MS-ready v2.5.8; McEachran et al. 2018, DOI 10.1186/s13321-018-0299-2",
        },
        "selection": {
            "sample_lines_1_based": list(SAMPLE_SELECTION_LINES),
            "sample_record_count": len(SAMPLE_SELECTION_LINES),
            "copied_full_sample": False,
            "copied_xls_or_sdf": False,
            "copied_qsar_supplement": False,
            "expected_outputs_regenerated": False,
        },
        "files": file_entries,
    }


def prepare(source_root: Path, output: Path) -> None:
    source_root = source_root.resolve()
    archive_manifest, indexed = load_source_manifest(source_root)
    files = destination_files(source_root, indexed)
    manifest_bytes = (json.dumps(build_manifest(files, indexed, archive_manifest), indent=2, sort_keys=True) + "\n").encode()
    files["manifest.json"] = ("generated_manifest", manifest_bytes)

    output.mkdir(parents=True, exist_ok=True)
    for relative, (_, data) in files.items():
        target = output / relative
        if target.exists():
            if target.read_bytes() != data:
                raise PreparationError(f"refusing to replace immutable fixture with different bytes: {target}")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    manifest_target = output / "manifest.json"
    if manifest_target.exists() and manifest_target.read_bytes() != manifest_bytes:
        raise PreparationError(f"refusing to replace immutable manifest: {manifest_target}")
    if not manifest_target.exists():
        manifest_target.write_bytes(manifest_bytes)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True, help="local extracted archive root; no network access")
    parser.add_argument("--output", type=Path, required=True, help="fixture output directory")
    args = parser.parse_args()
    try:
        prepare(args.source_root, args.output)
    except PreparationError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
