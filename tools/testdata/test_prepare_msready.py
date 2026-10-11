from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import prepare_msready as prep

REPOSITORY = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = REPOSITORY / "testdata" / "msready"
ARCHIVE_FIXTURES = FIXTURE_ROOT / "fixtures" / "source"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class PrepareMsreadyTests(unittest.TestCase):
    def source_tree(self, root: Path) -> None:
        """Create a minimal, fully hash-pinned local archive for tool tests."""
        inputs = {
            "README.txt": b"MIT archive; CC BY 4.0 and CC0 original-data terms.\n",
            "connected_rule_inventory.json": b'{"stages":[]}\n',
            "exclusion_list_v2_5_8.csv": b'ID,Structure(SMILES)\n83,ON(=O)O\n',
            "exclusion_list_v2_5_8.json": json.dumps(
                [{"ID": "83", "CANONICAL_SMILES": "ON(=O)O", "Formula": "H2NO3"}]
            ).encode(),
            "mrv/Metanode (#1222)__MarvinSketch (#2131).mrv": b"K rule fixture\n",
            "mrv/Metanode (#1222)__MarvinSketch (#2135).mrv": b"Li rule fixture\n",
            "mrv/Metanode (#1222)__MarvinSketch (#2136).mrv": b"Na rule fixture\n",
            "mrv/Tauto_Neut (#3652)__Metanode (#3383)__MarvinSketch (#3381).mrv": b"enamine fixture\n",
            "mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2136).mrv": b"ynol fixture\n",
            "mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2253).mrv": b"enol fixture\n",
            "mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2256).mrv": b"azide fixture\n",
            "workflow_configuration/Tauto_Neut (#3652)/N mesos (#3406)/Table Creator (#2133)/settings.xml": b"initial nitro settings\n",
            "workflow_configuration/Tauto_Neut (#3652)/N mesos (#3430)/Table Creator (#2158)/settings.xml": b"final nitro settings\n",
            "sample_input/Sample_input.smi": b"CCO\tcase-1\nCCN\tcase-2\nCCC\tcase-3\nC\tcase-4\nN\tcase-5\nO\tcase-6\nF\tcase-7\nCl\tcase-8\n",
            "sample_output/Sample_input_DiscardedStructures.csv": b"Molecule name,RowID,Original_SMILES,Error\n",
        }
        inventory = []
        for relpath, data in inputs.items():
            target = root / relpath
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            inventory.append({"path": relpath, "bytes": len(data), "sha256": digest(data)})
        (root / "manifest.json").write_text(
            json.dumps({"baseline": "test pinned archive", "files": inventory}, sort_keys=True),
            encoding="utf-8",
        )

    def tree_hashes(self, root: Path) -> dict[str, str]:
        return {
            path.relative_to(root).as_posix(): digest(path.read_bytes())
            for path in sorted(root.rglob("*"))
            if path.is_file()
        }

    def test_preparation_is_deterministic_and_pins_selected_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temp = Path(temporary)
            source = temp / "source"
            self.source_tree(source)
            first, second = temp / "first", temp / "second"
            prep.prepare(source, first)
            prep.prepare(source, second)
            self.assertEqual(self.tree_hashes(first), self.tree_hashes(second))
            manifest = json.loads((first / "manifest.json").read_text())
            self.assertEqual(manifest["selection"]["sample_record_count"], 8)
            selection = json.loads((first / "fixtures/native_sample_selection.json").read_text())
            self.assertEqual([r["line_number"] for r in selection["selected_records"]], list(range(1, 9)))
            self.assertTrue(all(len(r["line_sha256"]) == 64 for r in selection["selected_records"]))
            self.assertFalse(manifest["selection"]["copied_xls_or_sdf"])
            self.assertFalse(manifest["selection"]["copied_qsar_supplement"])

    def test_missing_source_fails_without_skip(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "source"
            self.source_tree(root)
            (root / "mrv/Metanode (#1222)__MarvinSketch (#2131).mrv").unlink()
            with self.assertRaisesRegex(prep.PreparationError, "required source is missing"):
                prep.prepare(root, Path(temporary) / "out")

    def test_corrupt_source_fails_checksum(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "source"
            self.source_tree(root)
            target = root / "exclusion_list_v2_5_8.json"
            target.write_bytes(target.read_bytes() + b"corrupt")
            with self.assertRaisesRegex(prep.PreparationError, "checksum mismatch"):
                prep.prepare(root, Path(temporary) / "out")

    def test_existing_fixture_bytes_cannot_be_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "source"
            output = Path(temporary) / "out"
            self.source_tree(root)
            prep.prepare(root, output)
            selected = output / "fixtures/source/README.txt"
            selected.write_bytes(b"locally modified\n")
            before = selected.read_bytes()
            with self.assertRaisesRegex(prep.PreparationError, "refusing to replace immutable fixture"):
                prep.prepare(root, output)
            self.assertEqual(selected.read_bytes(), before)

    def test_committed_catalogue_classifies_all_appearances_and_exclusions(self) -> None:
        catalogue_path = FIXTURE_ROOT / "rule_classification.json"
        catalogue = json.loads(catalogue_path.read_text(encoding="utf-8"))
        appearances = catalogue["appearances"]
        self.assertEqual(len(appearances), 115)
        self.assertEqual(sum(a["stage_id"] == "pre_component_disconnection" for a in appearances), 4)
        normalizations = [a for a in appearances if a["stage_id"] != "pre_component_disconnection"]
        self.assertEqual(len(normalizations), 111)
        self.assertEqual(sum(a["adoption_status"] == "adopted" for a in normalizations), 6)
        self.assertEqual(sum(a["adoption_status"] == "adopted" for a in appearances if a["stage_id"] == "pre_component_disconnection"), 3)
        self.assertEqual(len(catalogue["exclusions"]), 32)
        self.assertEqual(len({row["source_id"] for row in catalogue["exclusions"]}), 32)
        by_id = {row["source_id"]: row for row in catalogue["exclusions"]}
        self.assertEqual(by_id["58"]["selected_field"], "CANONICAL_SMILES")
        self.assertEqual(by_id["83"]["selected_field"], "CANONICAL_SMILES")
        self.assertEqual(by_id["29"]["selected_value"], by_id["155"]["selected_value"])
        self.assertNotEqual(by_id["29"]["source_id"], by_id["155"]["source_id"])
        self.assertNotIn("qsar_rules_2024_original_NOT_msready.txt", json.dumps(catalogue))

    def test_nitro_source_attribution_is_stage_specific(self) -> None:
        catalogue = json.loads((FIXTURE_ROOT / "rule_classification.json").read_text())
        initial = next(a for a in catalogue["appearances"] if a["appearance_id"] == "initial_nitro:03:03")
        final = next(a for a in catalogue["appearances"] if a["appearance_id"] == "final_nitro:03:03")
        self.assertEqual(initial["raw_reaction"], final["raw_reaction"])
        self.assertEqual(initial["source_path"], "workflow_configuration/Tauto_Neut (#3652)/N mesos (#3406)/Table Creator (#2133)/settings.xml")
        self.assertEqual(final["source_path"], "workflow_configuration/Tauto_Neut (#3652)/N mesos (#3430)/Table Creator (#2158)/settings.xml")
        manifest = json.loads((FIXTURE_ROOT / "manifest.json").read_text())
        copied_sources = {entry["source"]["path"]: entry["sha256"] for entry in manifest["files"]}
        for appearance in (initial, final):
            self.assertEqual(appearance["source_sha256"], copied_sources[appearance["source_path"]])
        self.assertNotEqual(initial["source_sha256"], final["source_sha256"])

    def test_exact_three_nitrate_templates_and_raw_source_remains_unchanged(self) -> None:
        nitrate = json.loads((FIXTURE_ROOT / "fixtures/nitrate_templates.json").read_text())
        self.assertEqual(len(nitrate["templates"]), 3)
        self.assertEqual([t["template_id"] for t in nitrate["templates"]], [
            "83_raw_h2no3", "83_nitric_acid", "83_nitrate"
        ])
        self.assertEqual(nitrate["source_smiles"], "ON(=O)O")
        self.assertEqual(nitrate["source_formula"], "H2NO3")
        self.assertFalse(nitrate["source_to_template_metadata"]["is_component_edit"])
        self.assertTrue(all(t["comparison_only"] for t in nitrate["templates"]))

    def test_fixture_manifest_hashes_and_fixture_tree_are_immutable(self) -> None:
        before = self.tree_hashes(FIXTURE_ROOT)
        manifest = json.loads((FIXTURE_ROOT / "manifest.json").read_text(encoding="utf-8"))
        for entry in manifest["files"]:
            data = (FIXTURE_ROOT / entry["path"]).read_bytes()
            self.assertEqual(len(data), entry["bytes"], entry["path"])
            self.assertEqual(digest(data), entry["sha256"], entry["path"])
            self.assertTrue(entry["source"].get("path"), entry["path"])
        after = self.tree_hashes(FIXTURE_ROOT)
        self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()
