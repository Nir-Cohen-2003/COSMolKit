use cosmolkit_msready::*;
use std::io::Write;
use std::process::{Command, Stdio};

const RAW: &str = include_str!("../src/data_exclusions_source.json");
const VARIANTS: &str = include_str!("../src/data_exclusion_variants.json");
const NITRATE: &str = include_str!("../src/data_nitrate_templates.json");
const RULES: &str = include_str!("../src/data/rules.json");

// Test-only checksum dependency: repository fixture preparation already requires
// Python3. Its standard-library hashlib reads immutable bytes, no environment
// setup, fixture refresh, chemistry oracle, external path or network access.
fn sha256(bytes: &[u8]) -> String {
    let mut child = Command::new("python3")
        .args([
            "-c",
            "import hashlib,sys;print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("repository Python3 fixture prerequisite");
    child.stdin.take().unwrap().write_all(bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    String::from_utf8(output.stdout).unwrap().trim().to_owned()
}
#[test]
fn packaged_source_hashes_and_profile_writer_options_are_frozen_and_runtime_owned() {
    let info = run_info();
    assert_eq!(info.profile_id, "cosmolkit_msready_v1");
    assert_eq!(info.profile_version, "1.0.0");
    assert_eq!(info.behavior, "Native");
    for (bytes, expected) in [
        (RAW.as_bytes(), info.exclusion_sha256),
        (VARIANTS.as_bytes(), info.exclusion_variants_sha256),
        (NITRATE.as_bytes(), info.nitrate_templates_sha256),
        (RULES.as_bytes(), info.rule_sha256),
    ] {
        assert_eq!(sha256(bytes), expected);
    }
    assert_eq!(
        RAW.as_bytes(),
        include_bytes!("../../../testdata/msready/fixtures/source/exclusion_list_v2_5_8.json")
    );
    assert_eq!(
        NITRATE.as_bytes(),
        include_bytes!("../../../testdata/msready/fixtures/nitrate_templates.json")
    );
    let p = info.writer_params;
    assert!(p.do_isomeric_smiles && p.canonical && p.ignore_atom_map_numbers);
    assert!(
        !p.do_kekule
            && !p.clean_stereo
            && !p.include_dative_bonds
            && !p.all_bonds_explicit
            && !p.all_hydrogens_explicit
    );
    assert!(p.rooted_at_atom.is_none());
    assert_eq!(info.inchi_engine_version, "1.07.5");
    for text in [
        include_str!("../src/pipeline.rs"),
        include_str!("../src/lib.rs"),
    ] {
        assert!(!text.contains("/data/dev/"));
        assert!(!text.contains("testdata/"));
        assert!(!text.contains("msready-reference"));
    }
}
#[test]
fn all_32_raw_source_rows_are_complete_pipeline_exclusions_including_alkali_splits() {
    let rows: serde_json::Value = serde_json::from_str(RAW).unwrap();
    let rows = rows.as_array().unwrap();
    assert_eq!(rows.len(), 32);
    for row in rows {
        let id = row["ID"].as_str().unwrap();
        let (field, text) = if let Some(s) = row["Structure(SMILES)"].as_str() {
            ("Structure(SMILES)", s)
        } else {
            (
                "CANONICAL_SMILES",
                row["CANONICAL_SMILES"].as_str().unwrap(),
            )
        };
        let result = standardize_smiles(text, &Default::default())
            .unwrap_or_else(|e| panic!("source row {id} {text}: {e:?}"));
        assert!(
            matches!(result.status, MsReadyStatus::NoEligibleComponent),
            "source row {id}"
        );
        assert_eq!(
            result.atom_dispositions.len(),
            result.original.topology.atoms.len()
        );
        if matches!(id, "180" | "181" | "182") {
            assert_eq!(
                result
                    .events
                    .iter()
                    .filter(|e| e.kind == MsReadyEventKind::AlkaliDisconnection)
                    .count(),
                1
            );
            assert!(
                result
                    .exclusions
                    .iter()
                    .any(|e| e.origin.atom_indices == vec![0])
            );
            assert!(
                result
                    .exclusions
                    .iter()
                    .any(|e| e.origin.atom_indices == vec![1])
            );
        } else {
            let exact = result
                .exclusions
                .iter()
                .find(|e| e.source_id == id)
                .unwrap();
            assert_eq!(exact.source_field, field);
            assert_eq!(exact.source_sha256, EXCLUSION_SHA256);
        }
    }
    let carbonate = standardize_smiles("O=C(O)O", &Default::default()).unwrap();
    let ids = carbonate
        .exclusions
        .iter()
        .map(|e| e.source_id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(ids, std::collections::BTreeSet::from(["29", "155"]));
    let row83 = rows.iter().find(|row| row["ID"] == "83").unwrap();
    assert_eq!(row83["Formula"], "H2NO3");
    assert!(row83.get("Structure(SMILES)").is_none());
}
#[test]
fn reviewed_variant_catalogue_has_no_runtime_proton_ladder_and_all_variants_are_tested() {
    let variants: serde_json::Value = serde_json::from_str(VARIANTS).unwrap();
    let variants = variants["variants"].as_array().unwrap();
    assert_eq!(variants.len(), 26);
    for variant in variants {
        let text = variant["smiles"].as_str().unwrap();
        let result = standardize_smiles(text, &Default::default())
            .unwrap_or_else(|e| panic!("variant {} {text}: {e:?}", variant["variant_id"]));
        assert!(
            matches!(result.status, MsReadyStatus::NoEligibleComponent),
            "{text}"
        );
    }
    let classification: serde_json::Value = serde_json::from_str(include_str!(
        "../../../testdata/msready/rule_classification.json"
    ))
    .unwrap();
    assert_eq!(classification["appearance_counts"]["total"], 115);
    assert_eq!(
        classification["appearance_counts"]["adopted_normalization"],
        6
    );
    assert_eq!(
        classification["appearance_counts"]["unadopted_normalization"],
        105
    );
    assert_eq!(
        classification["nitrate_templates"]
            .as_array()
            .unwrap()
            .len(),
        3
    );
    assert_eq!(
        classification["nitrate_adaptation_is_component_edit"],
        false
    );
    for template in classification["nitrate_templates"].as_array().unwrap() {
        let text = template["smiles"].as_str().unwrap();
        let result = standardize_smiles(text, &Default::default()).unwrap();
        assert!(result.events.is_empty());
        assert!(
            result
                .exclusions
                .iter()
                .any(|e| e.template_id.as_deref() == template["template_id"].as_str())
        );
    }
}
#[test]
fn immutable_selected_fixture_manifest_hashes_are_checked_without_preparation_or_mutation() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/msready");
    let manifest: serde_json::Value =
        serde_json::from_str(include_str!("../../../testdata/msready/manifest.json")).unwrap();
    let files = manifest["files"].as_array().unwrap();
    assert!(files.len() >= 10);
    for file in files {
        let path = root.join(file["path"].as_str().unwrap());
        let bytes = std::fs::read(path).unwrap();
        assert_eq!(bytes.len() as u64, file["bytes"].as_u64().unwrap());
        assert_eq!(sha256(&bytes), file["sha256"].as_str().unwrap());
    }
}
