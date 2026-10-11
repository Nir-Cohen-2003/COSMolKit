use cosmolkit_msready::*;

#[test]
fn selected_native_release_records_have_explicit_native_outcomes_not_historical_parity() {
    let selected: serde_json::Value = serde_json::from_str(include_str!(
        "../../../testdata/msready/fixtures/native_sample_selection.json"
    ))
    .unwrap();
    let records = selected["selected_records"].as_array().unwrap();
    assert_eq!(records.len(), 8);
    for (index, row) in records.iter().enumerate() {
        let text = row["smiles"].as_str().unwrap();
        let record_id = row["record_id"].as_str().unwrap();
        let params = MsReadyParams {
            record_id: Some(record_id.into()),
            ..Default::default()
        };
        let result = standardize_smiles(text, &params);
        if matches!(index, 3 | 4 | 5 | 7) {
            // Ge, Zr, Mg organometallic and Sm graphs are outside the frozen
            // retained-metal tuples (search-owner M includes Ge, excludes Si).
            // Native refuses them structurally, never omits them as success.
            let failure = result.unwrap_err();
            assert_eq!(
                failure.kind,
                MsReadyFailureKind::UnsupportedChemistry,
                "record {record_id}: {failure:?}"
            );
            assert_eq!(
                failure.evidence.input_link.record_id.as_deref(),
                Some(record_id)
            );
            assert_eq!(
                failure.evidence.input_link.source_smiles.as_deref(),
                Some(text)
            );
            assert!(failure.original.is_some());
            if index == 4 {
                assert!(!failure.evidence.exclusions.is_empty());
            }
        } else {
            let result = result.unwrap_or_else(|e| panic!("record {record_id}: {e:?}"));
            assert_eq!(result.input_link.record_id.as_deref(), Some(record_id));
            let MsReadyStatus::Products(products) = &result.status else {
                panic!("record {record_id} was silently omitted")
            };
            assert_eq!(products.len(), 1);
            let product = &products[0];
            assert_eq!(product.inchi_key.len(), 27);
            assert!(product.inchi.starts_with("InChI=1S/"));
            assert!(!product.canonical_smiles.contains(['@', '/', '\\']));
            if index == 0 {
                // Explicit Native difference: no carbon-presence filter. The
                // two eligible silicate anions remain one identity/count two;
                // archived Na+/Al3+ exclusions retain their own atom origins.
                assert_eq!(product.formula, "O3Si-2");
                assert_eq!(product.formal_charge, -2);
                assert_eq!(product.multiplicity, 2);
                assert_eq!(product.topology.atoms.len(), 4);
                assert_eq!(product.topology.bonds.len(), 3);
                assert!(
                    product
                        .isotope_composition
                        .iter()
                        .all(|e| e.atomic_number != 6)
                );
                assert!(!result.exclusions.is_empty());
            } else {
                // These pinned ordinary organic records have no adopted local
                // rewrite. Aromatic canonicalization is representational; no
                // CSV-derived identifier or historical string oracle is used.
                assert!(result.events.is_empty());
                assert_eq!(product.multiplicity, 1);
                assert_eq!(
                    product.topology.atoms.len(),
                    result.original.topology.atoms.len()
                );
                assert_eq!(
                    product.topology.bonds.len(),
                    result.original.topology.bonds.len()
                );
                assert!(result.exclusions.is_empty());
            }
            let again =
                standardize(&product.topology, &Default::default(), &Default::default()).unwrap();
            let MsReadyStatus::Products(again) = again.status else {
                panic!("non-idempotent release record")
            };
            assert_eq!(again[0].inchi_key, product.inchi_key);
            assert_eq!(again[0].canonical_smiles, product.canonical_smiles);
            assert_eq!(
                result.atom_dispositions.len(),
                result.original.topology.atoms.len()
            );
        }
    }
}

#[test]
fn stored_sample_identity_fields_are_not_the_oracle_and_fixture_records_are_not_refreshed() {
    let lines = include_str!("../../../testdata/msready/fixtures/native_sample_selection.smi");
    let selection: serde_json::Value = serde_json::from_str(include_str!(
        "../../../testdata/msready/fixtures/native_sample_selection.json"
    ))
    .unwrap();
    let rows = selection["selected_records"].as_array().unwrap();
    assert_eq!(lines.lines().count(), rows.len());
    for (line, row) in lines.lines().zip(rows) {
        let mut fields = line.split_whitespace();
        assert_eq!(fields.next(), row["smiles"].as_str());
        assert_eq!(fields.next(), row["record_id"].as_str());
        assert_eq!(row["line_sha256"].as_str().unwrap().len(), 64);
    }
    assert_eq!(selection["source_path"], "sample_input/Sample_input.smi");
    assert_eq!(
        selection["source_sha256"],
        "8ed7c503e7dfcf315e75a9d05300a4a9c6801ba5e50280cab2b8c52a37101f73"
    );
}
