use cosmolkit_model::{AtomId, MoleculeProperties};
use cosmolkit_msready::*;
use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

fn parse(text: &str) -> cosmolkit_smiles::SmilesRecord {
    parse_smiles(
        text,
        &SmilesParseParams {
            sanitize: false,
            remove_hydrogens: false,
            skip_cleanup: true,
            ..Default::default()
        },
    )
    .unwrap()
}
fn products(result: &MsReadyResult) -> &[MsReadyProduct] {
    if let MsReadyStatus::Products(p) = &result.status {
        p
    } else {
        &[]
    }
}
fn run(text: &str) -> MsReadyResult {
    standardize_smiles(text, &Default::default()).unwrap_or_else(|e| panic!("{text}: {e:?}"))
}
fn signature(result: &MsReadyResult) -> Vec<(String, String, String, i32, usize)> {
    products(result)
        .iter()
        .map(|p| {
            (
                p.inchi_key.clone(),
                p.canonical_smiles.clone(),
                p.formula.clone(),
                p.formal_charge,
                p.multiplicity,
            )
        })
        .collect()
}
fn assert_accounting(result: &MsReadyResult) {
    assert_eq!(
        result
            .atom_dispositions
            .iter()
            .map(|d| d.input_atom_row)
            .collect::<Vec<_>>(),
        (0..result.original.topology.atoms.len()).collect::<Vec<_>>()
    );
    let mut seen = std::collections::BTreeSet::new();
    for d in &result.atom_dispositions {
        assert!(seen.insert(d.input_atom_row));
        match &d.kind {
            MsReadyAtomDispositionKind::Retained {
                product_index,
                origin_index,
                output_atom_row,
            } => {
                let p = &products(result)[*product_index];
                let origin = &p.origins[*origin_index];
                assert_eq!(
                    origin.topology_mapping.atoms.old_to_new[d.input_atom_row],
                    Some(AtomId::new(*output_atom_row))
                );
                let input = &result.original.topology.atoms[d.input_atom_row];
                let output = &p.topology.atoms[*output_atom_row];
                assert_eq!(input.element(), output.element());
                assert_eq!(input.isotope(), output.isotope());
            }
            MsReadyAtomDispositionKind::FoldedHydrogen {
                product_index,
                origin_index,
                input_parent_row,
            } => {
                assert_eq!(
                    result.original.topology.atoms[d.input_atom_row].atomic_number(),
                    1
                );
                assert_eq!(
                    result.original.topology.atoms[d.input_atom_row].isotope(),
                    None
                );
                assert!(
                    products(result)[*product_index].origins[*origin_index]
                        .topology_mapping
                        .atoms
                        .old_to_new[*input_parent_row]
                        .is_some()
                );
            }
            MsReadyAtomDispositionKind::DeletedToReservoir => {
                assert_eq!(
                    result.original.topology.atoms[d.input_atom_row].atomic_number(),
                    1
                );
                assert_eq!(
                    result.original.topology.atoms[d.input_atom_row].isotope(),
                    None
                );
                assert!(result.events.iter().any(|e| e.carrier
                    == MsReadyProtonCarrier::ExplicitAtom {
                        input_atom_row: d.input_atom_row,
                        isotope: None
                    }
                    && e.kind == MsReadyEventKind::ReservoirDeprotonation));
            }
            MsReadyAtomDispositionKind::Excluded { exclusion_indices } => {
                assert!(!exclusion_indices.is_empty());
                for index in exclusion_indices {
                    assert!(
                        result.exclusions[*index]
                            .origin
                            .atom_indices
                            .contains(&d.input_atom_row)
                    );
                }
            }
        }
    }
    for p in products(result) {
        for origin in &p.origins {
            origin
                .topology_mapping
                .validate_for_counts(
                    result.original.topology.atoms.len(),
                    p.topology.atoms.len(),
                    result.original.topology.bonds.len(),
                    p.topology.bonds.len(),
                )
                .unwrap();
            for (row, source) in origin.topology_mapping.bonds.new_to_old.iter().enumerate() {
                if let Some(source) = source {
                    let input = &result.original.topology.bonds[source.index()];
                    let output = &p.topology.bonds[row];
                    // Proton/isotope transfers may change the H parent, but all
                    // heavy endpoints must remain the same stable input atoms.
                    if result.original.topology.atoms[input.begin().index()].atomic_number() != 1
                        && result.original.topology.atoms[input.end().index()].atomic_number() != 1
                    {
                        let begin = origin.topology_mapping.atoms.new_to_old
                            [output.begin().index()]
                        .unwrap();
                        let end =
                            origin.topology_mapping.atoms.new_to_old[output.end().index()].unwrap();
                        assert!(
                            (begin == input.begin() && end == input.end())
                                || (begin == input.end() && end == input.begin())
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn real_pipeline_idempotence_traversals_and_component_order_are_identity_invariant() {
    for text in [
        "CC[NH3+].[Cl-]",
        "C=C(O)C.CCO.CCO.O",
        "[NH3+]CC(=O)[O-]",
        "CC=C[NH2+]C",
        "CN(=O)=O",
        "C[N+]([O-])=[18O]",
        "[2H]OC(=C)C",
        "C[Hg]C.C[Sn](C)(C)C",
        "C[C@H](O)C(=O)O",
    ] {
        let result = run(text);
        assert_accounting(&result);
        for product in products(&result) {
            let repeated = standardize(
                &product.topology,
                &MoleculeProperties::default(),
                &Default::default(),
            )
            .unwrap();
            assert_eq!(products(&repeated).len(), 1);
            let p = &products(&repeated)[0];
            assert_eq!(
                (&p.inchi_key, &p.canonical_smiles, &p.formula),
                (
                    &product.inchi_key,
                    &product.canonical_smiles,
                    &product.formula
                )
            );
            assert_eq!(p.isotope_composition, product.isotope_composition);
            assert!((p.neutral_atom_mass - product.neutral_atom_mass).abs() < 1e-9);
            assert!(
                repeated
                    .events
                    .iter()
                    .all(|e| e.kind == MsReadyEventKind::LabelledReservoirBlocked)
            );
        }
    }
    for (left, right) in [
        (
            "CC[NH3+].O=C([O-])c1ccccc1.O",
            "O.[O-]C(=O)c1ccccc1.[NH3+]CC",
        ),
        ("CCO.CCO.CN", "NC.OCC.CCO"),
        ("CC(=O)O[Na].CN(=O)=O", "CN(=O)=O.[Na]OC(C)=O"),
    ] {
        assert_eq!(
            signature(&run(left)),
            signature(&run(right)),
            "{left} / {right}"
        );
    }
}

#[test]
fn dedup_origin_maps_align_to_representative_graph_not_fragment_row_order() {
    let result = run("CCO.OCC.CCO");
    assert_accounting(&result);
    let product = &products(&result)[0];
    assert_eq!(product.multiplicity, 3);
    assert_eq!(product.origins.len(), 3);
    assert_eq!(
        product
            .origins
            .iter()
            .map(|o| o.input_component)
            .collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
    for origin in &product.origins {
        for (output, input) in origin.topology_mapping.atoms.new_to_old.iter().enumerate() {
            let input = input.unwrap();
            assert_eq!(
                product.topology.atoms[output].element(),
                result.original.topology.atoms[input.index()].element()
            );
        }
    }
}

#[test]
fn atom_permutations_have_equal_normalized_identity_and_correct_original_ids() {
    for text in [
        "[H]OC([2H])([3H])C",
        "C=C(O)C.C[NH3+].O",
        "C[N+]([O-])=[18O]",
        "C[202Hg]C",
    ] {
        let record = parse(text);
        let count = record.topology.atoms.len();
        let order = (0..count).rev().map(AtomId::new).collect::<Vec<_>>();
        let (permuted, _) = record.topology.reordered_atoms(&order).unwrap();
        let left = standardize(&record.topology, &record.properties, &Default::default()).unwrap();
        let right = standardize(&permuted, &record.properties, &Default::default()).unwrap();
        assert_eq!(signature(&left), signature(&right));
        assert_accounting(&left);
        assert_accounting(&right);
    }
}

#[test]
fn explicit_folded_and_reservoir_hydrogens_each_have_one_original_atom_disposition() {
    for text in [
        "[H]OC([2H])([3H])C.O.[Na+]",
        "[H][NH2+]CC.[Cl-]",
        "[H][N+]([H])([H])CC.[Cl-]",
        "[H]O[H].CCO",
        "[2H][NH2+]CC(=O)[O-]",
    ] {
        let result = run(text);
        assert_accounting(&result);
        if text.starts_with("[H][N+]") {
            assert!(
                result
                    .atom_dispositions
                    .iter()
                    .any(|d| d.kind == MsReadyAtomDispositionKind::DeletedToReservoir)
            );
        }
    }
    let source = parse("[H]OC.C[Hg]C");
    let snapshot = source.clone();
    let result = standardize(&source.topology, &source.properties, &Default::default()).unwrap();
    assert_eq!(source, snapshot);
    assert_eq!(result.original.topology, snapshot.topology);
    assert_accounting(&result);
}

#[test]
fn tracked_and_explicit_isotope_products_match_with_formula_mass_and_parent_provenance() {
    for isotope in [2, 3] {
        let mut source = parse("C=C(O)C");
        source.topology.atoms[2].set_tracked_isotopic_hydrogens(vec![isotope]);
        let snapshot = source.clone();
        let tracked =
            standardize(&source.topology, &source.properties, &Default::default()).unwrap();
        let explicit = run(&format!("[{isotope}H]OC(=C)C"));
        assert_eq!(signature(&tracked), signature(&explicit));
        assert_eq!(source, snapshot);
        assert_accounting(&tracked);
        assert_accounting(&explicit);
        assert_eq!(
            products(&tracked)[0].isotope_composition,
            products(&explicit)[0].isotope_composition
        );
        assert!(
            (products(&tracked)[0].neutral_atom_mass - products(&explicit)[0].neutral_atom_mass)
                .abs()
                < 1e-9
        );
        assert!(
            products(&tracked)[0].origins[0]
                .topology_mapping
                .atoms
                .new_to_old
                .iter()
                .any(Option::is_none)
        );
    }
    let result = run("[2H][N+]([2H])([3H])C");
    let product = &products(&result)[0];
    assert_eq!(product.formal_charge, 1);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code == "labelled_hydrogen_not_removed_from_reservoir")
    );
    assert!(
        result
            .events
            .iter()
            .all(|e| e.kind == MsReadyEventKind::LabelledReservoirBlocked)
    );
    assert_accounting(&result);
}

#[test]
fn shared_input_budget_covers_components_disconnections_and_group_charge_rematches() {
    for (text, required_edits) in [
        ("C=C(O)C.C[NH3+].O", 2),
        ("CC=C[NH2+]C", 2),
        ("CC(=O)O[Na].C=C(O)C", 3),
    ] {
        let params = MsReadyParams {
            max_edits: 1,
            ..Default::default()
        };
        let failure = standardize_smiles(text, &params).unwrap_err();
        assert_eq!(
            failure.kind,
            MsReadyFailureKind::TransformationLimit,
            "{text}: {failure:?}"
        );
        assert!(failure.original.is_some());
        assert_eq!(
            failure.evidence.input_link.source_smiles.as_deref(),
            Some(text)
        );
        assert_eq!(
            failure
                .evidence
                .events
                .iter()
                .filter(|e| e.kind != MsReadyEventKind::ParserRadicalInitialization
                    && e.kind != MsReadyEventKind::LabelledReservoirBlocked)
                .count(),
            1
        );
        assert!(failure.evidence.location.component_index.is_some());
        assert!(failure.evidence.location.rule_id.is_some());
        assert!(
            standardize_smiles(
                text,
                &MsReadyParams {
                    max_edits: required_edits,
                    ..Default::default()
                }
            )
            .is_ok()
        );
    }
    let failure = standardize_smiles(
        "[2H][NH+]([3H])CC(=O)[O-]",
        &MsReadyParams {
            max_matches_per_sweep: 1,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(failure.kind, MsReadyFailureKind::TransformationLimit);
    let cuts = standardize_smiles(
        "CO[Li].CO[Na]",
        &MsReadyParams {
            max_matches_per_sweep: 1,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(cuts.kind, MsReadyFailureKind::TransformationLimit);
    assert_eq!(cuts.stage, MsReadyStage::Disconnection);
    assert!(
        cuts.evidence.events.is_empty(),
        "all candidate attempts must be bounded before constructing the first product"
    );
    let one = standardize_smiles(
        "C=C(O)C",
        &MsReadyParams {
            max_edits: 1,
            max_matches_per_sweep: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(products(&one)[0].inchi_key, "CSCPPACGZOOCGX-UHFFFAOYSA-N");
}

#[test]
fn separate_enol_components_start_separate_match_sweeps_but_share_input_edits() {
    for text in ["C=C(O)C.C=C(O)CC", "C=C(O)CC.C=C(O)C"] {
        let params = MsReadyParams {
            max_matches_per_sweep: 1,
            max_edits: 2,
            ..Default::default()
        };
        let result = standardize_smiles(text, &params).unwrap();
        assert_eq!(products(&result).len(), 2);
        assert_eq!(signature(&result), signature(&run("CC(=O)C.CC(=O)CC")));
        assert_eq!(
            result
                .events
                .iter()
                .filter(|e| e.kind == MsReadyEventKind::GroupRewrite)
                .count(),
            2
        );
        assert_accounting(&result);
        let source = parse(text);
        assert_eq!(
            signature(&standardize(&source.topology, &source.properties, &params).unwrap()),
            signature(&result)
        );
        let failure = standardize_smiles(
            text,
            &MsReadyParams {
                max_edits: 1,
                ..params
            },
        )
        .unwrap_err();
        assert_eq!(failure.kind, MsReadyFailureKind::TransformationLimit);
        assert_eq!(
            failure
                .evidence
                .events
                .iter()
                .filter(|e| e.kind == MsReadyEventKind::GroupRewrite)
                .count(),
            1,
            "edit counter must never reset across components"
        );
    }
}

#[test]
fn alkali_observation_and_actual_collection_share_a_finite_allowance_without_double_charging() {
    let params = MsReadyParams {
        max_matches_per_sweep: 1,
        max_edits: 2,
        ..Default::default()
    };
    for ligand in ["CO[Na]", "C[O-][Na+]"] {
        let one = standardize_smiles(ligand, &params).unwrap();
        assert_eq!(signature(&one), signature(&run("CO")));
        assert_eq!(
            one.events
                .iter()
                .filter(|e| e.kind == MsReadyEventKind::AlkaliDisconnection)
                .count(),
            1,
            "observation must not charge an actual match twice"
        );
        for text in [
            std::iter::repeat_n(ligand, 64)
                .collect::<Vec<_>>()
                .join("."),
            "[Li+][O-]C(=O)CC(=O)[O-][Na+]".into(),
        ] {
            let source = parse(&text);
            for failure in [
                standardize_smiles(&text, &params).unwrap_err(),
                standardize(&source.topology, &source.properties, &params).unwrap_err(),
            ] {
                assert_eq!(failure.kind, MsReadyFailureKind::TransformationLimit);
                assert_eq!(failure.stage, MsReadyStage::Disconnection);
                assert_eq!(
                    failure.evidence.location.rule_id.as_deref(),
                    Some("alkali_oxygen_disconnection")
                );
                assert_eq!(failure.evidence.location.atom_indices.len(), 2);
                assert_eq!(failure.evidence.location.bond_indices.len(), 1);
                assert!(
                    failure.evidence.events.is_empty(),
                    "bounded raw observation/actual collection must fail before chemical edits"
                );
                assert_eq!(failure.original.unwrap().topology, source.topology);
            }
        }
    }
}

#[test]
fn full_profile_repeats_charge_then_group_when_proton_edit_exposes_an_earlier_rule() {
    let result = run("CC=C[NH2+]C");
    assert_eq!(signature(&result), signature(&run("CCC=NC")));
    let kinds = result
        .events
        .iter()
        .map(|e| e.kind.clone())
        .collect::<Vec<_>>();
    assert_eq!(
        kinds,
        vec![
            MsReadyEventKind::ReservoirDeprotonation,
            MsReadyEventKind::GroupRewrite
        ]
    );
    assert_eq!(result.events[1].rule_id, "enamine_to_imine");
    assert_accounting(&result);
}

#[test]
fn formula_charge_and_isotope_accounting_distinguishes_internal_transfers_from_reservoirs() {
    for text in [
        "[NH3+]CC(=O)[O-]",
        "[2H][NH2+]CC(=O)[O-]",
        "CC[NH3+].CC(=O)[O-]",
        "C=C(O)C",
        "CC(=O)O[Na]",
    ] {
        let result = run(text);
        assert_accounting(&result);
        for event in &result.events {
            assert_eq!(
                event.total_charge_after - event.total_charge_before,
                event.delta_charge
            );
            if matches!(
                event.kind,
                MsReadyEventKind::GroupRewrite
                    | MsReadyEventKind::InternalProtonTransfer
                    | MsReadyEventKind::AlkaliDisconnection
            ) {
                assert_eq!(event.formula_before, event.formula_after);
                assert_eq!(event.delta_hydrogens, 0);
                assert_eq!(event.delta_charge, 0);
            }
            if matches!(
                event.kind,
                MsReadyEventKind::ReservoirDeprotonation | MsReadyEventKind::ReservoirProtonation
            ) {
                assert_eq!(event.delta_hydrogens, event.delta_charge);
                assert_eq!(
                    event
                        .atom_charge_deltas
                        .iter()
                        .map(|(_, q)| i32::from(*q))
                        .sum::<i32>(),
                    event.delta_charge
                );
            }
            event
                .topology_mapping
                .validate_for_counts(
                    result.original.topology.atoms.len(),
                    event.topology_mapping.atoms.new_to_old.len(),
                    result.original.topology.bonds.len(),
                    event.topology_mapping.bonds.new_to_old.len(),
                )
                .unwrap();
        }
        if text.starts_with("[NH3+]") || text.starts_with("[2H]") {
            assert!(
                result
                    .events
                    .iter()
                    .any(|e| e.kind == MsReadyEventKind::InternalProtonTransfer)
            );
            assert!(result.events.iter().all(|e| !matches!(
                e.kind,
                MsReadyEventKind::ReservoirDeprotonation | MsReadyEventKind::ReservoirProtonation
            )));
        }
    }
}

#[test]
fn failures_retain_whole_original_precise_owner_context_and_all_prior_exclusions_events() {
    let record = parse("O.C[NH3+].C=C(O)C");
    let snapshot = record.clone();
    let params = MsReadyParams {
        record_id: Some("record42".into()),
        original_smiles: Some("archival metadata, not actual parsed text".into()),
        max_edits: 1,
        ..Default::default()
    };
    let failure = standardize(&record.topology, &record.properties, &params).unwrap_err();
    assert_eq!(record, snapshot);
    assert_eq!(failure.original.unwrap().topology, snapshot.topology);
    assert!(!failure.evidence.exclusions.is_empty());
    assert!(!failure.evidence.events.is_empty());
    assert_eq!(
        failure.evidence.input_link.record_id.as_deref(),
        Some("record42")
    );
    let failure = standardize_smiles("O.CN(=O)(=O)C", &Default::default()).unwrap_err();
    assert_eq!(failure.kind, MsReadyFailureKind::InvalidValence);
    assert!(!failure.evidence.exclusions.is_empty());
    assert_eq!(failure.evidence.location.input_component, Some(1));
    assert_eq!(failure.evidence.location.atom_indices, vec![2]);
    assert!(failure.evidence.location.owner_category.is_some());
    let parse_failure = standardize_smiles("C>>O", &params).unwrap_err();
    assert_eq!(parse_failure.kind, MsReadyFailureKind::InvalidInput);
    assert_eq!(
        parse_failure.evidence.input_link.source_smiles.as_deref(),
        Some("C>>O")
    );
    assert_eq!(
        parse_failure.evidence.input_link.original_smiles,
        params.original_smiles
    );
    assert!(parse_failure.evidence.location.owner_category.is_some());
}

#[test]
fn graph_level_unsupported_annotations_and_metal_states_fail_without_changing_input() {
    for change in 0..5 {
        let mut source = parse("C[Hg]C");
        match change {
            0 => source.topology.atoms[1].set_radical_electrons(1),
            1 => source.topology.bonds[0].set_order(cosmolkit_types::BondOrder::Dative),
            2 => source.topology.bonds[0].set_order(cosmolkit_types::BondOrder::Double),
            3 => source.topology.atoms[1].set_explicit_hydrogens(1),
            _ => {
                let bond = &source.topology.bonds[0];
                source.topology.bonds[0] = cosmolkit_model::Bond::from_spec(
                    bond.id(),
                    cosmolkit_model::BondSpec::new(bond.begin(), bond.end(), bond.order())
                        .with_query(cosmolkit_model::QueryNode::Predicate(
                            cosmolkit_model::BondQueryPredicate::Any,
                        )),
                );
            }
        }
        let snapshot = source.clone();
        assert!(standardize(&source.topology, &source.properties, &Default::default()).is_err());
        assert_eq!(source, snapshot);
    }
    for invalid_text in ["CN([O-])=[18O]", "[2H][3H][NH+]CC(=O)[O-]"] {
        let invalid = standardize_smiles(invalid_text, &Default::default()).unwrap_err();
        assert_eq!(invalid.kind, MsReadyFailureKind::InvalidValence);
    }
}
