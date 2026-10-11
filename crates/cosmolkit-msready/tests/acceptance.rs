use cosmolkit_msready::*;

fn products(result: &MsReadyResult) -> &[MsReadyProduct] {
    match &result.status {
        MsReadyStatus::Products(products) => products,
        MsReadyStatus::NoEligibleComponent => &[],
    }
}
fn checked(text: &str) -> MsReadyResult {
    let result = standardize_smiles(text, &MsReadyParams::default())
        .unwrap_or_else(|e| panic!("{text}: {e:?}"));
    assert_eq!(result.input_link.source_smiles.as_deref(), Some(text));
    assert_eq!(
        result.atom_dispositions.len(),
        result.original.topology.atoms.len()
    );
    for product in products(&result) {
        assert_eq!(product.inchi_key.len(), 27);
        assert!(product.inchi.starts_with("InChI=1S/"));
        assert!(!product.canonical_smiles.contains(['@', '/', '\\']));
        assert!(
            !product
                .inchi
                .split('/')
                .any(|layer| matches!(layer.as_bytes().first(), Some(b'b' | b't' | b'm' | b's')))
        );
        assert_eq!(product.multiplicity, product.origins.len());
        product.topology.validate().unwrap();
        assert!(product.neutral_atom_mass.is_finite() && product.neutral_atom_mass > 0.0);
        assert_eq!(
            product.formal_charge,
            product
                .topology
                .atoms
                .iter()
                .map(|a| i32::from(a.formal_charge()))
                .sum::<i32>()
        );
        assert!(product.topology.stereo_groups.is_empty());
        for atom in &product.topology.atoms {
            assert_eq!(atom.chiral_tag(), cosmolkit_types::ChiralTag::Unspecified);
            assert_eq!(atom.chiral_permutation(), None);
            assert!(!atom.unknown_stereo());
            assert_eq!(atom.mol_parity(), None);
            assert_eq!(atom.mol_inversion_flag(), None);
            assert_eq!(atom.atom_map(), None);
            assert!(atom.tracked_isotopic_hydrogens().is_empty());
        }
        for bond in &product.topology.bonds {
            assert_eq!(bond.stereo(), cosmolkit_types::BondStereo::None);
            assert!(bond.stereo_atom_references().is_empty());
            assert_eq!(bond.direction(), cosmolkit_types::BondDirection::None);
            assert!(!bond.unknown_stereo());
        }
    }
    result
}
fn identity(text: &str) -> Vec<(String, String)> {
    products(&checked(text))
        .iter()
        .map(|p| (p.inchi_key.clone(), p.canonical_smiles.clone()))
        .collect()
}
fn only_product(result: &MsReadyResult) -> &MsReadyProduct {
    assert_eq!(products(result).len(), 1, "required retained product");
    &products(result)[0]
}

fn assert_composition(product: &MsReadyProduct, expected: &[(u8, u32)]) {
    let mut composition = product
        .isotope_composition
        .iter()
        .map(|entry| {
            assert_eq!(entry.isotope, None);
            (entry.atomic_number, entry.count)
        })
        .collect::<Vec<_>>();
    composition.sort();
    assert_eq!(composition, expected);
    assert_eq!(product.formal_charge, 0);
}

fn element_bonds(product: &MsReadyProduct) -> Vec<(u8, u8, cosmolkit_types::BondOrder)> {
    let mut bonds = product
        .topology
        .bonds
        .iter()
        .map(|bond| {
            let a = product.topology.atoms[bond.begin().index()].atomic_number();
            let b = product.topology.atoms[bond.end().index()].atomic_number();
            (a.min(b), a.max(b), bond.order())
        })
        .collect::<Vec<_>>();
    bonds.sort_by_key(|bond| (bond.0, bond.1));
    bonds
}

#[test]
fn all_twenty_section13_published_concept_cases_have_real_products_and_graph_assertions() {
    // Representative reference strings are graph expectations, not historical
    // byte or CSV-InChI oracles. This test explicitly covers all twenty rows.
    for (input, expected) in [
        ("CC[NH3+].[Cl-]", "CCN"),
        ("CC(=O)[O-].[Na+]", "CC(=O)O"),
        ("CC(=O)O[Na]", "CC(=O)O"),
        ("C[N+](C)(C)C.[Cl-]", "C[N+](C)(C)C"),
        ("CC[NH3+].O=C([O-])c1ccccc1", "CCN.O=C(O)c1ccccc1"),
        ("CCO.O", "CCO"),
        ("CCO.CCO", "CCO"),
        ("[Na+].[Cl-]", "[Na+].[Cl-]"),
        ("[NH3+]CC(=O)[O-]", "NCC(=O)O"),
        ("C=C(O)C", "CC(=O)C"),
        ("C[C@H](O)C(=O)O", "CC(O)C(=O)O"),
        ("C/C=C/C", "CC=CC"),
        ("C[13CH2]O", "C[13CH2]O"),
        ("C[Hg+].[Cl-]", "C[Hg+]"),
        ("*CCO", "*CCO"),
        ("[CH3]", "[CH3]"),
    ] {
        if input == "*CCO" {
            let failure = standardize_smiles(input, &Default::default()).unwrap_err();
            assert_eq!(failure.kind, MsReadyFailureKind::UnsupportedChemistry);
            continue;
        }
        assert_eq!(identity(input), identity(expected), "Section13 {input}");
    }
    // These four unchanged-chemistry rows need independent graph expectations,
    // not identity(input) == identity(the same input), which passes on exclusion.
    let mixture = checked("CCO.CN");
    assert_eq!(products(&mixture).len(), 2);
    assert!(mixture.exclusions.is_empty());
    let ethanol = products(&mixture)
        .iter()
        .find(|p| p.formula == "C2H6O")
        .unwrap();
    assert_composition(ethanol, &[(1, 6), (6, 2), (8, 1)]);
    assert_eq!(ethanol.topology.atoms.len(), 3);
    assert_eq!(ethanol.topology.bonds.len(), 2);
    assert_eq!(
        element_bonds(ethanol),
        vec![
            (6, 6, cosmolkit_types::BondOrder::Single),
            (6, 8, cosmolkit_types::BondOrder::Single)
        ]
    );
    let amine = products(&mixture)
        .iter()
        .find(|p| p.formula == "CH5N")
        .unwrap();
    assert_composition(amine, &[(1, 5), (6, 1), (7, 1)]);
    assert_eq!(amine.topology.atoms.len(), 2);
    assert_eq!(
        element_bonds(amine),
        vec![(6, 7, cosmolkit_types::BondOrder::Single)]
    );

    let covalent_chlorine = checked("CCl");
    let chlorine = only_product(&covalent_chlorine);
    assert!(covalent_chlorine.exclusions.is_empty());
    assert_eq!(chlorine.formula, "CH3Cl");
    assert_composition(chlorine, &[(1, 3), (6, 1), (17, 1)]);
    assert_eq!(chlorine.topology.atoms.len(), 2);
    assert_eq!(
        element_bonds(chlorine),
        vec![(6, 17, cosmolkit_types::BondOrder::Single)]
    );
    assert!(
        chlorine
            .topology
            .atoms
            .iter()
            .all(|a| a.formal_charge() == 0)
    );

    let connected_phosphate = checked("COP(=O)(O)O");
    let phosphate = only_product(&connected_phosphate);
    assert!(connected_phosphate.exclusions.is_empty());
    assert_eq!(phosphate.formula, "CH5O4P");
    assert_composition(phosphate, &[(1, 5), (6, 1), (8, 4), (15, 1)]);
    assert_eq!(phosphate.topology.atoms.len(), 6);
    assert_eq!(phosphate.topology.bonds.len(), 5);
    let phosphorus = phosphate
        .topology
        .atoms
        .iter()
        .find(|a| a.atomic_number() == 15)
        .unwrap();
    let p_neighbors = phosphate
        .topology
        .adjacency
        .neighbors_of(phosphorus.id().index());
    assert_eq!(p_neighbors.len(), 4);
    assert!(
        p_neighbors
            .iter()
            .all(|n| phosphate.topology.atoms[n.atom_index].atomic_number() == 8)
    );
    assert_eq!(
        p_neighbors
            .iter()
            .filter(|n| phosphate.topology.bonds[n.bond.index()].order()
                == cosmolkit_types::BondOrder::Double)
            .count(),
        1
    );
    assert_eq!(
        p_neighbors
            .iter()
            .filter(|n| phosphate.topology.bonds[n.bond.index()].order()
                == cosmolkit_types::BondOrder::Single)
            .count(),
        3
    );
    let carbon = phosphate
        .topology
        .atoms
        .iter()
        .find(|a| a.atomic_number() == 6)
        .unwrap();
    let c_neighbors = phosphate
        .topology
        .adjacency
        .neighbors_of(carbon.id().index());
    assert_eq!(c_neighbors.len(), 1);
    let bridge = c_neighbors[0].atom_index;
    assert_eq!(phosphate.topology.atoms[bridge].atomic_number(), 8);
    assert_eq!(
        phosphate.topology.bonds[c_neighbors[0].bond.index()].order(),
        cosmolkit_types::BondOrder::Single
    );
    let bridge_neighbors = phosphate.topology.adjacency.neighbors_of(bridge);
    assert_eq!(bridge_neighbors.len(), 2);
    assert!(
        bridge_neighbors
            .iter()
            .any(|n| n.atom_index == phosphorus.id().index()
                && phosphate.topology.bonds[n.bond.index()].order()
                    == cosmolkit_types::BondOrder::Single)
    );
    assert!(
        p_neighbors
            .iter()
            .filter(|n| n.atom_index != bridge)
            .all(|n| phosphate
                .topology
                .adjacency
                .neighbors_of(n.atom_index)
                .len()
                == 1)
    );
    assert!(
        phosphate
            .topology
            .atoms
            .iter()
            .all(|a| a.formal_charge() == 0)
    );

    let charged_nitro = checked("C[N+](=O)[O-]");
    let nitro = only_product(&charged_nitro);
    assert!(charged_nitro.exclusions.is_empty());
    assert_eq!(nitro.formula, "CH3NO2");
    assert_composition(nitro, &[(1, 3), (6, 1), (7, 1), (8, 2)]);
    assert_eq!(nitro.topology.atoms.len(), 4);
    assert_eq!(nitro.topology.bonds.len(), 3);
    let nitrogen = nitro
        .topology
        .atoms
        .iter()
        .find(|a| a.atomic_number() == 7)
        .unwrap();
    assert_eq!(nitrogen.formal_charge(), 1);
    let n_neighbors = nitro.topology.adjacency.neighbors_of(nitrogen.id().index());
    assert_eq!(n_neighbors.len(), 3);
    let mut terminal_states = n_neighbors
        .iter()
        .map(|n| {
            let atom = &nitro.topology.atoms[n.atom_index];
            assert_eq!(nitro.topology.adjacency.neighbors_of(n.atom_index).len(), 1);
            (
                atom.atomic_number(),
                atom.formal_charge(),
                nitro.topology.bonds[n.bond.index()].order(),
            )
        })
        .collect::<Vec<_>>();
    terminal_states.sort_by_key(|state| (state.0, state.1));
    assert_eq!(
        terminal_states,
        vec![
            (6, 0, cosmolkit_types::BondOrder::Single),
            (8, -1, cosmolkit_types::BondOrder::Single),
            (8, 0, cosmolkit_types::BondOrder::Double)
        ]
    );

    let repeats = checked("CCO.CCO");
    assert_eq!(products(&repeats)[0].multiplicity, 2);
    let mixture = checked("CC[NH3+].O=C([O-])c1ccccc1");
    assert_eq!(products(&mixture).len(), 2);
    assert_eq!(
        products(&checked("C/C=C\\C"))[0].inchi_key,
        products(&checked("C/C=C/C"))[0].inchi_key
    );
    assert!(matches!(
        checked("[Na+].[Cl-]").status,
        MsReadyStatus::NoEligibleComponent
    ));
    let radical = checked("[CH3]");
    assert_eq!(
        products(&radical)[0].topology.atoms[0].radical_electrons(),
        1
    );
    assert!(
        radical
            .warnings
            .iter()
            .any(|w| w.code == "stable_radical_preserved")
    );
    assert!(
        radical
            .events
            .iter()
            .any(|e| e.kind == MsReadyEventKind::ParserRadicalInitialization)
    );
    let isotope = checked("C[13CH2]O");
    assert_eq!(products(&isotope)[0].formula, "C[13C]H6O");
    assert_ne!(
        products(&isotope)[0].inchi_key,
        products(&checked("CCO"))[0].inchi_key
    );
}

#[test]
fn every_alkali_element_charge_pair_ligand_and_hydroxide_contract_is_integrated() {
    for metal in ["Li", "Na", "K"] {
        for input in [
            format!("CC(=O)O[{metal}]"),
            format!("CC(=O)[O-][{metal}+]"),
            format!("CC(=O)[O-].[{metal}+]"),
        ] {
            let result = checked(&input);
            assert_eq!(identity(&input), identity("CC(=O)O"));
            assert!(result.exclusions.iter().any(|e| {
                e.origin
                    .atom_indices
                    .iter()
                    .any(|row| result.original.topology.atoms[*row].element().symbol() == metal)
            }));
            for event in &result.events {
                if event.kind == MsReadyEventKind::AlkaliDisconnection {
                    assert_eq!(event.formula_before, event.formula_after);
                    assert_eq!(event.total_charge_before, event.total_charge_after);
                }
            }
        }
        assert_eq!(identity(&format!("CO[{metal}]")), identity("CO"));
        for input in [format!("[{metal}]O"), format!("[OH-][{metal}+]")] {
            let result = checked(&input);
            assert!(matches!(result.status, MsReadyStatus::NoEligibleComponent));
            assert!(
                result
                    .events
                    .iter()
                    .any(|e| e.kind == MsReadyEventKind::AlkaliDisconnection)
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
        }
        for input in [
            format!("CC(=O)O[{metal}+]"),
            format!("CC(=O)[O-][{metal}]"),
            format!("C[{metal}]"),
            format!("N[{metal}]"),
            format!("O=[{metal}]"),
            format!("C[O]([{metal}])C"),
            format!("[H][{metal}]O"),
        ] {
            assert!(
                standardize_smiles(&input, &Default::default()).is_err(),
                "must reject {input}"
            );
        }
    }
}

#[test]
fn retained_metal_connectivity_and_labelled_states_are_actual_products_not_identifier_roundtrips() {
    for text in [
        "C[Hg+].[Cl-]",
        "C[Hg]C",
        "C[Sn](C)(C)C",
        "C[202Hg]C",
        "C[120Sn](C)(C)C",
        "[23Na+]",
        "[7Li+]",
        "[41K+]",
    ] {
        let result = checked(text);
        let product = &products(&result)[0];
        let metal = product
            .topology
            .atoms
            .iter()
            .find(|a| matches!(a.atomic_number(), 3 | 11 | 19 | 50 | 80))
            .unwrap();
        let degree = product
            .topology
            .adjacency
            .neighbors_of(metal.id().index())
            .len();
        assert_eq!(
            degree,
            match metal.atomic_number() {
                50 => 4,
                80 if metal.formal_charge() == 1 => 1,
                80 => 2,
                _ => 0,
            }
        );
        for neighbor in product.topology.adjacency.neighbors_of(metal.id().index()) {
            assert_eq!(
                product.topology.atoms[neighbor.atom_index].atomic_number(),
                6
            );
            assert_eq!(
                product.topology.bonds[neighbor.bond.index()].order(),
                cosmolkit_types::BondOrder::Single
            );
        }
        for origin in &product.origins {
            origin
                .topology_mapping
                .validate_for_counts(
                    result.original.topology.atoms.len(),
                    product.topology.atoms.len(),
                    result.original.topology.bonds.len(),
                    product.topology.bonds.len(),
                )
                .unwrap();
        }
    }
    for text in [
        "C[Hg]",
        "C[Hg+2]",
        "C[Hg]Cl",
        "C[Sn](C)C",
        "[H][Sn](C)(C)C",
        "c1ccccc1[Hg+]",
        "C[Fe]C",
        "C[Cu]C",
        "C[Zn]C",
        "CO[Na]OC",
    ] {
        let failure = standardize_smiles(text, &Default::default()).unwrap_err();
        assert!(
            matches!(
                failure.kind,
                MsReadyFailureKind::UnsupportedChemistry | MsReadyFailureKind::InvalidValence
            ),
            "{text}: {failure:?}"
        );
        assert!(failure.original.is_some());
    }
}

#[test]
fn nitrate_templates_are_comparison_only_with_exact_source_provenance_and_no_fake_events() {
    for (text, template, formula, charge) in [
        ("ON(=O)O", "83_raw_h2no3", "H2NO3", 0),
        ("[OH][N+](=O)[O-]", "83_nitric_acid", "HNO3", 0),
        ("[O-][N+](=O)[O-]", "83_nitrate", "NO3-", -1),
    ] {
        for input in [
            text.to_owned(),
            format!("CCO.{text}"),
            format!("{text}.CCO"),
        ] {
            let result = checked(&input);
            let ex = result
                .exclusions
                .iter()
                .find(|e| e.source_id == "83")
                .unwrap();
            assert_eq!(ex.template_id.as_deref(), Some(template));
            assert_eq!(ex.source_sha256, EXCLUSION_SHA256);
            assert_eq!(ex.source_field, "CANONICAL_SMILES");
            assert!(ex.comparison_only);
            assert_eq!(ex.matched_formula, formula);
            assert_eq!(ex.matched_formal_charge, charge);
            assert_eq!(
                ex.adaptation_id.as_deref(),
                if template == "83_raw_h2no3" {
                    None
                } else {
                    Some(NITRATE_ADAPTATION_ID)
                }
            );
            assert!(result.events.iter().all(|e| {
                !e.location
                    .atom_indices
                    .iter()
                    .any(|row| ex.origin.atom_indices.contains(row))
            }));
        }
    }
    for text in [
        "[18OH][N+](=O)[O-]",
        "[18O-][N+](=O)[O-]",
        "O=N[O-]",
        "CO[N+](=O)[O-]",
        "C[N+](=O)[O-]",
        "[OH][N+](=O)[OH]",
    ] {
        let result = checked(text);
        assert!(result.exclusions.iter().all(|e| e.source_id != "83"));
    }
    for text in [
        "[18OH]N(=O)O",
        "O[N](=O)[O-]",
        "[15N](O)(=O)O",
        "[NH](O)(=O)O",
    ] {
        let failure = standardize_smiles(text, &Default::default()).unwrap_err();
        assert_eq!(
            failure.kind,
            MsReadyFailureKind::InvalidValence,
            "{text}: {failure:?}"
        );
    }
}

#[test]
fn labelled_alkali_hydroxide_and_explicit_h_nitrate_preserve_complete_atom_origins() {
    for (metal, isotope) in [("Li", 7), ("Na", 23), ("K", 41)] {
        let labelled = checked(&format!("CC(=O)O[{isotope}{metal}]"));
        assert_eq!(products(&labelled).len(), 2);
        assert!(
            products(&labelled)
                .iter()
                .any(|p| p.topology.atoms.len() == 1
                    && p.formal_charge == 1
                    && p.topology.atoms[0].isotope() == Some(isotope))
        );
        assert!(
            labelled.exclusions.is_empty(),
            "isotope labels must not be erased for metal exclusions"
        );
        let deuterated = checked(&format!("[{metal}]O[2H]"));
        let product = &products(&deuterated)[0];
        assert_eq!(product.topology.atoms.len(), 2);
        assert_eq!(product.formal_charge, -1);
        assert!(
            product
                .isotope_composition
                .iter()
                .any(|e| e.atomic_number == 1 && e.isotope == Some(2) && e.count == 1)
        );
        assert!(
            deuterated
                .exclusions
                .iter()
                .all(|e| e.origin.atom_indices == vec![0])
        );
    }
    for text in [
        "[H]ON(=O)O[H]",
        "O=N(O)O",
        "[H]O[N+](=O)[O-]",
        "[O-][N+]([O-])=O",
    ] {
        let result = checked(text);
        assert!(matches!(result.status, MsReadyStatus::NoEligibleComponent));
        assert!(result.events.is_empty());
        assert!(result.exclusions.iter().all(|e| e.source_id == "83"));
        assert!(
            result
                .atom_dispositions
                .iter()
                .all(|d| matches!(d.kind, MsReadyAtomDispositionKind::Excluded { .. }))
        );
        assert!(
            result
                .exclusions
                .iter()
                .all(|e| e.origin.atom_indices.len() == result.original.topology.atoms.len())
        );
    }
}

#[test]
fn text_radical_initialization_preserves_exact_raw_exclusions_and_explicit_cx_radicals() {
    for text in ["[Na]", "[Li]", "[K]", "[Ag]", "[Zn]"] {
        assert!(matches!(
            checked(text).status,
            MsReadyStatus::NoEligibleComponent
        ));
    }
    assert!(matches!(
        checked("Cl").status,
        MsReadyStatus::NoEligibleComponent
    ));
    let chlorine = checked("[Cl]");
    assert_eq!(
        products(&chlorine)[0].topology.atoms[0].radical_electrons(),
        1
    );
    assert!(chlorine.exclusions.is_empty());
    for text in ["[CH3] |^1:0|", "c1ccccc1C |^1:6|"] {
        let result = checked(text);
        assert!(
            products(&result).iter().any(|p| p
                .topology
                .atoms
                .iter()
                .any(|a| a.radical_electrons() == 1))
        );
    }
    for metal in ["Li", "Na", "K"] {
        for text in [
            format!("[CH2]O[{metal}]"),
            format!("C[O][{metal}] |^1:1|"),
            format!("[CH2][O-][{metal}+]"),
        ] {
            let failure = standardize_smiles(&text, &Default::default()).unwrap_err();
            assert!(
                failure
                    .evidence
                    .events
                    .iter()
                    .all(|e| e.kind != MsReadyEventKind::AlkaliDisconnection),
                "radical ligand must not be edited: {text}"
            );
        }
        // Neutral O with two single bonds is not a radical under the owner;
        // the earlier C[O][Na] negative expectation was chemically incorrect.
        assert_eq!(identity(&format!("C[O][{metal}]")), identity("CO"));
    }
    for text in [
        "CC(=O)[O-][Li+].CC(=O)[O-][Na+].CO[K]",
        "CO[Li].CC(=O)[O-][K+].CC(=O)[O-][K+]",
    ] {
        let result = checked(text);
        assert_eq!(
            result
                .events
                .iter()
                .filter(|e| e.kind == MsReadyEventKind::AlkaliDisconnection)
                .count(),
            3
        );
        assert!(
            result
                .original
                .topology
                .atoms
                .iter()
                .all(|a| a.radical_electrons() == 0)
        );
        assert!(
            result
                .events
                .iter()
                .all(|e| e.kind != MsReadyEventKind::ParserRadicalInitialization)
        );
    }
}
