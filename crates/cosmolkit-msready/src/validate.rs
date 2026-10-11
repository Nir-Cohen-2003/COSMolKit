use cosmolkit_core::{SanitizeOperations, SanitizeParams, sanitize_topology};
use cosmolkit_model::{BondOrder, Element, MoleculeProperties, TopologyBlock};

use crate::{
    MsReadyFailure, MsReadyFailureKind as Kind, MsReadyParams, MsReadyStage as Stage, profile,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValidationClass {
    Ordinary,
    RawNitrateExclusionOnly,
    RecognizedAlkaliRepair,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RetainedMetalKind {
    HgCationMonovalent,
    HgNeutralDivalent,
    SnNeutralTetravalent,
    AlkaliSpectator,
}

fn failure(kind: Kind, stage: Stage, context: impl Into<String>) -> MsReadyFailure {
    MsReadyFailure {
        evidence: Default::default(),
        kind,
        stage,
        context: context.into(),
        original: None,
    }
}

pub fn validate_syntax(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &MsReadyParams,
) -> Result<(), MsReadyFailure> {
    validate_syntax_inner(topology, properties, params).map_err(|mut error| {
        error.original = Some(crate::MsReadyOriginal {
            topology: topology.clone(),
            properties: properties.clone(),
        });
        error
    })
}

fn validate_syntax_inner(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &MsReadyParams,
) -> Result<(), MsReadyFailure> {
    if params.max_edits == 0 || params.max_matches_per_sweep == 0 {
        return Err(failure(
            Kind::InvalidInput,
            Stage::Input,
            "transformation limits must be nonzero",
        ));
    }
    if topology.atoms.is_empty() {
        return Err(failure(
            Kind::InvalidInput,
            Stage::Input,
            "empty molecule topology",
        ));
    }
    topology
        .validate()
        .map_err(|error| failure(Kind::InvalidInput, Stage::Input, error.to_string()))?;
    drop(
        properties
            .property_records(true, true)
            .map_err(|error| failure(Kind::InvalidInput, Stage::Input, error.to_string()))?,
    );
    for (i, atom) in topology.atoms.iter().enumerate() {
        if atom.atomic_number() == 0 {
            return Err(failure(
                Kind::UnsupportedChemistry,
                Stage::Input,
                format!("wildcard/query atom at {i}"),
            )
            .with_atoms(vec![i]));
        }
        if atom.template_attachment_order().is_some() || atom.pdb_residue_info().is_some() {
            return Err(failure(
                Kind::UnsupportedChemistry,
                Stage::Input,
                format!("unsupported attachment or atom annotation at {i}"),
            )
            .with_atoms(vec![i]));
        }
    }
    for (i, bond) in topology.bonds.iter().enumerate() {
        if bond.query().is_some() {
            return Err(failure(
                Kind::UnsupportedChemistry,
                Stage::Input,
                format!("query bond at {i}"),
            )
            .with_bonds(vec![i]));
        }
        if !matches!(
            bond.order(),
            BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
        ) {
            return Err(failure(
                Kind::UnsupportedChemistry,
                Stage::Input,
                format!("unsupported bond order {:?} at {i}", bond.order()),
            )
            .with_bonds(vec![i]));
        }
    }
    if !topology.substance_groups.is_empty() {
        return Err(failure(
            Kind::UnsupportedChemistry,
            Stage::Input,
            "polymer/substance-group state is outside the finite MSReady input model",
        ));
    }
    Ok(())
}

pub fn validate_input(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &MsReadyParams,
) -> Result<ValidationClass, MsReadyFailure> {
    validate_input_inner(topology, properties, params).map_err(|mut error| {
        error.original = Some(crate::MsReadyOriginal {
            topology: topology.clone(),
            properties: properties.clone(),
        });
        error
    })
}

fn validate_input_inner(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &MsReadyParams,
) -> Result<ValidationClass, MsReadyFailure> {
    validate_syntax(topology, properties, params)?;
    if crate::exclusions::is_exact_raw_nitrate(topology)? {
        return Ok(ValidationClass::RawNitrateExclusionOnly);
    }
    if crate::disconnect::alkali_oxygen_candidate(topology)?.is_some() {
        crate::disconnect::validate_pending_repairs(topology, true, &[])?;
        return Ok(ValidationClass::RecognizedAlkaliRepair);
    }
    let metals = retained_metal_states(topology)?;
    selected_sanitize(topology)?;
    if metals.is_empty() {
        return Ok(ValidationClass::Ordinary);
    }
    Ok(ValidationClass::Ordinary)
}

pub fn retained_metal_states(
    topology: &TopologyBlock,
) -> Result<Vec<(usize, RetainedMetalKind)>, MsReadyFailure> {
    let inferred_hydrogens = cosmolkit_core::assign_valence_with_options_for_topology(
        topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| failure(Kind::OwnerError, Stage::Validation, error.to_string()))?
    .implicit_hydrogens;
    // Reuse the search owner's canonical metal vocabulary; profile policy,
    // not owner permissive valence, decides which of these states may remain.
    let metal_query = cosmolkit_search::convert_complex_name_to_query("M")
        .map_err(|error| failure(Kind::OwnerError, Stage::Validation, error.to_string()))?;
    let coordinates = cosmolkit_model::CoordinateBlock::default();
    let target = cosmolkit_search::SearchTarget::new(topology, &coordinates, &[], None, None);
    let context = cosmolkit_search::build_topology_query_match_context(topology)
        .map_err(|error| failure(Kind::OwnerError, Stage::Validation, error.to_string()))?;
    let mut accepted = Vec::new();
    for (i, atom) in topology.atoms.iter().enumerate() {
        if !matches!(
            atom.element(),
            Element::HG | Element::SN | Element::LI | Element::NA | Element::K
        ) {
            if cosmolkit_search::atom_matches_query_with_context(
                atom,
                &metal_query,
                &target,
                &context,
            )
            .map_err(|error| failure(Kind::OwnerError, Stage::Validation, error.to_string()))?
            {
                return Err(profile::unsupported_at(
                    i,
                    format!(
                        "unmodeled retained metal {} at atom {i}",
                        atom.element().symbol()
                    ),
                ));
            }
            continue;
        }
        let neighbors = topology.adjacency.neighbors_of(i);
        let isotope_h = atom.tracked_isotopic_hydrogens().len();
        let heavy_neighbors: Vec<_> = neighbors
            .iter()
            .filter(|n| topology.atoms[n.atom_index].atomic_number() != 1)
            .collect();
        let inferred = inferred_hydrogens
            .get(i)
            .copied()
            .ok_or_else(|| {
                failure(
                    Kind::InvariantViolation,
                    Stage::Validation,
                    format!("missing implicit H assignment at atom {i}"),
                )
            })?
            .max(0) as usize;
        let h_count = atom.explicit_hydrogens() as usize + neighbors.len() - heavy_neighbors.len()
            + isotope_h
            + inferred;
        if atom.is_aromatic() || atom.radical_electrons() != 0 || h_count != 0 {
            return Err(profile::unsupported_at(
                i,
                format!(
                    "metal {} at atom {i} has unsupported aromatic/radical/H state",
                    atom.element().symbol()
                ),
            ));
        }
        for neighbor in &heavy_neighbors {
            let bond = &topology.bonds[neighbor.bond.index()];
            if !matches!(bond.order(), BondOrder::Single)
                || bond.is_aromatic()
                || bond.query().is_some()
            {
                return Err(profile::unsupported_at(
                    i,
                    format!(
                        "metal {} at atom {i} has unsupported bond",
                        atom.element().symbol()
                    ),
                ));
            }
        }
        let state = match atom.element() {
            Element::HG
                if atom.formal_charge() == 1
                    && neighbors.len() == 1
                    && heavy_neighbors.len() == 1
                    && carbon_only(topology, &heavy_neighbors) =>
            {
                RetainedMetalKind::HgCationMonovalent
            }
            Element::HG
                if atom.formal_charge() == 0
                    && neighbors.len() == 2
                    && heavy_neighbors.len() == 2
                    && carbon_only(topology, &heavy_neighbors) =>
            {
                RetainedMetalKind::HgNeutralDivalent
            }
            Element::SN
                if atom.formal_charge() == 0
                    && neighbors.len() == 4
                    && heavy_neighbors.len() == 4
                    && carbon_only(topology, &heavy_neighbors) =>
            {
                RetainedMetalKind::SnNeutralTetravalent
            }
            Element::LI | Element::NA | Element::K
                if atom.formal_charge() == 1 && neighbors.is_empty() =>
            {
                RetainedMetalKind::AlkaliSpectator
            }
            _ => {
                return Err(profile::unsupported_at(
                    i,
                    format!(
                        "unmodeled retained metal state {} charge {} degree {} at atom {i}",
                        atom.element().symbol(),
                        atom.formal_charge(),
                        neighbors.len()
                    ),
                ));
            }
        };
        if !matches!(state, RetainedMetalKind::AlkaliSpectator) {
            for neighbor in &heavy_neighbors {
                let carbon = &topology.atoms[neighbor.atom_index];
                if carbon.formal_charge() != 0
                    || carbon.is_aromatic()
                    || carbon.radical_electrons() != 0
                {
                    return Err(profile::unsupported_at(
                        i,
                        format!(
                            "metal atom {i} is attached to a non-neutral/nonaromatic carbon state"
                        ),
                    ));
                }
            }
        }
        accepted.push((i, state));
    }
    Ok(accepted)
}

fn carbon_only(topology: &TopologyBlock, neighbors: &[&cosmolkit_model::NeighborRef]) -> bool {
    neighbors
        .iter()
        .all(|n| topology.atoms[n.atom_index].element() == Element::C)
}

fn selected_sanitize(
    topology: &TopologyBlock,
) -> Result<cosmolkit_core::SanitizeAssignment, MsReadyFailure> {
    let selected = SanitizeOperations::PROPERTIES
        | SanitizeOperations::SYMM_RINGS
        | SanitizeOperations::KEKULIZE
        | SanitizeOperations::FIND_RADICALS
        | SanitizeOperations::SET_AROMATICITY
        | SanitizeOperations::SET_CONJUGATION
        | SanitizeOperations::SET_HYBRIDIZATION
        | SanitizeOperations::CLEANUP_CHIRALITY
        | SanitizeOperations::ADJUST_HS;
    sanitize_topology(
        topology,
        &SanitizeParams {
            operations: selected,
        },
    )
    .map_err(|error| {
        let kind = match error {
            cosmolkit_core::SanitizeError::Properties { .. } => Kind::InvalidValence,
            _ => Kind::OwnerError,
        };
        failure(kind, Stage::Validation, error.to_string())
    })
}

pub(crate) fn strict_owner_valence(
    topology: &TopologyBlock,
) -> Result<cosmolkit_core::ValenceAssignment, MsReadyFailure> {
    let assignment = sanitize_topology(
        topology,
        &SanitizeParams {
            operations: SanitizeOperations::PROPERTIES,
        },
    )
    .map_err(|error| {
        let kind = match error {
            cosmolkit_core::SanitizeError::Properties { .. } => Kind::InvalidValence,
            _ => Kind::OwnerError,
        };
        let mut failure = failure(kind, Stage::Validation, error.to_string());
        failure.evidence.location.owner_category =
            Some(format!("cosmolkit_core::SanitizeError::{error:?}"));
        let mut source = std::error::Error::source(&error);
        while let Some(underlying) = source {
            if let Some(cosmolkit_core::ValenceError::InvalidValence { atom, .. }) =
                underlying.downcast_ref::<cosmolkit_core::ValenceError>()
            {
                failure.evidence.location.atom_indices = vec![atom.index()];
                break;
            }
            source = underlying.source();
        }
        failure
    })?;
    assignment.final_valence.ok_or_else(|| {
        failure(
            Kind::InvariantViolation,
            Stage::Validation,
            "strict owner validation returned no final valence assignment",
        )
    })
}

#[cfg(test)]
mod validation_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, BondSpec, TopologyBlock};
    fn graph(specs: Vec<AtomSpec>, edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = specs
            .into_iter()
            .enumerate()
            .map(|(i, spec)| cosmolkit_model::Atom::from_spec(AtomId::new(i), spec))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                cosmolkit_model::Bond::from_spec(
                    cosmolkit_model::BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn a(e: Element, q: i8) -> AtomSpec {
        AtomSpec::new(e).with_formal_charge(q)
    }
    #[test]
    fn retained_other_metals_are_unsupported_after_exact_exclusion_comparison() {
        for smiles in ["C[Fe]C", "C[Cu]C", "C[Zn]C"] {
            let graph = cosmolkit_smiles::parse_smiles(smiles, &Default::default())
                .unwrap()
                .topology;
            assert!(crate::match_exact_exclusions(&graph).unwrap().is_empty());
            assert_eq!(
                validate_input(&graph, &Default::default(), &Default::default())
                    .unwrap_err()
                    .kind,
                Kind::UnsupportedChemistry,
                "{smiles}"
            );
        }
        for smiles in ["[Zn]", "[Ag]", "[Na]"] {
            let graph = cosmolkit_smiles::parse_smiles(smiles, &Default::default())
                .unwrap()
                .topology;
            assert!(
                !crate::match_exact_exclusions(&graph).unwrap().is_empty(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn raw_nitrate_exemption_uses_exact_all_atom_and_hydrogen_state() {
        let raw = cosmolkit_smiles::parse_smiles("ON(=O)O", &Default::default())
            .unwrap()
            .topology;
        let n = raw
            .atoms
            .iter()
            .position(|atom| atom.element() == Element::N)
            .unwrap();
        let single_o = raw
            .bonds
            .iter()
            .find(|bond| bond.order() == BondOrder::Single)
            .unwrap()
            .begin()
            .index();
        let double_o = raw
            .atoms
            .iter()
            .enumerate()
            .find(|(row, atom)| {
                atom.element() == Element::O
                    && raw.bonds[raw.adjacency.neighbors_of(*row)[0].bond.index()].order()
                        == BondOrder::Double
            })
            .unwrap()
            .0;
        for case in 0..8 {
            let mut altered = raw.clone();
            match case {
                0 => altered.atoms[n].set_isotope(Some(15)),
                1 => altered.atoms[n].set_explicit_hydrogens(1),
                2 => altered.atoms[n].set_radical_electrons(1),
                3 => altered.atoms[n].set_aromatic(true),
                4 => altered.atoms[n].set_tracked_isotopic_hydrogens(vec![2]),
                5 => altered.atoms[double_o].set_explicit_hydrogens(1),
                6 => altered.atoms[single_o].set_no_implicit(true),
                7 => altered.atoms[single_o].set_tracked_isotopic_hydrogens(vec![3]),
                _ => unreachable!(),
            }
            assert!(
                crate::match_exact_exclusions(&altered).unwrap().is_empty(),
                "case {case}"
            );
            let result = validate_input(&altered, &Default::default(), &Default::default());
            // Aromatic N changes owner valence rounding and instead fails
            // selected kekulization (SanitizeError::Kekulize -> OwnerError).
            let expected = if case == 3 {
                Kind::OwnerError
            } else {
                Kind::InvalidValence
            };
            assert_eq!(result.unwrap_err().kind, expected, "case {case}");
        }
    }

    #[test]
    fn syntax_gate_rejects_empty_wildcard_and_zero_limits() {
        let params = MsReadyParams::default();
        let empty_error = validate_syntax(
            &TopologyBlock::default(),
            &MoleculeProperties::default(),
            &params,
        )
        .unwrap_err();
        assert_eq!(empty_error.kind, Kind::InvalidInput);
        assert_eq!(
            empty_error.original.unwrap().topology,
            TopologyBlock::default()
        );
        let wildcard = graph(vec![a(Element::from_atomic_number(0).unwrap(), 0)], &[]);
        assert_eq!(
            validate_syntax(&wildcard, &MoleculeProperties::default(), &params)
                .unwrap_err()
                .kind,
            Kind::UnsupportedChemistry
        );
        let mut limits = params;
        limits.max_edits = 0;
        assert_eq!(
            validate_syntax(
                &graph(vec![a(Element::C, 0)], &[]),
                &MoleculeProperties::default(),
                &limits
            )
            .unwrap_err()
            .kind,
            Kind::InvalidInput
        );
    }

    #[test]
    fn accepts_exact_hg_sn_tuples_and_rejects_neighbors_outside_domain() {
        let cases = [
            (
                graph(vec![a(Element::C, 0), a(Element::HG, 1)], &[(0, 1)]),
                RetainedMetalKind::HgCationMonovalent,
            ),
            (
                graph(
                    vec![a(Element::C, 0), a(Element::HG, 0), a(Element::C, 0)],
                    &[(0, 1), (1, 2)],
                ),
                RetainedMetalKind::HgNeutralDivalent,
            ),
            (
                graph(
                    vec![
                        a(Element::C, 0),
                        a(Element::SN, 0),
                        a(Element::C, 0),
                        a(Element::C, 0),
                        a(Element::C, 0),
                    ],
                    &[(0, 1), (1, 2), (1, 3), (1, 4)],
                ),
                RetainedMetalKind::SnNeutralTetravalent,
            ),
        ];
        for (t, want) in cases {
            assert_eq!(retained_metal_states(&t).unwrap()[0].1, want);
        }
        let bad_cases = [
            graph(vec![a(Element::C, 0), a(Element::HG, 0)], &[(0, 1)]),
            graph(
                vec![a(Element::C, 0), a(Element::HG, 1), a(Element::C, 0)],
                &[(0, 1), (1, 2)],
            ),
            graph(vec![a(Element::C, 0), a(Element::HG, 2)], &[(0, 1)]),
            graph(
                vec![a(Element::C, 0), a(Element::HG, 0), a(Element::O, 0)],
                &[(0, 1), (1, 2)],
            ),
            graph(
                vec![
                    a(Element::C, 0),
                    a(Element::SN, 0),
                    a(Element::C, 0),
                    a(Element::C, 0),
                ],
                &[(0, 1), (1, 2), (1, 3)],
            ),
            graph(
                vec![
                    a(Element::C, 0),
                    a(Element::SN, 1),
                    a(Element::C, 0),
                    a(Element::C, 0),
                    a(Element::C, 0),
                ],
                &[(0, 1), (1, 2), (1, 3), (1, 4)],
            ),
            graph(
                vec![
                    a(Element::C, 0),
                    a(Element::SN, 0),
                    a(Element::C, 0),
                    a(Element::C, 0),
                    a(Element::O, 0),
                ],
                &[(0, 1), (1, 2), (1, 3), (1, 4)],
            ),
            graph(
                vec![a(Element::H, 0), a(Element::HG, 0), a(Element::C, 0)],
                &[(0, 1), (1, 2)],
            ),
            graph(
                vec![
                    a(Element::H, 0),
                    a(Element::SN, 0),
                    a(Element::C, 0),
                    a(Element::C, 0),
                    a(Element::C, 0),
                ],
                &[(0, 1), (1, 2), (1, 3), (1, 4)],
            ),
        ];
        for bad in bad_cases {
            assert!(
                retained_metal_states(&bad).is_err(),
                "unexpectedly accepted {bad:?}"
            );
        }
    }
    #[test]
    fn raw_neutral_nitrate_is_exclusion_only_and_its_charge_variants_are_not_repairs() {
        let nitrate = graph(
            vec![
                a(Element::N, 0),
                a(Element::O, 0).with_explicit_hydrogens(1),
                a(Element::O, 0).with_explicit_hydrogens(1),
                a(Element::O, 0),
            ],
            &[(0, 1), (0, 2), (0, 3)],
        );
        let mut nitrate = nitrate;
        nitrate.bonds[2].set_order(BondOrder::Double);
        assert!(crate::exclusions::is_exact_raw_nitrate(&nitrate).unwrap());
        assert_eq!(
            validate_input(
                &nitrate,
                &MoleculeProperties::default(),
                &MsReadyParams::default()
            ),
            Ok(ValidationClass::RawNitrateExclusionOnly)
        );
        let mut charged = nitrate.clone();
        charged.atoms[0].set_formal_charge(1);
        assert!(!crate::exclusions::is_exact_raw_nitrate(&charged).unwrap());
    }

    #[test]
    fn validation_keeps_supported_graph_bonds_and_reports_unmodeled_metal_tuples() {
        let accepted = graph(vec![a(Element::C, 0), a(Element::HG, 1)], &[(0, 1)]);
        let before = accepted.clone();
        let result = validate_input(
            &accepted,
            &MoleculeProperties::default(),
            &MsReadyParams::default(),
        );
        assert!(
            result.is_ok(),
            "accepted Hg tuple failed validation: {result:?}"
        );
        assert_eq!(accepted, before);
        assert_eq!(accepted.bonds[0].order(), BondOrder::Single);
        let unmodeled = graph(vec![a(Element::C, 0), a(Element::HG, 2)], &[(0, 1)]);
        let error = retained_metal_states(&unmodeled).unwrap_err();
        assert_eq!(error.kind, Kind::UnsupportedChemistry);
        assert!(error.context.contains("degree 1"));
    }

    #[test]
    fn ignores_only_the_exact_bonded_charged_alkali_oxygen_precheck() {
        let t = graph(
            vec![
                a(Element::NA, 1),
                a(Element::O, -1).with_explicit_hydrogens(1),
            ],
            &[(0, 1)],
        );
        assert_eq!(
            crate::disconnect::alkali_oxygen_candidate(&t)
                .unwrap()
                .map(|candidate| candidate.bond_index),
            Some(0)
        );
        let t = graph(vec![a(Element::NA, 0), a(Element::O, -1)], &[(0, 1)]);
        assert!(
            crate::disconnect::alkali_oxygen_candidate(&t)
                .unwrap()
                .is_none()
        );

        let neutral = graph(vec![a(Element::NA, 0), a(Element::O, 0)], &[(0, 1)]);
        assert_eq!(
            validate_input(
                &neutral,
                &MoleculeProperties::default(),
                &MsReadyParams::default(),
            )
            .unwrap(),
            ValidationClass::RecognizedAlkaliRepair
        );
        let rejected = graph(vec![a(Element::NA, 1), a(Element::O, 0)], &[(0, 1)]);
        assert!(
            validate_input(
                &rejected,
                &MoleculeProperties::default(),
                &MsReadyParams::default(),
            )
            .is_err()
        );
    }
}
