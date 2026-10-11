use cosmolkit_core::get_molecule_fragments;
use cosmolkit_model::{
    AtomId, BondId, CoordinateBlock, MoleculeProperties, TopologyBlock, TopologyMapping,
};

use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};

#[derive(Clone, Debug, PartialEq)]
pub struct MappedComponent {
    pub topology: TopologyBlock,
    /// Original input atom rows in component topology order.
    pub input_atom_indices: Vec<usize>,
    /// Original-input-to-component atom and bond mapping.
    pub topology_mapping: TopologyMapping,
}

/// Materialize source-ordered graph components through core's shared copier.
pub fn split_mapped_components(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<Vec<MappedComponent>, MsReadyFailure> {
    let source_atom_rows = (0..topology.atoms.len()).collect::<Vec<_>>();
    split_mapped_components_from_source_rows(topology, coordinates, properties, &source_atom_rows)
}

/// Split an edited topology and compose each component map through the
/// original-input-to-current-topology mapping.
pub fn split_mapped_components_from_source_mapping(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    original_to_current: &TopologyMapping,
) -> Result<Vec<MappedComponent>, MsReadyFailure> {
    original_to_current
        .validate_for_counts(
            original_to_current.atoms.old_to_new.len(),
            topology.atoms.len(),
            original_to_current.bonds.old_to_new.len(),
            topology.bonds.len(),
        )
        .map_err(|error| MsReadyFailure {
            evidence: Default::default(),
            kind: Kind::InvariantViolation,
            stage: Stage::Validation,
            context: format!("invalid original-to-current topology mapping: {error}"),
            original: None,
        })?;
    let source_atom_rows = original_to_current
        .atoms
        .new_to_old
        .iter()
        .map(|id| {
            id.map(|id| id.index()).ok_or_else(|| MsReadyFailure {
                evidence: Default::default(),
                kind: Kind::InvariantViolation,
                stage: Stage::Validation,
                context: "current atom has no original source row".to_owned(),
                original: None,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut components = split_mapped_components_from_source_rows(
        topology,
        coordinates,
        properties,
        &source_atom_rows,
    )?;
    for component in &mut components {
        component.topology_mapping = compose_component_mapping(
            original_to_current,
            &component.topology_mapping,
            topology,
            &component.topology,
        )?;
    }
    Ok(components)
}

fn compose_component_mapping(
    original_to_current: &TopologyMapping,
    current_to_component: &TopologyMapping,
    current: &TopologyBlock,
    component: &TopologyBlock,
) -> Result<TopologyMapping, MsReadyFailure> {
    current_to_component
        .validate_for_counts(
            current.atoms.len(),
            component.atoms.len(),
            current.bonds.len(),
            component.bonds.len(),
        )
        .map_err(|error| MsReadyFailure {
            evidence: Default::default(),
            kind: Kind::InvariantViolation,
            stage: Stage::Validation,
            context: format!("invalid current-to-component mapping: {error}"),
            original: None,
        })?;
    let mut atoms_old_to_new = vec![None; original_to_current.atoms.old_to_new.len()];
    let mut atoms_new_to_old = Vec::with_capacity(component.atoms.len());
    for current_id in current_to_component.atoms.new_to_old.iter().flatten() {
        let original_id =
            original_to_current.atoms.new_to_old[current_id.index()].ok_or_else(|| {
                MsReadyFailure {
                    evidence: Default::default(),
                    kind: Kind::InvariantViolation,
                    stage: Stage::Validation,
                    context: "component atom has no original source row".to_owned(),
                    original: None,
                }
            })?;
        atoms_old_to_new[original_id.index()] = Some(AtomId::new(atoms_new_to_old.len()));
        atoms_new_to_old.push(Some(original_id));
    }
    let mut bonds_old_to_new = vec![None; original_to_current.bonds.old_to_new.len()];
    let mut bonds_new_to_old = Vec::with_capacity(component.bonds.len());
    for current_id in current_to_component.bonds.new_to_old.iter().flatten() {
        let original_id =
            original_to_current.bonds.new_to_old[current_id.index()].ok_or_else(|| {
                MsReadyFailure {
                    evidence: Default::default(),
                    kind: Kind::InvariantViolation,
                    stage: Stage::Validation,
                    context: "component bond has no original source row".to_owned(),
                    original: None,
                }
            })?;
        bonds_old_to_new[original_id.index()] = Some(BondId::new(bonds_new_to_old.len()));
        bonds_new_to_old.push(Some(original_id));
    }
    let mapping = TopologyMapping {
        atoms: cosmolkit_model::AtomMapping {
            old_to_new: atoms_old_to_new,
            new_to_old: atoms_new_to_old,
        },
        bonds: cosmolkit_model::BondMapping {
            old_to_new: bonds_old_to_new,
            new_to_old: bonds_new_to_old,
        },
    };
    mapping
        .validate_for_counts(
            original_to_current.atoms.old_to_new.len(),
            component.atoms.len(),
            original_to_current.bonds.old_to_new.len(),
            component.bonds.len(),
        )
        .map_err(|error| MsReadyFailure {
            evidence: Default::default(),
            kind: Kind::InvariantViolation,
            stage: Stage::Validation,
            context: format!("invalid composed component mapping: {error}"),
            original: None,
        })?;
    Ok(mapping)
}

/// Split an edited graph and compose core component rows with a caller's
/// current-topology-to-original-input atom map.
pub fn split_mapped_components_from_source_rows(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    source_atom_rows: &[usize],
) -> Result<Vec<MappedComponent>, MsReadyFailure> {
    if source_atom_rows.len() != topology.atoms.len() {
        return Err(MsReadyFailure {
            evidence: Default::default(),
            kind: Kind::InvariantViolation,
            stage: Stage::Validation,
            context: format!(
                "source atom map has {} rows for {} topology atoms",
                source_atom_rows.len(),
                topology.atoms.len()
            ),
            original: None,
        });
    }
    let fragments = get_molecule_fragments(topology, coordinates, properties, false, false)
        .map_err(|e| MsReadyFailure {
            evidence: Default::default(),
            kind: Kind::OwnerError,
            stage: Stage::Validation,
            context: e.to_string(),
            original: None,
        })?;
    fragments
        .into_iter()
        .map(|fragment| {
            let input_atom_indices = fragment
                .component_atoms()
                .iter()
                .map(|id| source_atom_rows[id.index()])
                .collect();
            let topology_mapping = fragment.topology_mapping().clone();
            let (topology, _, _) = fragment.into_parts();
            Ok(MappedComponent {
                topology,
                input_atom_indices,
                topology_mapping,
            })
        })
        .collect()
}

#[cfg(test)]
mod component_mapping_tests {
    use super::*;
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

    #[test]
    fn core_fragment_mapping_covers_each_original_atom_once_in_source_order() {
        let source = parse_smiles("CC(=O)[O-].[Na+]", &SmilesParseParams::default()).unwrap();
        let parts =
            split_mapped_components(&source.topology, &source.coordinates, &source.properties)
                .unwrap();
        assert_eq!(parts.len(), 2);
        let mut mapped: Vec<_> = parts
            .iter()
            .flat_map(|part| part.input_atom_indices.iter().copied())
            .collect();
        mapped.sort_unstable();
        assert_eq!(mapped, (0..source.topology.atoms.len()).collect::<Vec<_>>());
        assert_eq!(parts[0].topology.atoms.len(), 4);
        assert_eq!(parts[1].topology.atoms.len(), 1);
        for part in parts {
            assert_eq!(
                part.topology_mapping.atoms.new_to_old.len(),
                part.topology.atoms.len()
            );
            assert_eq!(
                part.topology_mapping.bonds.new_to_old.len(),
                part.topology.bonds.len()
            );
        }
    }

    #[test]
    fn disconnection_and_core_component_maps_compose_to_original_atom_rows() {
        let source = parse_smiles("CC(=O)O[Na]", &SmilesParseParams::default()).unwrap();
        let run = crate::disconnect_all_alkali_oxygen(&source.topology).unwrap();
        let parts = split_mapped_components_from_source_mapping(
            &run.topology,
            &source.coordinates,
            &source.properties,
            &run.topology_mapping,
        )
        .unwrap();
        assert_eq!(parts.len(), 2);
        let mut original_rows: Vec<_> = parts
            .iter()
            .flat_map(|part| part.input_atom_indices.iter().copied())
            .collect();
        original_rows.sort_unstable();
        assert_eq!(
            original_rows,
            (0..source.topology.atoms.len()).collect::<Vec<_>>()
        );
        for part in &parts {
            for (local_row, old_row) in part.input_atom_indices.iter().copied().enumerate() {
                assert_eq!(
                    part.topology_mapping.atoms.new_to_old[local_row]
                        .unwrap()
                        .index(),
                    old_row
                );
            }
        }
        let mapped_bond_rows: std::collections::BTreeSet<_> = parts
            .iter()
            .flat_map(|part| {
                part.topology_mapping
                    .bonds
                    .new_to_old
                    .iter()
                    .flatten()
                    .map(|id| id.index())
            })
            .collect();
        assert_eq!(mapped_bond_rows.len(), source.topology.bonds.len() - 1);
        let cut_bond = run.edits[0].bond_index;
        assert!(run.topology_mapping.bonds.old_to_new[cut_bond].is_none());
        assert!(!mapped_bond_rows.contains(&cut_bond));
    }
}
