use std::collections::BTreeSet;

use cosmolkit_core::{RemoveHsParams, ValenceAssignment, remove_hydrogens_with_params};
use cosmolkit_model::{
    AtomId, BondId, CoordinateBlock, MoleculeProperties, TopologyBatchEdit, TopologyBlock,
    TopologyMapping,
};
use cosmolkit_types::BondOrder;

use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ProtonCarrier {
    Implicit,
    Bracket,
    ExplicitAtom {
        atom: usize,
        bond: usize,
        isotope: Option<u16>,
    },
    TrackedIsotope {
        index: usize,
        isotope: u16,
    },
}

impl ProtonCarrier {
    pub(crate) const fn is_labelled(self) -> bool {
        matches!(
            self,
            Self::ExplicitAtom {
                isotope: Some(_),
                ..
            } | Self::TrackedIsotope { .. }
        )
    }

    pub(crate) fn preference(self) -> (u8, usize, u16) {
        match self {
            Self::Bracket => (0, 0, 0),
            Self::ExplicitAtom {
                atom,
                isotope: None,
                ..
            } => (1, atom, 0),
            Self::Implicit => (2, 0, 0),
            Self::ExplicitAtom {
                atom,
                isotope: Some(isotope),
                ..
            } => (3, atom, isotope),
            Self::TrackedIsotope { index, isotope } => (4, index, isotope),
        }
    }
}

pub(crate) use crate::types::{
    MsReadyHydrogenDisposition as HydrogenDisposition,
    MsReadyHydrogenDispositionKind as HydrogenDispositionKind,
};

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct HydrogenFoldRun {
    pub topology: TopologyBlock,
    pub topology_mapping: TopologyMapping,
    pub dispositions: Vec<HydrogenDisposition>,
}

pub(crate) fn available_carriers(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    donor: usize,
) -> Result<Vec<ProtonCarrier>, MsReadyFailure> {
    // BEGIN MSREADY REFERENCE §2
    // MSReady❗❗: Keep explicit hydrogen atoms distinct from bracket hydrogen counts and inferred hydrogens.
    // END MSREADY REFERENCE
    let atom = topology.atoms.get(donor).ok_or_else(|| {
        failure(
            Kind::InvariantViolation,
            format!("proton donor atom row {donor} is out of range"),
        )
    })?;
    let mut carriers = Vec::new();
    if atom.explicit_hydrogens() != 0 {
        carriers.push(ProtonCarrier::Bracket);
    }
    for neighbor in topology.adjacency.neighbors_of(donor) {
        let hydrogen = &topology.atoms[neighbor.atom_index];
        let bond = &topology.bonds[neighbor.bond.index()];
        if hydrogen.atomic_number() != 1
            || hydrogen.formal_charge() != 0
            || hydrogen.radical_electrons() != 0
            || topology.adjacency.neighbors_of(neighbor.atom_index).len() != 1
            || bond.order() != BondOrder::Single
            || bond.is_aromatic()
            || bond.query().is_some()
        {
            continue;
        }
        carriers.push(ProtonCarrier::ExplicitAtom {
            atom: neighbor.atom_index,
            bond: neighbor.bond.index(),
            isotope: hydrogen.isotope(),
        });
    }
    for (index, isotope) in atom
        .tracked_isotopic_hydrogens()
        .iter()
        .copied()
        .enumerate()
    {
        carriers.push(ProtonCarrier::TrackedIsotope { index, isotope });
    }
    let implicit = valence
        .implicit_hydrogens
        .get(donor)
        .copied()
        .ok_or_else(|| {
            failure(
                Kind::OwnerError,
                format!("missing inferred H count at atom row {donor}"),
            )
        })?;
    if implicit < 0 {
        return Err(failure(
            Kind::OwnerError,
            format!("negative inferred H count at atom row {donor}"),
        ));
    }
    let implicit = implicit as u32;
    let ordinary_implicit = implicit.saturating_sub(atom.tracked_isotopic_hydrogens().len() as u32);
    if !atom.no_implicit() && ordinary_implicit != 0 {
        carriers.push(ProtonCarrier::Implicit);
    }
    carriers.sort_by_key(|carrier| carrier.preference());
    carriers.dedup();
    Ok(carriers)
}

pub(crate) fn transfer_carrier(
    edit: &mut TopologyBatchEdit,
    topology: &TopologyBlock,
    donor: usize,
    acceptor: usize,
    carrier: ProtonCarrier,
) -> Result<(), MsReadyFailure> {
    match carrier {
        ProtonCarrier::Implicit => {}
        ProtonCarrier::Bracket => {
            let donor_atom = edit
                .atom_mut(AtomId::new(donor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            donor_atom.set_explicit_hydrogens(
                donor_atom
                    .explicit_hydrogens()
                    .checked_sub(1)
                    .ok_or_else(|| {
                        failure(Kind::InvariantViolation, "bracket-H donor count underflow")
                    })?,
            );
            let acceptor_atom = edit
                .atom_mut(AtomId::new(acceptor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            acceptor_atom.set_explicit_hydrogens(
                acceptor_atom
                    .explicit_hydrogens()
                    .checked_add(1)
                    .ok_or_else(|| {
                        failure(
                            Kind::InvariantViolation,
                            "bracket-H acceptor count overflow",
                        )
                    })?,
            );
        }
        ProtonCarrier::ExplicitAtom {
            atom,
            bond,
            isotope,
        } => {
            let source_h = topology.atoms.get(atom).ok_or_else(|| {
                failure(
                    Kind::InvariantViolation,
                    "selected explicit H atom is out of range",
                )
            })?;
            if source_h.atomic_number() != 1 || source_h.isotope() != isotope {
                return Err(failure(
                    Kind::InvariantViolation,
                    "selected explicit H carrier no longer matches source state",
                ));
            }
            let bond_id = BondId::new(bond);
            let bond_mut = edit
                .bond_mut(bond_id)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let hydrogen_id = AtomId::new(atom);
            let donor_id = AtomId::new(donor);
            let acceptor_id = AtomId::new(acceptor);
            if bond_mut.begin() == donor_id && bond_mut.end() == hydrogen_id {
                bond_mut.set_endpoints(acceptor_id, hydrogen_id);
            } else if bond_mut.end() == donor_id && bond_mut.begin() == hydrogen_id {
                bond_mut.set_endpoints(hydrogen_id, acceptor_id);
            } else {
                return Err(failure(
                    Kind::InvariantViolation,
                    "explicit H bond does not join selected donor and H atom",
                ));
            }
        }
        ProtonCarrier::TrackedIsotope { index, isotope } => {
            let donor_atom = edit
                .atom_mut(AtomId::new(donor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let mut donor_isotopes = donor_atom.tracked_isotopic_hydrogens().to_vec();
            if donor_isotopes.get(index).copied() != Some(isotope) {
                return Err(failure(
                    Kind::InvariantViolation,
                    "selected tracked isotope H slot no longer matches source",
                ));
            }
            donor_isotopes.remove(index);
            donor_atom.set_tracked_isotopic_hydrogens(donor_isotopes);
            let acceptor_atom = edit
                .atom_mut(AtomId::new(acceptor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let mut acceptor_isotopes = acceptor_atom.tracked_isotopic_hydrogens().to_vec();
            acceptor_isotopes.push(isotope);
            acceptor_atom.set_tracked_isotopic_hydrogens(acceptor_isotopes);
        }
    }
    Ok(())
}

pub(crate) fn remove_reservoir_carrier(
    edit: &mut TopologyBatchEdit,
    topology: &TopologyBlock,
    donor: usize,
    carrier: ProtonCarrier,
) -> Result<Option<usize>, MsReadyFailure> {
    match carrier {
        ProtonCarrier::Implicit => Ok(None),
        ProtonCarrier::Bracket => {
            let atom = edit
                .atom_mut(AtomId::new(donor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let remaining = atom.explicit_hydrogens().checked_sub(1).ok_or_else(|| {
                failure(
                    Kind::InvariantViolation,
                    "bracket-H reservoir donor count underflow",
                )
            })?;
            atom.set_explicit_hydrogens(remaining);
            Ok(None)
        }
        ProtonCarrier::ExplicitAtom {
            atom,
            bond,
            isotope: None,
        } => {
            let bond = topology.bonds.get(bond).ok_or_else(|| {
                failure(
                    Kind::InvariantViolation,
                    "selected explicit-H bond is out of range",
                )
            })?;
            let donor_id = AtomId::new(donor);
            let hydrogen_id = AtomId::new(atom);
            if !((bond.begin() == donor_id && bond.end() == hydrogen_id)
                || (bond.end() == donor_id && bond.begin() == hydrogen_id))
            {
                return Err(failure(
                    Kind::InvariantViolation,
                    "selected explicit-H bond does not join donor and H atom",
                ));
            }
            edit.remove_atom(hydrogen_id)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            Ok(Some(atom))
        }
        ProtonCarrier::ExplicitAtom {
            isotope: Some(_), ..
        }
        | ProtonCarrier::TrackedIsotope { .. } => Err(failure(
            Kind::InvariantViolation,
            "labelled D/T cannot be exchanged with the ordinary proton reservoir",
        )),
    }
}

pub(crate) fn add_reservoir_proton(
    edit: &mut TopologyBatchEdit,
    acceptor: usize,
) -> Result<(), MsReadyFailure> {
    let atom = edit
        .atom_mut(AtomId::new(acceptor))
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    if atom.no_implicit() {
        atom.set_explicit_hydrogens(atom.explicit_hydrogens().checked_add(1).ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                "bracket-H reservoir acceptor count overflow",
            )
        })?);
    }
    Ok(())
}

/// Materialize tracked labelled H only on a transient identifier/comparison graph.
/// Tracked H replaces inferred H (or augments a no-implicit carrier), never
/// bracket H. Adding its bond lets the owner recompute the remaining ordinary H.
pub(crate) fn materialize_tracked_isotopic_hydrogens(
    topology: &TopologyBlock,
) -> Result<TopologyBlock, MsReadyFailure> {
    let mut edit = topology
        .begin_batch_edit()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    for (row, atom) in topology.atoms.iter().enumerate() {
        if atom.tracked_isotopic_hydrogens().is_empty() {
            continue;
        }
        edit.atom_mut(AtomId::new(row))
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
            .set_tracked_isotopic_hydrogens(Vec::new());
        for &isotope in atom.tracked_isotopic_hydrogens() {
            let hydrogen = edit.add_atom(
                cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::H)
                    .with_isotope(isotope)
                    .with_no_implicit(true),
            );
            edit.add_bond(cosmolkit_model::BondSpec::new(
                AtomId::new(row),
                hydrogen,
                BondOrder::Single,
            ))
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        }
    }
    let (graph, _) = edit
        .finish()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    Ok(graph)
}

/// Fold only source-eligible ordinary explicit H atoms through core RemoveHs.
pub(crate) fn fold_ordinary_explicit_hydrogens(
    original: &TopologyBlock,
    topology: &TopologyBlock,
    original_to_current: &TopologyMapping,
) -> Result<HydrogenFoldRun, MsReadyFailure> {
    // BEGIN MSREADY REFERENCE §9
    // MSReady❗❗: Fold ordinary explicit H atoms when useful for pattern matching, then unfold or recompute them as the selected rules require. Bracket H counts, aromatic [nH], hydrides and isotope-labeled H need separate handling. Final ordinary SMILES can use implicit H even if intermediate structures contain explicit H atoms.
    // END MSREADY REFERENCE
    original_to_current
        .validate_for_counts(
            original.atoms.len(),
            topology.atoms.len(),
            original.bonds.len(),
            topology.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    let mut current_hydrogens = Vec::new();
    for (input_row, input_atom) in original.atoms.iter().enumerate() {
        if input_atom.atomic_number() != 1 {
            continue;
        }
        let Some(current_id) = original_to_current.atoms.old_to_new[input_row] else {
            continue;
        };
        let current_row = current_id.index();
        let current_atom = &topology.atoms[current_row];
        let parent = topology
            .adjacency
            .neighbors_of(current_row)
            .iter()
            .find(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() != 1)
            .map(|neighbor| neighbor.atom_index);
        current_hydrogens.push((input_row, current_row, parent, current_atom.isotope()));
    }

    let mut params = RemoveHsParams::default();
    params.remove_isotopes = false;
    params.remove_and_track_isotopes = false;
    params.remove_mapped = false;
    params.sanitize = false;
    let removed = remove_hydrogens_with_params(
        topology.clone(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &params,
    )
    .map_err(|error| {
        failure(
            Kind::OwnerError,
            format!("core explicit-H folding failed: {error}"),
        )
    })?;
    let mapping = crate::rewrite::compose_mappings(
        original_to_current,
        &removed.mapping,
        original,
        &removed.topology,
    )?;
    let valence = cosmolkit_core::assign_valence_for_topology(
        &removed.topology,
        cosmolkit_core::ValenceModel::RdkitLike,
    )
    .map_err(|error| failure(Kind::InvalidValence, error.to_string()))?;
    let _ = valence;
    let mut dispositions = Vec::with_capacity(current_hydrogens.len());
    let mut removed_current_rows = BTreeSet::new();
    for (input_row, current_row, parent_row, isotope) in current_hydrogens {
        if removed.mapping.atoms.old_to_new[current_row].is_none() {
            removed_current_rows.insert(current_row);
            let parent_row = parent_row.ok_or_else(|| {
                failure(
                    Kind::InvariantViolation,
                    format!("folded H input row {input_row} has no heavy-atom parent"),
                )
            })?;
            let input_parent_row = original_to_current.atoms.new_to_old[parent_row]
                .ok_or_else(|| {
                    failure(Kind::InvariantViolation, "folded H parent has no input map")
                })?
                .index();
            dispositions.push(HydrogenDisposition {
                input_atom_row: input_row,
                kind: HydrogenDispositionKind::FoldedInto { input_parent_row },
            });
        } else {
            let output_row = mapping.atoms.old_to_new[input_row]
                .ok_or_else(|| {
                    failure(Kind::InvariantViolation, "retained H lost its composed map")
                })?
                .index();
            let kind = isotope.map_or(
                HydrogenDispositionKind::Retained {
                    output_atom_row: output_row,
                },
                |isotope| HydrogenDispositionKind::RetainedIsotope {
                    output_atom_row: output_row,
                    isotope,
                },
            );
            dispositions.push(HydrogenDisposition {
                input_atom_row: input_row,
                kind,
            });
        }
    }
    if removed_current_rows.len()
        != current_hydrogens_count_removed(&removed.mapping.atoms.old_to_new)
    {
        return Err(failure(
            Kind::InvariantViolation,
            "hydrogen folding removal disposition count does not match core mapping",
        ));
    }
    Ok(HydrogenFoldRun {
        topology: removed.topology,
        topology_mapping: mapping,
        dispositions,
    })
}

fn current_hydrogens_count_removed(mapping: &[Option<AtomId>]) -> usize {
    mapping.iter().filter(|mapped| mapped.is_none()).count()
}

fn failure(kind: Kind, context: impl Into<String>) -> MsReadyFailure {
    MsReadyFailure {
        evidence: Default::default(),
        kind,
        stage: Stage::Normalization,
        context: context.into(),
        original: None,
    }
}
