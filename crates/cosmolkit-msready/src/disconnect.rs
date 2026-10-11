use cosmolkit_model::{
    AtomId, BondId, BondOrder, CoordinateBlock, MoleculeProperties, TopologyBlock, TopologyMapping,
};
use cosmolkit_types::Element;

use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlkaliDisconnection {
    /// Atom rows in the input topology (the disconnection never removes atoms).
    pub metal_atom: usize,
    pub oxygen_atom: usize,
    /// Bond row in the input topology for a single disconnection.
    pub bond_index: usize,
    pub metal_charge_delta: i8,
    pub oxygen_charge_delta: i8,
    pub total_charge_before: i32,
    pub total_charge_after: i32,
    pub formula_before: String,
    pub formula_after: String,
    /// Core-assigned total hydrogen counts in input atom-row order.
    pub hydrogen_counts_by_input_atom: Vec<u32>,
    /// Atom/bond mapping from the input graph to the disconnected graph.
    pub topology_mapping: TopologyMapping,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AlkaliDisconnectionRun {
    pub topology: TopologyBlock,
    pub edits: Vec<AlkaliDisconnection>,
    pub topology_mapping: TopologyMapping,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct AlkaliOxygenCandidate {
    pub metal_atom: usize,
    pub oxygen_atom: usize,
    pub bond_index: usize,
    pub metal_charge_delta: i8,
    pub oxygen_charge_delta: i8,
}

/// Find only the frozen Li/Na/K–O pairs and ligand states.
///
/// The source K/Li/Na MRVs depict one ordinary R–M bond and a separated M+;
/// Native scope permits the explicit C/P/S-oxygen and hydroxide adaptations
/// encoded here, not a general metal cleanup. Atom-row order chooses among
/// multiple eligible bonds; callers repeat this deterministic edit to a fixed
/// point without using molecular identity or atom labels to guess chemistry.
pub(crate) fn alkali_oxygen_candidate(
    topology: &TopologyBlock,
) -> Result<Option<AlkaliOxygenCandidate>, MsReadyFailure> {
    Ok(alkali_oxygen_candidates(topology)?.into_iter().next())
}

#[derive(Clone, Copy)]
enum CandidateCollectionStep {
    MetalExamined,
    ValenceAssigned,
    CandidateBuffered,
}

pub(crate) fn alkali_oxygen_candidates(
    topology: &TopologyBlock,
) -> Result<Vec<AlkaliOxygenCandidate>, MsReadyFailure> {
    alkali_oxygen_candidates_limited(topology, usize::MAX)
}

pub(crate) fn alkali_oxygen_candidates_limited(
    topology: &TopologyBlock,
    remaining: usize,
) -> Result<Vec<AlkaliOxygenCandidate>, MsReadyFailure> {
    alkali_oxygen_candidates_observed(topology, remaining, |_| {})
}

fn alkali_oxygen_candidates_observed(
    topology: &TopologyBlock,
    remaining: usize,
    mut observe: impl FnMut(CandidateCollectionStep),
) -> Result<Vec<AlkaliOxygenCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: mrv/Metanode (#1222)__MarvinSketch (#2131).mrv; SHA-256 81a8813b374a5b7f9997addac27a022d1f98e2c9d7cb0b46e595d6156b51eb1f
    // <bondArray><bond id="b1" atomRefs2="a2 a1" order="1"></bond></bondArray>
    // <atomArray atomID="a1" elementType="K" formalCharge="1" x2="-5.627368450164795" y2="3.418105363845825"></atomArray><bondArray></bondArray>
    // END MSREADY MRV SOURCE
    // BEGIN MSREADY MRV SOURCE: mrv/Metanode (#1222)__MarvinSketch (#2135).mrv; SHA-256 551ede15d62ba8a13b700ae937344a1e5512c64737d4594369927a279f758c5b
    // <bondArray><bond id="b1" atomRefs2="a2 a1" order="1"></bond></bondArray>
    // <atomArray atomID="a1" elementType="Li" formalCharge="1" x2="-5.627368450164795" y2="3.418105363845825"></atomArray><bondArray></bondArray>
    // END MSREADY MRV SOURCE
    // BEGIN MSREADY MRV SOURCE: mrv/Metanode (#1222)__MarvinSketch (#2136).mrv; SHA-256 53bdbfaa7e9a514a2d1d60fc848a9303c18cfe4842bf9454593ca574c53c37dd
    // <bondArray><bond id="b1" atomRefs2="a2 a1" order="1"></bond></bondArray>
    // <atomArray atomID="a1" elementType="Na" formalCharge="1" x2="-5.627368450164795" y2="3.418105363845825"></atomArray><bondArray></bondArray>
    // END MSREADY MRV SOURCE
    // MSReady❗❗: Behavior: Native guarded O/C/P/S or hydroxide ligands and
    // exact (0,0)/(+1,-1) charge pairs adapt the schematic R-M disconnection;
    // source R is not an arbitrary deleted ligand, and Indigo parity is unproved.
    // Complexity: adjacency scans and one lazy O(V+E) owner valence assignment
    // per collection, reused for every metal. Stop on the first excess match
    // before buffering it; no candidate products are built by this collector.
    // Source matcher cost is unknown from static MRV, so no source
    // complexity-equivalence claim is made.
    topology
        .validate()
        .map_err(|error| failure(Kind::InvalidInput, error.to_string()))?;
    let mut candidates = Vec::new();
    let mut valence = None;
    for (metal_index, metal) in topology.atoms.iter().enumerate() {
        if !matches!(metal.element(), Element::LI | Element::NA | Element::K) {
            continue;
        }
        observe(CandidateCollectionStep::MetalExamined);
        let neighbors = topology.adjacency.neighbors_of(metal_index);
        if neighbors.len() != 1
            || metal.is_aromatic()
            || metal.radical_electrons() != 0
            || metal.explicit_hydrogens() != 0
            || !metal.tracked_isotopic_hydrogens().is_empty()
        {
            continue;
        }
        let neighbor = &neighbors[0];
        let oxygen_index = neighbor.atom_index;
        let oxygen = &topology.atoms[oxygen_index];
        let bond = &topology.bonds[neighbor.bond.index()];
        if oxygen.element() != Element::O
            || oxygen.is_aromatic()
            || oxygen.radical_electrons() != 0
            || bond.order() != BondOrder::Single
            || bond.is_aromatic()
            || bond.query().is_some()
        {
            continue;
        }
        let (metal_charge_delta, oxygen_charge_delta) =
            match (metal.formal_charge(), oxygen.formal_charge()) {
                (0, 0) => (1, -1),
                (1, -1) => (0, 0),
                _ => continue,
            };
        let ligand_neighbors: Vec<_> = topology
            .adjacency
            .neighbors_of(oxygen_index)
            .iter()
            .filter(|item| item.atom_index != metal_index)
            .collect();
        if ligand_neighbors.iter().any(|item| {
            let atom = &topology.atoms[item.atom_index];
            let bond = &topology.bonds[item.bond.index()];
            bond.order() != BondOrder::Single
                || bond.is_aromatic()
                || bond.query().is_some()
                || (atom.atomic_number() == 1
                    && (atom.formal_charge() != 0 || atom.radical_electrons() != 0))
        }) {
            continue;
        }
        let heavy_neighbors: Vec<_> = ligand_neighbors
            .iter()
            .filter(|item| topology.atoms[item.atom_index].atomic_number() != 1)
            .collect();
        let organic_ligand = heavy_neighbors.len() == 1
            && matches!(
                topology.atoms[heavy_neighbors[0].atom_index].element(),
                Element::C | Element::P | Element::S
            )
            && topology.atoms[heavy_neighbors[0].atom_index].radical_electrons() == 0
            && (!topology.atoms[heavy_neighbors[0].atom_index].is_aromatic()
                || topology.atoms[heavy_neighbors[0].atom_index].element() == Element::C);
        if !organic_ligand && !heavy_neighbors.is_empty() {
            continue;
        }
        // Isotopic hydrogen labels are represented either by actual H-neighbor
        // atoms or the model's tracked-isotope carrier, and neither is erased.
        if valence.is_none() {
            observe(CandidateCollectionStep::ValenceAssigned);
            valence = Some(
                cosmolkit_core::assign_valence_with_options_for_topology(
                    topology,
                    cosmolkit_core::ValenceModel::RdkitLike,
                    false,
                )
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?,
            );
        }
        let valence = valence.as_ref().expect("assignment initialized above");
        let metal_hydrogens = atom_hydrogen_count(topology, valence, metal_index)?;
        if metal_hydrogens != 0 {
            continue;
        }
        let oxygen_hydrogens = atom_hydrogen_count(topology, valence, oxygen_index)?;
        let hydroxide_ligand = heavy_neighbors.is_empty() && oxygen_hydrogens == 1;
        if organic_ligand && oxygen_hydrogens == 0 || hydroxide_ligand {
            if candidates.len() >= remaining {
                let mut error = failure(
                    Kind::TransformationLimit,
                    format!("alkali disconnection candidate allowance {remaining} exceeded"),
                );
                error.evidence.location.rule_id = Some("alkali_oxygen_disconnection".into());
                error.evidence.location.atom_indices = vec![metal_index, oxygen_index];
                error.evidence.location.bond_indices = vec![neighbor.bond.index()];
                return Err(error);
            }
            observe(CandidateCollectionStep::CandidateBuffered);
            candidates.push(AlkaliOxygenCandidate {
                metal_atom: metal_index,
                oxygen_atom: oxygen_index,
                bond_index: neighbor.bond.index(),
                metal_charge_delta,
                oxygen_charge_delta,
            });
        }
    }
    Ok(candidates)
}

/// Apply the first eligible frozen alkali–oxygen disconnection, if present.
/// Neutral `(0,0)` inputs must pass the selected pre-edit core validation;
/// only exact `(1,-1)` depictions defer that check to post-edit fragments.
pub fn disconnect_alkali_oxygen(
    topology: &TopologyBlock,
) -> Result<Option<(TopologyBlock, AlkaliDisconnection)>, MsReadyFailure> {
    disconnect_alkali_oxygen_with_budget(topology, None)
}

fn disconnect_alkali_oxygen_with_budget(
    topology: &TopologyBlock,
    mut budget: Option<&mut crate::rewrite::TransformationBudget>,
) -> Result<Option<(TopologyBlock, AlkaliDisconnection)>, MsReadyFailure> {
    crate::validate::validate_syntax(
        topology,
        &MoleculeProperties::default(),
        &crate::MsReadyParams::default(),
    )?;
    let remaining = budget
        .as_deref()
        .map_or(usize::MAX, |b| b.remaining_matches());
    let candidates = alkali_oxygen_candidates_limited(topology, remaining)?;
    if let Some(budget) = budget.as_deref_mut() {
        budget
            .record_candidates(candidates.len())
            .map_err(|mut error| {
                error.evidence.location.rule_id = Some("alkali_oxygen_disconnection".into());
                error.evidence.location.atom_indices = candidates
                    .iter()
                    .flat_map(|c| [c.metal_atom, c.oxygen_atom])
                    .collect();
                error.evidence.location.bond_indices =
                    candidates.iter().map(|c| c.bond_index).collect();
                error
            })?;
    }
    let Some(candidate) = candidates.into_iter().next() else {
        return Ok(None);
    };
    validate_pending_repairs(topology, true, &[])?;

    let before_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let before_hydrogens = topology
        .atoms
        .iter()
        .enumerate()
        .map(|(row, _)| atom_hydrogen_count(topology, &before_valence, row))
        .collect::<Result<Vec<_>, _>>()?;
    let before_formula = formula(topology, &before_valence)?;
    let before_charge = total_charge(topology);

    let metal_id = AtomId::new(candidate.metal_atom);
    let oxygen_id = AtomId::new(candidate.oxygen_atom);
    let bond_id = BondId::new(candidate.bond_index);
    let mut edit = topology
        .begin_batch_edit()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    edit.remove_bond(bond_id)
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    if candidate.metal_charge_delta != 0 {
        edit.atom_mut(metal_id)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
            .set_formal_charge(1);
        edit.atom_mut(oxygen_id)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
            .set_formal_charge(-1);
    }
    let (mut changed, mut mapping) = edit
        .finish()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;

    // If disconnecting the neutral metal makes an implicit hydroxide H vanish,
    // transport that same H into the explicit count. This is representation
    // transport, not proton transfer; the per-atom H and isotope invariants
    // below are authoritative.
    let after_initial_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &changed,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let after_initial_oxygen_h =
        atom_hydrogen_count(&changed, &after_initial_valence, candidate.oxygen_atom)?;
    if after_initial_oxygen_h > before_hydrogens[candidate.oxygen_atom] {
        return Err(failure(
            Kind::InvariantViolation,
            "disconnection increased oxygen hydrogen count",
        ));
    }
    let missing_oxygen_h = before_hydrogens[candidate.oxygen_atom] - after_initial_oxygen_h;
    if missing_oxygen_h != 0 {
        let oxygen = &changed.atoms[candidate.oxygen_atom];
        let inferred =
            after_initial_valence.implicit_hydrogens[candidate.oxygen_atom].max(0) as u32;
        let explicit_neighbors = changed
            .adjacency
            .neighbors_of(candidate.oxygen_atom)
            .iter()
            .filter(|item| changed.atoms[item.atom_index].atomic_number() == 1)
            .count() as u32;
        let tracked_isotopic = oxygen.tracked_isotopic_hydrogens().len() as u32;
        let required_explicit = before_hydrogens[candidate.oxygen_atom]
            .checked_sub(inferred + explicit_neighbors + tracked_isotopic)
            .ok_or_else(|| {
                failure(
                    Kind::InvariantViolation,
                    "post-edit implicit/neighbour H exceeds the original oxygen H count",
                )
            })?;
        let required_explicit = u8::try_from(required_explicit).map_err(|_| {
            failure(
                Kind::InvariantViolation,
                "oxygen explicit-H transport exceeds model capacity",
            )
        })?;
        let old_explicit = oxygen.explicit_hydrogens();
        let mut h_edit = changed
            .begin_batch_edit()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        h_edit
            .atom_mut(oxygen_id)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
            .set_explicit_hydrogens(required_explicit);
        let (next, h_mapping) = h_edit
            .finish()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        changed = next;
        mapping = compose_mappings(&mapping, &h_mapping)?;
        let explicit_delta = i16::from(required_explicit) - i16::from(old_explicit);
        if explicit_delta < 0 {
            return Err(failure(
                Kind::InvariantViolation,
                "implicit-H transport unexpectedly reduced explicit hydrogen count",
            ));
        }
    }

    // Validate all repaired fragments now, deferring only exact pending charged
    // pairs in other fragments (or in the same connected ligand), and exact
    // exclusion-only components. The invariant assignment is not a bypass:
    // strict owner checks run on those bounded validation graphs first.
    validate_pending_repairs(
        &changed,
        true,
        &[candidate.metal_atom, candidate.oxygen_atom],
    )?;
    let after_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &changed,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let after_formula = formula(&changed, &after_valence)?;
    let after_hydrogens = changed
        .atoms
        .iter()
        .enumerate()
        .map(|(row, _)| atom_hydrogen_count(&changed, &after_valence, row))
        .collect::<Result<Vec<_>, _>>()?;
    let after_charge = total_charge(&changed);
    if before_formula != after_formula
        || before_charge != after_charge
        || before_hydrogens != after_hydrogens
    {
        return Err(failure(
            Kind::InvariantViolation,
            format!(
                "alkali disconnection changed formula/charge/H: {before_formula} ({before_charge}) -> {after_formula} ({after_charge})"
            ),
        ));
    }
    for (row, (before, after)) in topology.atoms.iter().zip(&changed.atoms).enumerate() {
        let mut expected = before.clone();
        if row == candidate.metal_atom || row == candidate.oxygen_atom {
            expected.set_formal_charge(after.formal_charge());
        }
        if row == candidate.oxygen_atom {
            expected.set_explicit_hydrogens(after.explicit_hydrogens());
        }
        if expected != *after {
            return Err(failure(
                Kind::InvariantViolation,
                format!("alkali disconnection changed protected atom state at row {row}"),
            ));
        }
    }
    for (old_row, mapped) in mapping.bonds.old_to_new.iter().enumerate() {
        match mapped {
            None if old_row == candidate.bond_index => {}
            None => {
                return Err(failure(
                    Kind::InvariantViolation,
                    format!("alkali disconnection removed unrelated bond row {old_row}"),
                ));
            }
            Some(new_id) => {
                let before = &topology.bonds[old_row];
                let after = &changed.bonds[new_id.index()];
                let before_properties = before
                    .property_records(true, true)
                    .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
                let after_properties = after
                    .property_records(true, true)
                    .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
                if before.begin() != after.begin()
                    || before.end() != after.end()
                    || before.order() != after.order()
                    || before.is_aromatic() != after.is_aromatic()
                    || before.is_conjugated() != after.is_conjugated()
                    || before.direction() != after.direction()
                    || before.stereo() != after.stereo()
                    || before.stereo_atoms() != after.stereo_atoms()
                    || before.unknown_stereo() != after.unknown_stereo()
                    || before.query() != after.query()
                    || !before_properties.eq(after_properties)
                {
                    return Err(failure(
                        Kind::InvariantViolation,
                        format!("alkali disconnection changed retained bond row {old_row}"),
                    ));
                }
            }
        }
    }
    mapping
        .validate_for_counts(
            topology.atoms.len(),
            changed.atoms.len(),
            topology.bonds.len(),
            changed.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;

    if let Some(budget) = budget {
        budget.commit_edit().map_err(|mut error| {
            error.evidence.location.rule_id = Some("alkali_oxygen_disconnection".into());
            error.evidence.location.atom_indices =
                vec![candidate.metal_atom, candidate.oxygen_atom];
            error.evidence.location.bond_indices = vec![candidate.bond_index];
            error
        })?;
    }
    Ok(Some((
        changed,
        AlkaliDisconnection {
            metal_atom: candidate.metal_atom,
            oxygen_atom: candidate.oxygen_atom,
            bond_index: candidate.bond_index,
            metal_charge_delta: candidate.metal_charge_delta,
            oxygen_charge_delta: candidate.oxygen_charge_delta,
            total_charge_before: before_charge,
            total_charge_after: after_charge,
            formula_before: before_formula,
            formula_after: after_formula,
            hydrogen_counts_by_input_atom: before_hydrogens,
            topology_mapping: mapping,
        },
    )))
}

/// Apply every eligible bond in deterministic source-atom order and return one
/// composed input-to-output mapping plus source-row-addressed edit records.
pub fn disconnect_all_alkali_oxygen(
    topology: &TopologyBlock,
) -> Result<AlkaliDisconnectionRun, MsReadyFailure> {
    disconnect_all_alkali_oxygen_observed(topology, None, |_, _| Ok(()))
}

pub(crate) fn disconnect_all_alkali_oxygen_observed(
    topology: &TopologyBlock,
    mut budget: Option<&mut crate::rewrite::TransformationBudget>,
    mut observe: impl FnMut(&TopologyBlock, &AlkaliDisconnection) -> Result<(), MsReadyFailure>,
) -> Result<AlkaliDisconnectionRun, MsReadyFailure> {
    let mut changed = topology.clone();
    let mut edits = Vec::new();
    let mut input_atom_rows: Vec<_> = (0..topology.atoms.len()).collect();
    let mut input_bond_rows: Vec<_> = (0..topology.bonds.len()).collect();
    loop {
        let Some((next, mut edit)) =
            disconnect_alkali_oxygen_with_budget(&changed, budget.as_deref_mut()).map_err(
                |mut error| {
                    error.evidence.location.atom_indices = error
                        .evidence
                        .location
                        .atom_indices
                        .iter()
                        .filter_map(|row| input_atom_rows.get(*row).copied())
                        .collect();
                    error.evidence.location.bond_indices = error
                        .evidence
                        .location
                        .bond_indices
                        .iter()
                        .filter_map(|row| input_bond_rows.get(*row).copied())
                        .collect();
                    error
                },
            )?
        else {
            break;
        };
        edit.metal_atom = input_atom_rows[edit.metal_atom];
        edit.oxygen_atom = input_atom_rows[edit.oxygen_atom];
        edit.bond_index = input_bond_rows[edit.bond_index];
        let next_atom_rows = edit
            .topology_mapping
            .atoms
            .new_to_old
            .iter()
            .map(|old| old.map(|id| input_atom_rows[id.index()]))
            .collect::<Vec<_>>();
        let next_bond_rows = edit
            .topology_mapping
            .bonds
            .new_to_old
            .iter()
            .map(|old| old.map(|id| input_bond_rows[id.index()]))
            .collect::<Vec<_>>();
        input_atom_rows = next_atom_rows.into_iter().flatten().collect();
        input_bond_rows = next_bond_rows.into_iter().flatten().collect();
        edit.topology_mapping =
            mapping_from_rows(topology, &next, &input_atom_rows, &input_bond_rows)?;
        observe(&next, &edit)?;
        changed = next;
        edits.push(edit);
    }
    validate_pending_repairs(&changed, false, &[])?;
    let topology_mapping =
        mapping_from_rows(topology, &changed, &input_atom_rows, &input_bond_rows)?;
    Ok(AlkaliDisconnectionRun {
        topology: changed,
        edits,
        topology_mapping,
    })
}

/// Strict validation with only frozen charged-alkali and neutral-nitro
/// representation repairs deferred to their ordered actual stages.
/// Exclusion-only states are recognized by the whole-component exact matcher,
/// never by an element/formula shortcut. Validation copies do not replace input.
pub(crate) fn validate_pending_repairs(
    topology: &TopologyBlock,
    allow_pending: bool,
    repaired_rows: &[usize],
) -> Result<(), MsReadyFailure> {
    let parts = crate::components::split_mapped_components(
        topology,
        &CoordinateBlock::default(),
        &MoleculeProperties::default(),
    )?;
    for (index, part) in parts.iter().enumerate() {
        let annotate = |mut error: MsReadyFailure| {
            error.evidence.location.component_index = Some(index);
            error.evidence.location.atom_indices =
                if error.evidence.location.atom_indices.is_empty() {
                    part.input_atom_indices.clone()
                } else {
                    error
                        .evidence
                        .location
                        .atom_indices
                        .iter()
                        .filter_map(|row| part.input_atom_indices.get(*row).copied())
                        .collect()
                };
            error
        };
        let candidates = alkali_oxygen_candidates(&part.topology).map_err(annotate)?;
        let repaired_fragment = part
            .input_atom_indices
            .iter()
            .any(|row| repaired_rows.contains(row));
        if candidates.is_empty()
            && !repaired_fragment
            && !crate::exclusions::match_exact_exclusions(&part.topology)?.is_empty()
        {
            continue;
        }
        if !allow_pending && !candidates.is_empty() {
            return Err(failure(
                Kind::InvariantViolation,
                "completed disconnection batch contains pending repairs",
            ));
        }
        let mut strict_edit = part
            .topology
            .begin_batch_edit()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        for candidate in &candidates {
            if candidate.metal_charge_delta == 0 {
                strict_edit
                    .remove_bond(BondId::new(candidate.bond_index))
                    .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            }
        }
        let (strict_graph, _) = strict_edit
            .finish()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        // Neutral nitro is repaired later by the ordered actual group stage;
        // only its exact adopted guards may defer strict N valence here.
        let strict_graph = crate::rewrite::nitro_validation_graph(&strict_graph)?;
        crate::validate::strict_owner_valence(&strict_graph).map_err(|mut error| {
            error.context = format!("disconnection component {index}: {}", error.context);
            annotate(error)
        })?;
        // Pending neutral pairs passed ordinary pre-valence above. For domain
        // checks only, project their declared metal spectator states so no
        // unrelated retained metal can hide behind the repair deferral.
        let mut domain_edit = strict_graph
            .begin_batch_edit()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        for candidate in alkali_oxygen_candidates(&strict_graph)? {
            domain_edit
                .remove_bond(BondId::new(candidate.bond_index))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            domain_edit
                .atom_mut(AtomId::new(candidate.metal_atom))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
                .set_formal_charge(1);
            domain_edit
                .atom_mut(AtomId::new(candidate.oxygen_atom))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
                .set_formal_charge(-1);
        }
        let (domain_graph, _) = domain_edit
            .finish()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        crate::validate::retained_metal_states(&domain_graph).map_err(annotate)?;
        if candidates.is_empty() {
            crate::validate::validate_input(
                &strict_graph,
                &MoleculeProperties::default(),
                &Default::default(),
            )
            .map_err(annotate)?;
        }
    }
    Ok(())
}

fn mapping_from_rows(
    old: &TopologyBlock,
    new: &TopologyBlock,
    input_atom_rows: &[usize],
    input_bond_rows: &[usize],
) -> Result<TopologyMapping, MsReadyFailure> {
    let mut atoms_old_to_new = vec![None; old.atoms.len()];
    let atoms_new_to_old: Vec<_> = input_atom_rows
        .iter()
        .copied()
        .enumerate()
        .map(|(new_row, old_row)| {
            atoms_old_to_new[old_row] = Some(AtomId::new(new_row));
            Some(AtomId::new(old_row))
        })
        .collect();
    let mut bonds_old_to_new = vec![None; old.bonds.len()];
    let bonds_new_to_old: Vec<_> = input_bond_rows
        .iter()
        .copied()
        .enumerate()
        .map(|(new_row, old_row)| {
            bonds_old_to_new[old_row] = Some(BondId::new(new_row));
            Some(BondId::new(old_row))
        })
        .collect();
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
            old.atoms.len(),
            new.atoms.len(),
            old.bonds.len(),
            new.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    Ok(mapping)
}

fn compose_mappings(
    first: &TopologyMapping,
    second: &TopologyMapping,
) -> Result<TopologyMapping, MsReadyFailure> {
    let atoms_old_to_new = first
        .atoms
        .old_to_new
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|id| second.atoms.old_to_new.get(id.index()).copied().flatten())
        })
        .collect();
    let atoms_new_to_old = second
        .atoms
        .new_to_old
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|id| first.atoms.new_to_old.get(id.index()).copied().flatten())
        })
        .collect();
    let bonds_old_to_new = first
        .bonds
        .old_to_new
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|id| second.bonds.old_to_new.get(id.index()).copied().flatten())
        })
        .collect();
    let bonds_new_to_old = second
        .bonds
        .new_to_old
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|id| first.bonds.new_to_old.get(id.index()).copied().flatten())
        })
        .collect();
    Ok(TopologyMapping {
        atoms: cosmolkit_model::AtomMapping {
            old_to_new: atoms_old_to_new,
            new_to_old: atoms_new_to_old,
        },
        bonds: cosmolkit_model::BondMapping {
            old_to_new: bonds_old_to_new,
            new_to_old: bonds_new_to_old,
        },
    })
}

fn atom_hydrogen_count(
    topology: &TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
    atom: usize,
) -> Result<u32, MsReadyFailure> {
    let count = cosmolkit_core::total_hydrogen_count(topology, valence, AtomId::new(atom), true)
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    count
        .checked_add(topology.atoms[atom].tracked_isotopic_hydrogens().len() as u32)
        .ok_or_else(|| failure(Kind::InvariantViolation, "hydrogen count overflow"))
}

fn formula(
    topology: &TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
) -> Result<String, MsReadyFailure> {
    cosmolkit_descriptors::molecular_formula_with_valence(topology, true, false, Some(valence))
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))
}

fn total_charge(topology: &TopologyBlock) -> i32 {
    topology
        .atoms
        .iter()
        .map(|atom| i32::from(atom.formal_charge()))
        .sum()
}

fn failure(kind: Kind, context: impl Into<String>) -> MsReadyFailure {
    MsReadyFailure {
        evidence: Default::default(),
        kind,
        stage: Stage::Validation,
        context: context.into(),
        original: None,
    }
}

#[cfg(test)]
mod component_tests {
    use super::*;
    use crate::{match_exact_exclusions, validate_input};
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondQueryPredicate, BondSpec, QueryNode};
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

    fn parsed(smiles: &str) -> TopologyBlock {
        parse_smiles(smiles, &SmilesParseParams::default())
            .unwrap_or_else(|error| panic!("{smiles}: {error}"))
            .topology
    }

    #[test]
    fn candidate_collector_stops_at_remaining_allowance_before_buffering_or_repeated_owner_work() {
        for ligand in ["CO[Na]", "C[O-][Na+]"] {
            let source = parsed(
                &std::iter::repeat_n(ligand, 64)
                    .collect::<Vec<_>>()
                    .join("."),
            );
            for remaining in [0, 1, 63, 64] {
                let mut examined = 0;
                let mut assignments = 0;
                let mut buffered = 0;
                let result =
                    alkali_oxygen_candidates_observed(&source, remaining, |step| match step {
                        CandidateCollectionStep::MetalExamined => examined += 1,
                        CandidateCollectionStep::ValenceAssigned => assignments += 1,
                        CandidateCollectionStep::CandidateBuffered => buffered += 1,
                    });
                assert_eq!(
                    buffered, remaining,
                    "{ligand}: collection must not buffer beyond the cap"
                );
                assert_eq!(
                    examined,
                    (remaining + 1).min(64),
                    "{ligand}: stop on the first excess candidate"
                );
                assert_eq!(
                    assignments, 1,
                    "{ligand}: reuse the topology owner assignment"
                );
                if remaining < 64 {
                    assert_eq!(result.unwrap_err().kind, Kind::TransformationLimit);
                } else {
                    assert_eq!(result.unwrap().len(), 64, "the exact boundary is inclusive");
                }
            }
        }
    }

    #[test]
    fn complete_repair_batch_preserves_pending_charged_and_exclusion_only_components() {
        for smiles in [
            "CC(=O)[O-][Li+].CC(=O)[O-][Na+]",
            "CC(=O)O[Li].CC(=O)[O-][Na+]",
            "CC(=O)[O-][Li+].CC(=O)O[Na]",
            "CC(=O)[O-][Li+].CC(=O)[O-][Na+].ON(=O)O.P(F)(F)(F)(F)(F)F.[Zn]",
            "[Li+][O-]C(=O)CC(=O)[O-][Na+]",
        ] {
            let source = parsed(smiles);
            let snapshot = source.clone();
            let result = disconnect_all_alkali_oxygen(&source)
                .unwrap_or_else(|error| panic!("{smiles}: {error}"));
            assert_eq!(result.edits.len(), 2, "{smiles}");
            assert_eq!(source, snapshot);
            for event in &result.edits {
                assert_eq!(event.formula_before, event.formula_after);
                assert_eq!(event.total_charge_before, event.total_charge_after);
            }
            let parts = crate::components::split_mapped_components(
                &result.topology,
                &Default::default(),
                &Default::default(),
            )
            .unwrap();
            for part in parts {
                if crate::match_exact_exclusions(&part.topology)
                    .unwrap()
                    .is_empty()
                {
                    crate::validate::strict_owner_valence(&part.topology).unwrap();
                    crate::validate_input(&part.topology, &Default::default(), &Default::default())
                        .unwrap();
                }
            }
        }
    }

    #[test]
    fn disconnection_defers_only_guarded_neutral_nitro_to_its_actual_rule_stage() {
        for (smiles, expected_cuts, expected_nitros) in [
            ("CN(=O)=O", 0, 1),
            ("CN(=O)=O.CN(=O)=O", 0, 2),
            ("CC(=O)[O-][Li+].CN(=O)=O", 1, 1),
            ("CC(=O)O[Li].CC(=O)[O-][Na+].CN(=O)=O.CN(=O)=O", 2, 2),
        ] {
            let source = parsed(smiles);
            let snapshot = source.clone();
            let result = disconnect_all_alkali_oxygen(&source)
                .unwrap_or_else(|error| panic!("{smiles}: {error}"));
            assert_eq!(source, snapshot);
            assert_eq!(result.edits.len(), expected_cuts);
            let neutral_nitros = result
                .topology
                .atoms
                .iter()
                .filter(|atom| atom.element() == Element::N && atom.formal_charge() == 0)
                .count();
            assert_eq!(
                neutral_nitros, expected_nitros,
                "validation must not silently normalize analysis nitro"
            );
            let normalized =
                crate::rewrite::normalize_group_rules(&result.topology, &Default::default())
                    .unwrap();
            assert_eq!(
                normalized.edits.len(),
                expected_nitros,
                "later actual rule logs each repair"
            );
            crate::validate::strict_owner_valence(&normalized.topology).unwrap();
        }
        for smiles in ["CN(=O)(=O)C", "C[NH](=O)=O", "O[15N](=O)O"] {
            let source = parsed(smiles);
            assert_eq!(
                disconnect_all_alkali_oxygen(&source).unwrap_err().kind,
                Kind::InvalidValence,
                "{smiles}"
            );
        }
        let mut radical = parsed("CN(=O)=O");
        radical.atoms[1].set_radical_electrons(1);
        assert_eq!(
            disconnect_all_alkali_oxygen(&radical).unwrap_err().kind,
            Kind::InvalidValence
        );
    }

    #[test]
    fn pending_repairs_do_not_bypass_unrelated_invalid_or_unmodeled_states() {
        for (smiles, kind) in [
            ("CC(=O)[O-][Li+].C[Fe]C", Kind::UnsupportedChemistry),
            ("CC(=O)[O-][Li+].CN(=O)(=O)C", Kind::InvalidValence),
        ] {
            let source = parsed(smiles);
            assert_eq!(
                disconnect_all_alkali_oxygen(&source).unwrap_err().kind,
                kind,
                "{smiles}"
            );
        }
    }

    fn assert_repair(smiles: &str, deltas: (i8, i8)) -> TopologyBlock {
        let source = parsed(smiles);
        let snapshot = source.clone();
        let (changed, event) = disconnect_alkali_oxygen(&source)
            .unwrap_or_else(|error| panic!("{smiles}: {error}"))
            .unwrap_or_else(|| panic!("expected supported disconnection: {smiles}"));
        assert_eq!(source, snapshot, "input mutation: {smiles}");
        assert_eq!(
            (event.metal_charge_delta, event.oxygen_charge_delta),
            deltas
        );
        assert_eq!(event.total_charge_before, event.total_charge_after);
        assert_eq!(event.total_charge_before, total_charge(&source));
        assert_eq!(event.total_charge_after, total_charge(&changed));
        assert_eq!(event.formula_before, event.formula_after, "{smiles}");
        assert_eq!(
            event.hydrogen_counts_by_input_atom.len(),
            source.atoms.len()
        );
        assert_eq!(source.atoms.len(), changed.atoms.len());
        assert_eq!(source.bonds.len(), changed.bonds.len() + 1);
        assert!(event.topology_mapping.bonds.old_to_new[event.bond_index].is_none());
        assert!(
            event
                .topology_mapping
                .atoms
                .old_to_new
                .iter()
                .all(Option::is_some)
        );
        let valence = crate::validate::strict_owner_valence(&changed)
            .unwrap_or_else(|error| panic!("{smiles}: {error}"));
        for row in 0..source.atoms.len() {
            assert_eq!(
                event.hydrogen_counts_by_input_atom[row],
                atom_hydrogen_count(&changed, &valence, row).unwrap(),
                "{smiles}, atom {row}"
            );
        }
        let parts = crate::components::split_mapped_components(
            &changed,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
        )
        .unwrap();
        assert_eq!(parts.len(), 2, "{smiles}");
        for part in parts {
            crate::validate::strict_owner_valence(&part.topology)
                .unwrap_or_else(|error| panic!("{smiles}: {error}"));
        }
        for (before, after) in source.atoms.iter().zip(&changed.atoms) {
            assert_eq!(before.element(), after.element(), "{smiles}");
            assert_eq!(before.isotope(), after.isotope(), "{smiles}");
            assert_eq!(
                before.tracked_isotopic_hydrogens(),
                after.tracked_isotopic_hydrogens(),
                "{smiles}"
            );
            assert_eq!(
                before.radical_electrons(),
                after.radical_electrons(),
                "{smiles}"
            );
        }
        changed
    }

    fn assert_not_repaired(smiles: &str) {
        let source = parsed(smiles);
        let snapshot = source.clone();
        match disconnect_alkali_oxygen(&source) {
            Ok(None) | Err(_) => {}
            Ok(Some((_, event))) => {
                panic!("unsupported graph was disconnected: {smiles}: {event:?}")
            }
        }
        assert_eq!(source, snapshot, "input mutation: {smiles}");
    }

    fn assert_not_repaired_topology(source: &TopologyBlock) {
        let snapshot = source.clone();
        match disconnect_alkali_oxygen(source) {
            Ok(None) | Err(_) => {}
            Ok(Some((_, event))) => panic!("unsupported graph was disconnected: {event:?}"),
        }
        assert_eq!(*source, snapshot);
    }

    #[test]
    fn all_three_metals_cover_both_charge_pairs_and_both_ligand_classes() {
        for metal in ["Li", "Na", "K"] {
            assert_repair(&format!("CC(=O)O[{metal}]"), (1, -1));
            assert_repair(&format!("CC(=O)[O-][{metal}+]"), (0, 0));
            assert_repair(&format!("CO[{metal}]"), (1, -1));
            assert_repair(&format!("O[{metal}]"), (1, -1));
            assert_repair(&format!("[OH-][{metal}+]"), (0, 0));
        }
        assert_repair("O=P(O)(O)O[Li]", (1, -1));
        assert_repair("O=S(=O)(O)O[Na]", (1, -1));
    }

    #[test]
    fn unsupported_charge_pairs_and_connectivity_are_not_repaired() {
        for metal in ["Li", "Na", "K"] {
            for smiles in [
                format!("CC(=O)O[{metal}+]"),
                format!("CC(=O)[O-][{metal}]"),
                format!("O[{metal}]O"),
                format!("C[{metal}]"),
                format!("N[{metal}]"),
                format!("O=[{metal}]"),
                format!("[H][{metal}]"),
                format!("CC(=O)[OH][{metal}]"),
            ] {
                assert_not_repaired(&smiles);
            }
        }
        assert_not_repaired("CO([Li])P");
        assert_repair("c1ccccc1O[Li]", (1, -1));
    }

    #[test]
    fn isotope_hydrogens_and_metal_labels_are_preserved_and_not_erased_for_exclusions() {
        for smiles in [
            "[2H]O[Li]",
            "[2H]O[Na]",
            "[2H]O[K]",
            "[6Li]O",
            "[23Na]O",
            "[41K]O",
        ] {
            let changed = assert_repair(smiles, (1, -1));
            let components = crate::components::split_mapped_components(
                &changed,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            )
            .unwrap();
            let metal = components
                .iter()
                .find(|part| {
                    part.topology
                        .atoms
                        .iter()
                        .any(|a| matches!(a.element(), Element::LI | Element::NA | Element::K))
                })
                .unwrap();
            if smiles.starts_with("[2H]") {
                assert!(
                    match_exact_exclusions(&metal.topology)
                        .unwrap()
                        .iter()
                        .any(|entry| matches!(entry.source_id.as_str(), "70" | "80" | "67"))
                );
                let hydroxide = components
                    .iter()
                    .find(|part| {
                        part.topology
                            .atoms
                            .iter()
                            .any(|a| a.element() == Element::O)
                    })
                    .unwrap();
                assert!(
                    match_exact_exclusions(&hydroxide.topology)
                        .unwrap()
                        .is_empty(),
                    "isotopic H was erased: {smiles}"
                );
            } else {
                assert!(
                    match_exact_exclusions(&metal.topology).unwrap().is_empty(),
                    "isotopic metal matched unlabelled ion: {smiles}"
                );
                let hydroxide = components
                    .iter()
                    .find(|part| {
                        part.topology
                            .atoms
                            .iter()
                            .any(|a| a.element() == Element::O)
                    })
                    .unwrap();
                assert!(
                    match_exact_exclusions(&hydroxide.topology)
                        .unwrap()
                        .iter()
                        .any(|entry| entry.source_id == "143")
                );
            }
        }
    }

    #[test]
    fn ionic_spectators_and_cut_hydroxides_retain_both_original_exclusion_origins() {
        for (metal, source_id) in [("Li", "70"), ("Na", "80"), ("K", "67")] {
            let dotted = parsed(&format!("CC(=O)[O-].[{metal}+]"));
            assert!(disconnect_alkali_oxygen(&dotted).unwrap().is_none());
            let parts = crate::components::split_mapped_components(
                &dotted,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            )
            .unwrap();
            let ion = parts
                .iter()
                .find(|part| part.topology.atoms.len() == 1)
                .unwrap();
            assert!(
                match_exact_exclusions(&ion.topology)
                    .unwrap()
                    .iter()
                    .any(|entry| entry.source_id == source_id)
            );

            let run = disconnect_all_alkali_oxygen(&parsed(&format!("O[{metal}]"))).unwrap();
            assert_eq!(run.edits.len(), 1);
            let original_rows: Vec<_> = run
                .topology_mapping
                .atoms
                .new_to_old
                .iter()
                .map(|id| id.unwrap().index())
                .collect();
            let cut_parts = crate::components::split_mapped_components_from_source_rows(
                &run.topology,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
                &original_rows,
            )
            .unwrap();
            assert_eq!(cut_parts.len(), 2);
            let matches: Vec<_> = cut_parts
                .iter()
                .flat_map(|part| {
                    match_exact_exclusions(&part.topology)
                        .unwrap()
                        .into_iter()
                        .map(move |entry| (part.input_atom_indices.clone(), entry))
                })
                .collect();
            assert!(
                matches
                    .iter()
                    .any(|(rows, entry)| rows.len() == 1 && entry.source_id == "143")
            );
            assert!(
                matches
                    .iter()
                    .any(|(rows, entry)| rows.len() == 1 && entry.source_id == source_id)
            );
        }
    }

    #[test]
    fn repeated_permuted_disconnections_compose_original_atom_and_bond_rows() {
        for smiles in ["O[Li].O[Na].O[K]", "[K]O.[Na]O.[Li]O"] {
            let source = parsed(smiles);
            let run = disconnect_all_alkali_oxygen(&source).unwrap();
            assert_eq!(run.edits.len(), 3, "{smiles}");
            assert_eq!(source.bonds.len(), run.topology.bonds.len() + 3);
            assert!(
                run.topology_mapping
                    .atoms
                    .old_to_new
                    .iter()
                    .all(Option::is_some)
            );
            let edited_bonds: std::collections::BTreeSet<_> =
                run.edits.iter().map(|edit| edit.bond_index).collect();
            assert_eq!(edited_bonds.len(), 3);
            for event in &run.edits {
                assert_eq!(event.total_charge_before, event.total_charge_after);
                assert_eq!(event.formula_before, event.formula_after);
                assert!(event.topology_mapping.bonds.old_to_new[event.bond_index].is_none());
                assert!(
                    event
                        .topology_mapping
                        .atoms
                        .old_to_new
                        .iter()
                        .all(Option::is_some)
                );
            }
            let parts = crate::components::split_mapped_components(
                &run.topology,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            )
            .unwrap();
            assert_eq!(parts.len(), 6);
            let all_rows: std::collections::BTreeSet<_> = parts
                .iter()
                .flat_map(|part| part.input_atom_indices.iter().copied())
                .collect();
            assert_eq!(all_rows, (0..source.atoms.len()).collect());
            let exclusion_origins: Vec<_> = parts
                .iter()
                .flat_map(|part| {
                    match_exact_exclusions(&part.topology)
                        .unwrap()
                        .into_iter()
                        .map(|entry| (entry.source_id, part.input_atom_indices.clone()))
                })
                .collect();
            assert_eq!(exclusion_origins.len(), 6);
        }
    }

    #[test]
    fn retained_hg_sn_domain_is_preserved_and_unmodeled_states_are_not_edited() {
        for smiles in [
            "C[Hg+]",
            "C[Hg]C",
            "C[Sn](C)(C)C",
            "[13C][Hg+]",
            "C[203Hg]C",
            "[119Sn](C)(C)(C)C",
        ] {
            let source = parsed(smiles);
            let snapshot = source.clone();
            assert!(disconnect_alkali_oxygen(&source).unwrap().is_none());
            assert_eq!(source, snapshot);
        }
        for smiles in ["C[Hg]", "C[Hg+2]", "C[Hg]Cl", "C[Sn](C)C"] {
            let source = parsed(smiles);
            assert!(disconnect_alkali_oxygen(&source).unwrap().is_none());
            assert!(
                validate_input(
                    &source,
                    &MoleculeProperties::default(),
                    &crate::MsReadyParams::default()
                )
                .is_err()
            );
        }
    }

    #[test]
    fn graph_level_radical_and_dative_multiple_query_bonds_are_rejected() {
        let mut radical_metal = parsed("CC(=O)O[Li]");
        let li = radical_metal
            .atoms
            .iter()
            .position(|atom| atom.element() == Element::LI)
            .unwrap();
        radical_metal.atoms[li].set_radical_electrons(1);
        assert_not_repaired_topology(&radical_metal);
        let mut radical_oxygen = parsed("CC(=O)O[Li]");
        let oxygen = radical_oxygen
            .atoms
            .iter()
            .position(|atom| {
                atom.element() == Element::O
                    && radical_oxygen
                        .adjacency
                        .neighbors_of(atom.id().index())
                        .iter()
                        .any(|n| radical_oxygen.atoms[n.atom_index].element() == Element::LI)
            })
            .unwrap();
        radical_oxygen.atoms[oxygen].set_radical_electrons(1);
        assert_not_repaired_topology(&radical_oxygen);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::O)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::LI)),
        ];
        for order in [
            BondOrder::Dative,
            BondOrder::Double,
            BondOrder::Zero,
            BondOrder::Unspecified,
        ] {
            let bond = Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), order),
            );
            let topology =
                TopologyBlock::try_from_parts(atoms.clone(), vec![bond], vec![], vec![]).unwrap();
            assert_not_repaired_topology(&topology);
        }
        let query_bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_query(QueryNode::predicate(BondQueryPredicate::Any)),
        );
        let query = TopologyBlock::try_from_parts(atoms, vec![query_bond], vec![], vec![]).unwrap();
        assert_not_repaired_topology(&query);
    }
}
