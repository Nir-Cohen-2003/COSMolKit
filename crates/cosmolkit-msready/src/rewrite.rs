use std::collections::BTreeSet;

use cosmolkit_model::{
    AtomId, AtomMapping, BondDirection, BondId, BondMapping, BondStereo, CoordinateBlock,
    MoleculeProperties, TopologyBlock, TopologyMapping,
};
use cosmolkit_smiles::{SmilesRecord, SmilesWriteParams, write_smiles_with_params};
use cosmolkit_types::ChiralTag;

use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyParams, MsReadyStage as Stage};

use crate::rules::{
    HydrogenCarrier, HydrogenTransfer, RuleCandidate, RuleId, assign_for_match,
    find_candidates_limited, hydrogen_count, validate_rule_data,
};

/// One input-wide edit budget shared by group and proton stages. Candidate
/// counts are reset only at an explicitly started sweep, never per rematch.
pub(crate) struct TransformationBudget {
    max_edits: usize,
    max_matches: usize,
    edits: usize,
    matches: usize,
}

impl TransformationBudget {
    pub(crate) fn new(params: &MsReadyParams) -> Self {
        Self {
            max_edits: params.max_edits,
            max_matches: params.max_matches_per_sweep,
            edits: 0,
            matches: 0,
        }
    }

    pub(crate) fn begin_sweep(&mut self) {
        self.matches = 0;
    }

    pub(crate) fn remaining_matches(&self) -> usize {
        self.max_matches.saturating_sub(self.matches)
    }

    pub(crate) fn record_candidates(&mut self, count: usize) -> Result<(), MsReadyFailure> {
        let next = self
            .matches
            .checked_add(count)
            .ok_or_else(|| failure(Kind::TransformationLimit, "candidate count overflow"))?;
        if next > self.max_matches {
            return Err(failure(
                Kind::TransformationLimit,
                format!("candidate sweep limit {} exceeded", self.max_matches),
            ));
        }
        self.matches = next;
        Ok(())
    }

    pub(crate) fn commit_edit(&mut self) -> Result<(), MsReadyFailure> {
        if self.edits >= self.max_edits {
            return Err(failure(
                Kind::TransformationLimit,
                format!("input edit limit {} reached", self.max_edits),
            ));
        }
        self.edits += 1;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RuleEdit {
    pub rule_id: &'static str,
    pub appearance_order: usize,
    pub matched_atom_rows: Vec<usize>,
    pub matched_bond_rows: Vec<usize>,
    pub hydrogen_transfer: Option<HydrogenTransfer>,
    pub formula_before: String,
    pub formula_after: String,
    pub total_charge_before: i32,
    pub total_charge_after: i32,
    /// Original-input-to-topology mapping immediately after this edit.
    pub topology_mapping: TopologyMapping,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RuleRewriteRun {
    pub topology: TopologyBlock,
    /// Original-input-to-final-topology atom and bond mapping.
    pub topology_mapping: TopologyMapping,
    pub edits: Vec<RuleEdit>,
}

#[derive(Clone)]
struct ValidCandidate {
    canonical: String,
    tie_key: (Vec<usize>, Vec<usize>, String),
    candidate: RuleCandidate,
    topology: TopologyBlock,
    step_mapping: TopologyMapping,
    formula_before: String,
    formula_after: String,
    charge_before: i32,
    charge_after: i32,
}

/// Apply the six frozen source appearances until the graph reaches a fixed point.
pub(crate) fn normalize_group_rules(
    topology: &TopologyBlock,
    params: &MsReadyParams,
) -> Result<RuleRewriteRun, MsReadyFailure> {
    let mut budget = TransformationBudget::new(params);
    normalize_group_rules_with_budget(topology, params, &mut budget)
}

pub(crate) fn normalize_group_rules_with_budget(
    topology: &TopologyBlock,
    params: &MsReadyParams,
    budget: &mut TransformationBudget,
) -> Result<RuleRewriteRun, MsReadyFailure> {
    normalize_group_rules_observed(topology, params, budget, |_| Ok(()))
}

pub(crate) fn normalize_group_rules_observed(
    topology: &TopologyBlock,
    params: &MsReadyParams,
    budget: &mut TransformationBudget,
    observe: impl FnMut(&RuleEdit) -> Result<(), MsReadyFailure>,
) -> Result<RuleRewriteRun, MsReadyFailure> {
    normalize_group_rules_with_finder_and_observer(
        topology,
        params,
        budget,
        find_candidates_limited,
        observe,
    )
}

fn normalize_group_rules_with_finder(
    topology: &TopologyBlock,
    params: &MsReadyParams,
    budget: &mut TransformationBudget,
    mut find: impl FnMut(&TopologyBlock, RuleId, usize) -> Result<Vec<RuleCandidate>, MsReadyFailure>,
) -> Result<RuleRewriteRun, MsReadyFailure> {
    normalize_group_rules_with_finder_and_observer(topology, params, budget, &mut find, |_| Ok(()))
}

fn normalize_group_rules_with_finder_and_observer(
    topology: &TopologyBlock,
    params: &MsReadyParams,
    budget: &mut TransformationBudget,
    mut find: impl FnMut(&TopologyBlock, RuleId, usize) -> Result<Vec<RuleCandidate>, MsReadyFailure>,
    mut observe: impl FnMut(&RuleEdit) -> Result<(), MsReadyFailure>,
) -> Result<RuleRewriteRun, MsReadyFailure> {
    validate_rule_data()?;
    crate::validate::validate_syntax(topology, &MoleculeProperties::default(), params)?;
    let mut current = topology.clone();
    let mut input_to_current = identity_mapping(&current);
    let mut events = Vec::new();
    let mut visited = BTreeSet::new();
    visited.insert(canonical_graph_key(&current)?);

    loop {
        let edits_at_sweep_start = events.len();
        budget.begin_sweep();
        for rule in RuleId::ORDER {
            loop {
                let annotate = |mut error: MsReadyFailure| {
                    error.evidence.location.rule_id = Some(rule.id().to_owned());
                    error
                };
                let matches = find(&current, rule, budget.remaining_matches()).map_err(annotate)?;
                budget.record_candidates(matches.len()).map_err(annotate)?;
                if matches.is_empty() {
                    break;
                }
                let current_state = canonical_graph_key(&current)?;
                let mut valid_candidates = Vec::with_capacity(matches.len());
                for candidate in matches {
                    let (
                        next,
                        step_mapping,
                        formula_before,
                        formula_after,
                        charge_before,
                        charge_after,
                    ) = apply_candidate(&current, &candidate).map_err(|mut error| {
                        error.evidence.location.rule_id = Some(rule.id().to_owned());
                        error.evidence.location.atom_indices = candidate
                            .atom_rows
                            .iter()
                            .filter_map(|row| {
                                input_to_current.atoms.new_to_old[*row].map(|id| id.index())
                            })
                            .collect();
                        error.evidence.location.bond_indices = candidate
                            .bond_rows
                            .iter()
                            .filter_map(|row| {
                                input_to_current.bonds.new_to_old[*row].map(|id| id.index())
                            })
                            .collect();
                        error
                    })?;
                    let state = canonical_graph_key(&next)?;
                    // An exact no-op is a fixed point, not a cycle. Charged nitro
                    // role selection may only improve the canonical representation.
                    if state == current_state
                        || (matches!(rule, RuleId::NitroInitial | RuleId::NitroFinal)
                            && current.atoms[candidate.atom_rows[0]].formal_charge() == 1
                            && state >= current_state)
                    {
                        continue;
                    }
                    let canonical = state;
                    let tie_key = (
                        candidate.atom_rows.clone(),
                        candidate.bond_rows.clone(),
                        format!("{:?}", candidate.hydrogen_transfer),
                    );
                    valid_candidates.push(ValidCandidate {
                        canonical,
                        tie_key,
                        candidate,
                        topology: next,
                        step_mapping,
                        formula_before,
                        formula_after,
                        charge_before,
                        charge_after,
                    });
                }
                valid_candidates.sort_by(|left, right| {
                    left.canonical
                        .cmp(&right.canonical)
                        .then_with(|| left.tie_key.cmp(&right.tie_key))
                });
                let Some(selected) = valid_candidates.into_iter().next() else {
                    break;
                };
                if visited.contains(&selected.canonical) {
                    return Err(failure(
                        Kind::TransformationCycle,
                        format!("{} revisits an earlier distinct state", rule.id()),
                    ));
                }
                budget.commit_edit().map_err(|mut error| {
                    error.evidence.location.rule_id = Some(rule.id().to_owned());
                    error.evidence.location.atom_indices = selected
                        .candidate
                        .atom_rows
                        .iter()
                        .filter_map(|row| {
                            input_to_current.atoms.new_to_old[*row].map(|id| id.index())
                        })
                        .collect();
                    error.evidence.location.bond_indices = selected
                        .candidate
                        .bond_rows
                        .iter()
                        .filter_map(|row| {
                            input_to_current.bonds.new_to_old[*row].map(|id| id.index())
                        })
                        .collect();
                    error
                })?;
                visited.insert(selected.canonical.clone());
                input_to_current = compose_mappings(
                    &input_to_current,
                    &selected.step_mapping,
                    topology,
                    &selected.topology,
                )?;
                let matched_atom_rows =
                    map_candidate_rows(&input_to_current, &selected.candidate.atom_rows)?;
                let matched_bond_rows =
                    map_candidate_bond_rows(&input_to_current, &selected.candidate.bond_rows)?;
                let appearance_order = RuleId::ORDER
                    .iter()
                    .position(|appearance| *appearance == rule)
                    .expect("rule belongs to the frozen appearance sequence")
                    + 1;
                events.push(RuleEdit {
                    rule_id: rule.id(),
                    appearance_order,
                    matched_atom_rows,
                    matched_bond_rows,
                    hydrogen_transfer: selected.candidate.hydrogen_transfer,
                    formula_before: selected.formula_before,
                    formula_after: selected.formula_after,
                    total_charge_before: selected.charge_before,
                    total_charge_after: selected.charge_after,
                    topology_mapping: input_to_current.clone(),
                });
                observe(events.last().expect("committed rule event"))?;
                current = selected.topology;
            }
        }
        if events.len() == edits_at_sweep_start {
            break;
        }
    }

    crate::validate::validate_input(&current, &MoleculeProperties::default(), params)?;
    Ok(RuleRewriteRun {
        topology: current,
        topology_mapping: input_to_current,
        edits: events,
    })
}

/// A bounded validation-only copy for the already-declared neutral nitro
/// repair ordering. Reuse the adopted matcher and checked edit, never relax
/// unrelated valence or return/log these transient repairs as analysis edits.
pub(crate) fn nitro_validation_graph(
    topology: &TopologyBlock,
) -> Result<TopologyBlock, MsReadyFailure> {
    let mut graph = topology.clone();
    loop {
        let candidates = find_candidates_limited(&graph, RuleId::NitroInitial, usize::MAX)?;
        let mut valid = Vec::new();
        for candidate in candidates {
            if graph.atoms[candidate.atom_rows[0]].formal_charge() != 0 {
                continue;
            }
            let (next, _, _, _, _, _) = apply_candidate(&graph, &candidate)?;
            valid.push((canonical_graph_key(&next)?, next));
        }
        valid.sort_by(|left, right| left.0.cmp(&right.0));
        let Some((_, next)) = valid.into_iter().next() else {
            break;
        };
        graph = next;
    }
    Ok(graph)
}

fn apply_candidate(
    topology: &TopologyBlock,
    candidate: &RuleCandidate,
) -> Result<(TopologyBlock, TopologyMapping, String, String, i32, i32), MsReadyFailure> {
    let before_valence = assign_for_match(topology)?;
    let before_formula = formula_signature(topology, &before_valence)?;
    let before_charge = total_charge(topology);
    let before_hydrogens = topology
        .atoms
        .iter()
        .enumerate()
        .map(|(row, _)| hydrogen_count(topology, &before_valence, row))
        .collect::<Result<Vec<_>, _>>()?;
    let before_isotopes = isotope_inventory(topology);
    let before_radicals = topology
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons())
        .collect::<Vec<_>>();

    let mut edit = topology
        .begin_batch_edit()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let mut touched_bonds = BTreeSet::new();
    for change in &candidate.bond_changes {
        let bond_id = BondId::new(change.row);
        if !touched_bonds.insert(change.row) {
            return Err(failure(
                Kind::InvariantViolation,
                format!(
                    "{} edits bond row {} more than once",
                    candidate.rule.id(),
                    change.row
                ),
            ));
        }
        let bond = edit
            .bond_mut(bond_id)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        bond.set_order(change.order);
        if let Some((begin, end)) = change.endpoints {
            bond.set_endpoints(AtomId::new(begin), AtomId::new(end));
        }
        clear_edited_bond_stereo(bond)?;
        bond.clear_computed_props()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    }
    let mut touched_atoms = BTreeSet::new();
    for change in &candidate.charge_changes {
        let atom = edit
            .atom_mut(AtomId::new(change.row))
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        atom.set_formal_charge(change.charge);
        touched_atoms.insert(change.row);
    }
    if let Some(transfer) = candidate.hydrogen_transfer {
        apply_hydrogen_transfer(&mut edit, transfer, &mut touched_atoms, &mut touched_bonds)?;
    }
    for row in touched_atoms {
        edit.atom_mut(AtomId::new(row))
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
            .clear_computed_props()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    }
    let (changed, mapping) = edit
        .finish()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    mapping
        .validate_for_counts(
            topology.atoms.len(),
            changed.atoms.len(),
            topology.bonds.len(),
            changed.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    if topology.atoms.len() != changed.atoms.len() || topology.bonds.len() != changed.bonds.len() {
        return Err(failure(
            Kind::InvariantViolation,
            "a directed normalization rule added or removed an atom or bond row",
        ));
    }
    validate_candidate_state(
        topology,
        &changed,
        candidate,
        &mapping,
        &before_isotopes,
        &before_radicals,
    )?;

    let after_valence = assign_for_match(&changed)?;
    let after_formula = formula_signature(&changed, &after_valence)?;
    let after_hydrogens = changed
        .atoms
        .iter()
        .enumerate()
        .map(|(row, _)| hydrogen_count(&changed, &after_valence, row))
        .collect::<Result<Vec<_>, _>>()?;
    let after_charge = total_charge(&changed);
    if before_formula != after_formula || before_charge != after_charge {
        return Err(failure(
            Kind::InvariantViolation,
            format!(
                "{} changed formula/charge: {before_formula} ({before_charge}) -> {after_formula} ({after_charge})",
                candidate.rule.id()
            ),
        ));
    }
    validate_hydrogen_delta(
        &before_hydrogens,
        &after_hydrogens,
        candidate.hydrogen_transfer,
    )?;
    Ok((
        changed,
        mapping,
        before_formula,
        after_formula,
        before_charge,
        after_charge,
    ))
}

fn validate_candidate_state(
    before: &TopologyBlock,
    after: &TopologyBlock,
    candidate: &RuleCandidate,
    mapping: &TopologyMapping,
    before_isotopes: &(Vec<(u8, Option<u16>)>, Vec<u16>),
    before_radicals: &[u8],
) -> Result<(), MsReadyFailure> {
    if &isotope_inventory(after) != before_isotopes
        || after
            .atoms
            .iter()
            .map(|atom| atom.radical_electrons())
            .collect::<Vec<_>>()
            != before_radicals
    {
        return Err(failure(
            Kind::InvariantViolation,
            format!("{} changed isotope or radical state", candidate.rule.id()),
        ));
    }
    let mut touched_atoms = BTreeSet::new();
    for change in &candidate.bond_changes {
        let bond = &before.bonds[change.row];
        touched_atoms.extend([bond.begin().index(), bond.end().index()]);
        if let Some((begin, end)) = change.endpoints {
            touched_atoms.extend([begin, end]);
        }
    }
    for change in &candidate.charge_changes {
        touched_atoms.insert(change.row);
    }
    if let Some(transfer) = candidate.hydrogen_transfer {
        touched_atoms.extend([transfer.donor, transfer.acceptor]);
        if let HydrogenCarrier::ExplicitAtom { atom, .. } = transfer.carrier {
            touched_atoms.insert(atom);
        }
    }
    let mut charge_changes = vec![None; before.atoms.len()];
    for change in &candidate.charge_changes {
        charge_changes[change.row] = Some(change.charge);
    }
    let mut expected_tracked: Vec<_> = before
        .atoms
        .iter()
        .map(|atom| atom.tracked_isotopic_hydrogens().to_vec())
        .collect();
    if let Some(HydrogenTransfer {
        donor,
        acceptor,
        carrier: HydrogenCarrier::TrackedIsotope { index, isotope },
    }) = candidate.hydrogen_transfer
    {
        let value = expected_tracked[donor].get(index).copied().ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                "tracked isotope H slot disappeared",
            )
        })?;
        if value != isotope {
            return Err(failure(
                Kind::InvariantViolation,
                "tracked isotope H slot changed before rule application",
            ));
        }
        expected_tracked[donor].remove(index);
        expected_tracked[acceptor].push(isotope);
    }
    let mut explicit_h_delta = vec![0i16; before.atoms.len()];
    if let Some(HydrogenTransfer {
        donor,
        acceptor,
        carrier: HydrogenCarrier::Bracket,
    }) = candidate.hydrogen_transfer
    {
        explicit_h_delta[donor] = -1;
        explicit_h_delta[acceptor] = 1;
    }
    for (old_row, old_atom) in before.atoms.iter().enumerate() {
        let new_id = mapping.atoms.old_to_new[old_row].ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                format!("rule removed atom row {old_row}"),
            )
        })?;
        let new_atom = &after.atoms[new_id.index()];
        let mut expected = old_atom.clone();
        if let Some(charge) = charge_changes[old_row] {
            expected.set_formal_charge(charge);
        }
        if explicit_h_delta[old_row] != 0 {
            let adjusted = i16::from(expected.explicit_hydrogens()) + explicit_h_delta[old_row];
            expected.set_explicit_hydrogens(u8::try_from(adjusted).map_err(|_| {
                failure(
                    Kind::InvariantViolation,
                    "bracket H transfer exceeds atom count",
                )
            })?);
        }
        if old_row < expected_tracked.len()
            && expected.tracked_isotopic_hydrogens() != expected_tracked[old_row]
        {
            expected.set_tracked_isotopic_hydrogens(expected_tracked[old_row].clone());
        }
        if touched_atoms.contains(&old_row) {
            expected
                .clear_computed_props()
                .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
        }
        if expected != *new_atom {
            return Err(failure(
                Kind::InvariantViolation,
                format!(
                    "{} changed protected atom state at original row {old_row}",
                    candidate.rule.id()
                ),
            ));
        }
    }
    for (old_row, old_bond) in before.bonds.iter().enumerate() {
        let new_id = mapping.bonds.old_to_new[old_row].ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                format!("rule removed bond row {old_row}"),
            )
        })?;
        let mut expected = old_bond.clone();
        if let Some(change) = candidate
            .bond_changes
            .iter()
            .find(|change| change.row == old_row)
        {
            expected.set_order(change.order);
            if let Some((begin, end)) = change.endpoints {
                expected.set_endpoints(AtomId::new(begin), AtomId::new(end));
            }
            clear_edited_bond_stereo(&mut expected)?;
        }
        if let Some(HydrogenTransfer {
            carrier: HydrogenCarrier::ExplicitAtom { bond, atom, .. },
            donor,
            acceptor,
        }) = candidate.hydrogen_transfer
            && bond == old_row
        {
            let new_endpoint =
                if old_bond.begin().index() == donor || old_bond.end().index() == donor {
                    acceptor
                } else {
                    return Err(failure(
                        Kind::InvariantViolation,
                        "explicit hydrogen bond is not attached to the declared donor",
                    ));
                };
            let (begin, end) = if old_bond.begin().index() == donor {
                (new_endpoint, atom)
            } else {
                (atom, new_endpoint)
            };
            expected.set_endpoints(AtomId::new(begin), AtomId::new(end));
            clear_edited_bond_stereo(&mut expected)?;
        }
        if expected != after.bonds[new_id.index()] {
            return Err(failure(
                Kind::InvariantViolation,
                format!(
                    "{} changed protected bond state at original row {old_row}",
                    candidate.rule.id()
                ),
            ));
        }
    }
    Ok(())
}

fn apply_hydrogen_transfer(
    edit: &mut cosmolkit_model::TopologyBatchEdit,
    transfer: HydrogenTransfer,
    touched_atoms: &mut BTreeSet<usize>,
    touched_bonds: &mut BTreeSet<usize>,
) -> Result<(), MsReadyFailure> {
    let donor = AtomId::new(transfer.donor);
    let acceptor = AtomId::new(transfer.acceptor);
    match transfer.carrier {
        HydrogenCarrier::Implicit => {}
        HydrogenCarrier::Bracket => {
            let donor_atom = edit
                .atom_mut(donor)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let remaining = donor_atom
                .explicit_hydrogens()
                .checked_sub(1)
                .ok_or_else(|| {
                    failure(
                        Kind::InvariantViolation,
                        "bracket H donor has no explicit H slot",
                    )
                })?;
            donor_atom.set_explicit_hydrogens(remaining);
            let acceptor_atom = edit
                .atom_mut(acceptor)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let added = acceptor_atom
                .explicit_hydrogens()
                .checked_add(1)
                .ok_or_else(|| {
                    failure(
                        Kind::InvariantViolation,
                        "bracket H acceptor count overflow",
                    )
                })?;
            acceptor_atom.set_explicit_hydrogens(added);
            touched_atoms.extend([transfer.donor, transfer.acceptor]);
        }
        HydrogenCarrier::ExplicitAtom { atom, bond, .. } => {
            if !touched_bonds.insert(bond) {
                return Err(failure(
                    Kind::InvariantViolation,
                    "explicit H transfer conflicts with a second edit to the same bond",
                ));
            }
            let bond = edit
                .bond_mut(BondId::new(bond))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let h_id = AtomId::new(atom);
            let (begin, end) = if bond.begin() == donor && bond.end() == h_id {
                (acceptor, h_id)
            } else if bond.end() == donor && bond.begin() == h_id {
                (h_id, acceptor)
            } else {
                return Err(failure(
                    Kind::InvariantViolation,
                    "explicit H bond no longer joins donor and selected H atom",
                ));
            };
            bond.set_endpoints(begin, end);
            clear_edited_bond_stereo(bond)?;
            bond.clear_computed_props()
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            touched_atoms.insert(transfer.donor);
            touched_atoms.insert(transfer.acceptor);
        }
        HydrogenCarrier::TrackedIsotope { index, isotope } => {
            let donor_atom = edit
                .atom_mut(donor)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let mut donor_isotopes = donor_atom.tracked_isotopic_hydrogens().to_vec();
            if donor_isotopes.get(index).copied() != Some(isotope) {
                return Err(failure(
                    Kind::InvariantViolation,
                    "selected tracked isotope H does not match the source donor slot",
                ));
            }
            donor_isotopes.remove(index);
            donor_atom.set_tracked_isotopic_hydrogens(donor_isotopes);
            let acceptor_atom = edit
                .atom_mut(acceptor)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
            let mut acceptor_isotopes = acceptor_atom.tracked_isotopic_hydrogens().to_vec();
            acceptor_isotopes.push(isotope);
            acceptor_atom.set_tracked_isotopic_hydrogens(acceptor_isotopes);
            touched_atoms.extend([transfer.donor, transfer.acceptor]);
        }
    }
    Ok(())
}

fn clear_edited_bond_stereo(bond: &mut cosmolkit_model::Bond) -> Result<(), MsReadyFailure> {
    bond.set_stereo(BondStereo::None)
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    bond.set_stereo_atoms(None);
    bond.set_source_stereo_atom_references(Vec::new());
    bond.set_direction(BondDirection::None);
    bond.set_unknown_stereo(false);
    Ok(())
}

fn validate_hydrogen_delta(
    before: &[u32],
    after: &[u32],
    transfer: Option<HydrogenTransfer>,
) -> Result<(), MsReadyFailure> {
    if before.len() != after.len() {
        return Err(failure(
            Kind::InvariantViolation,
            "rule changed the atom-row count",
        ));
    }
    for (row, (old, new)) in before.iter().zip(after).enumerate() {
        let expected = match transfer {
            Some(transfer) if row == transfer.donor => old.checked_sub(1),
            Some(transfer) if row == transfer.acceptor => old.checked_add(1),
            _ => Some(*old),
        }
        .ok_or_else(|| failure(Kind::InvariantViolation, "rule H delta overflows"))?;
        if *new != expected {
            return Err(failure(
                Kind::InvariantViolation,
                format!("rule changed H count at row {row}: {old} -> {new}, expected {expected}"),
            ));
        }
    }
    Ok(())
}

pub(crate) fn isotope_inventory(topology: &TopologyBlock) -> (Vec<(u8, Option<u16>)>, Vec<u16>) {
    let mut atom_isotopes: Vec<_> = topology
        .atoms
        .iter()
        .map(|atom| (atom.atomic_number(), atom.isotope()))
        .collect();
    atom_isotopes.sort_unstable();
    let mut tracked_hydrogen_isotopes: Vec<_> = topology
        .atoms
        .iter()
        .flat_map(|atom| atom.tracked_isotopic_hydrogens().iter().copied())
        .collect();
    tracked_hydrogen_isotopes.sort_unstable();
    (atom_isotopes, tracked_hydrogen_isotopes)
}

pub(crate) fn formula_counts(
    topology: &TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
) -> Result<std::collections::BTreeMap<u8, u32>, MsReadyFailure> {
    let mut counts = std::collections::BTreeMap::<u8, u32>::new();
    for atom in &topology.atoms {
        if atom.atomic_number() != 1 {
            *counts.entry(atom.atomic_number()).or_default() += 1;
        }
    }
    let mut hydrogen_count_total = 0u32;
    for (row, atom) in topology.atoms.iter().enumerate() {
        if atom.atomic_number() != 1 {
            hydrogen_count_total = hydrogen_count_total
                .checked_add(hydrogen_count(topology, valence, row)?)
                .ok_or_else(|| failure(Kind::InvariantViolation, "formula H count overflow"))?;
        } else if topology
            .adjacency
            .neighbors_of(row)
            .iter()
            .all(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 1)
        {
            hydrogen_count_total = hydrogen_count_total
                .checked_add(1)
                .ok_or_else(|| failure(Kind::InvariantViolation, "formula H count overflow"))?;
        }
    }
    if hydrogen_count_total != 0 {
        counts.insert(1, hydrogen_count_total);
    }
    Ok(counts)
}

pub(crate) fn formula_signature(
    topology: &TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
) -> Result<String, MsReadyFailure> {
    Ok(formula_counts(topology, valence)?
        .into_iter()
        .map(|(atomic_number, count)| format!("{atomic_number}:{count}"))
        .collect::<Vec<_>>()
        .join(";"))
}

pub(crate) fn total_charge(topology: &TopologyBlock) -> i32 {
    topology
        .atoms
        .iter()
        .map(|atom| i32::from(atom.formal_charge()))
        .sum()
}

fn identity_mapping(topology: &TopologyBlock) -> TopologyMapping {
    TopologyMapping {
        atoms: AtomMapping {
            old_to_new: (0..topology.atoms.len())
                .map(|row| Some(AtomId::new(row)))
                .collect(),
            new_to_old: (0..topology.atoms.len())
                .map(|row| Some(AtomId::new(row)))
                .collect(),
        },
        bonds: BondMapping {
            old_to_new: (0..topology.bonds.len())
                .map(|row| Some(BondId::new(row)))
                .collect(),
            new_to_old: (0..topology.bonds.len())
                .map(|row| Some(BondId::new(row)))
                .collect(),
        },
    }
}

pub(crate) fn compose_mappings(
    original_to_current: &TopologyMapping,
    current_to_next: &TopologyMapping,
    original: &TopologyBlock,
    next: &TopologyBlock,
) -> Result<TopologyMapping, MsReadyFailure> {
    compose_mapping_counts(
        original_to_current,
        current_to_next,
        original.atoms.len(),
        original.bonds.len(),
        next.atoms.len(),
        next.bonds.len(),
    )
}

pub(crate) fn compose_mapping_counts(
    original_to_current: &TopologyMapping,
    current_to_next: &TopologyMapping,
    original_atoms: usize,
    original_bonds: usize,
    next_atoms: usize,
    next_bonds: usize,
) -> Result<TopologyMapping, MsReadyFailure> {
    original_to_current
        .validate_for_counts(
            original_atoms,
            current_to_next.atoms.old_to_new.len(),
            original_bonds,
            current_to_next.bonds.old_to_new.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    current_to_next
        .validate_for_counts(
            current_to_next.atoms.old_to_new.len(),
            next_atoms,
            current_to_next.bonds.old_to_new.len(),
            next_bonds,
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    let mut atoms_old_to_new = vec![None; original_atoms];
    for (old, current_id) in original_to_current.atoms.old_to_new.iter().enumerate() {
        if let Some(current_id) = current_id {
            atoms_old_to_new[old] = current_to_next.atoms.old_to_new[current_id.index()];
        }
    }
    let atoms_new_to_old = current_to_next
        .atoms
        .new_to_old
        .iter()
        .map(|current_id| {
            current_id
                .and_then(|current_id| original_to_current.atoms.new_to_old[current_id.index()])
        })
        .collect();
    let mut bonds_old_to_new = vec![None; original_bonds];
    for (old, current_id) in original_to_current.bonds.old_to_new.iter().enumerate() {
        if let Some(current_id) = current_id {
            bonds_old_to_new[old] = current_to_next.bonds.old_to_new[current_id.index()];
        }
    }
    let bonds_new_to_old = current_to_next
        .bonds
        .new_to_old
        .iter()
        .map(|current_id| {
            current_id
                .and_then(|current_id| original_to_current.bonds.new_to_old[current_id.index()])
        })
        .collect();
    let mapping = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atoms_old_to_new,
            new_to_old: atoms_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bonds_old_to_new,
            new_to_old: bonds_new_to_old,
        },
    };
    mapping
        .validate_for_counts(original_atoms, next_atoms, original_bonds, next_bonds)
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    Ok(mapping)
}

pub(crate) fn canonical_graph_key(topology: &TopologyBlock) -> Result<String, MsReadyFailure> {
    let mut graph = crate::hydrogens::materialize_tracked_isotopic_hydrogens(topology)?;
    for atom in &mut graph.atoms {
        atom.set_chiral_tag(ChiralTag::Unspecified);
        atom.set_chiral_permutation(None);
        atom.set_unknown_stereo(false);
        atom.set_mol_parity(None);
        atom.set_mol_inversion_flag(None);
        atom.set_atom_map(None);
    }
    for bond in &mut graph.bonds {
        clear_edited_bond_stereo(bond)?;
    }
    graph.stereo_groups.clear();
    let record = SmilesRecord {
        topology: graph,
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default(),
    };
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        do_kekule: false,
        canonical: true,
        clean_stereo: false,
        rooted_at_atom: None,
        all_bonds_explicit: false,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: true,
    };
    let notation = write_smiles_with_params(&record, &params).map_err(|error| {
        failure(
            Kind::OwnerError,
            format!("canonical rule-state key: {error}"),
        )
    })?;
    String::from_utf8(notation.as_bytes().to_vec()).map_err(|error| {
        failure(
            Kind::OwnerError,
            format!("canonical SMILES is not UTF-8: {error}"),
        )
    })
}

fn map_candidate_rows(
    original_to_current: &TopologyMapping,
    rows: &[usize],
) -> Result<Vec<usize>, MsReadyFailure> {
    rows.iter()
        .map(|row| {
            let id = original_to_current
                .atoms
                .new_to_old
                .get(*row)
                .copied()
                .flatten()
                .ok_or_else(|| {
                    failure(
                        Kind::InvariantViolation,
                        "candidate row has no original mapping",
                    )
                })?;
            Ok(id.index())
        })
        .collect()
}

fn map_candidate_bond_rows(
    original_to_current: &TopologyMapping,
    rows: &[usize],
) -> Result<Vec<usize>, MsReadyFailure> {
    rows.iter()
        .map(|row| {
            original_to_current
                .bonds
                .new_to_old
                .get(*row)
                .copied()
                .flatten()
                .map(|id| id.index())
                .ok_or_else(|| {
                    failure(
                        Kind::InvariantViolation,
                        "candidate bond row has no original mapping",
                    )
                })
        })
        .collect()
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

#[cfg(test)]
mod rule_tests {
    use super::*;
    use crate::rules::find_candidates;
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

    fn parse(smiles: &str) -> TopologyBlock {
        parse_smiles(smiles, &SmilesParseParams::default())
            .unwrap_or_else(|error| panic!("cannot parse {smiles}: {error}"))
            .topology
    }

    #[test]
    fn enol_rule_rewrites_once_and_composes_original_mapping() {
        let input = parse("CC(C)=C(O)C");
        let result = normalize_group_rules(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.edits.len(), 1);
        assert_eq!(result.edits[0].rule_id, "enol_to_carbonyl");
        assert_eq!(result.edits[0].appearance_order, 3);
        assert_eq!(
            result.edits[0].formula_before,
            result.edits[0].formula_after
        );
        assert_eq!(
            result.edits[0].total_charge_before,
            result.edits[0].total_charge_after
        );
        assert_eq!(
            result.topology_mapping.atoms.old_to_new.len(),
            input.atoms.len()
        );
        assert_eq!(
            result.topology_mapping.bonds.old_to_new.len(),
            input.bonds.len()
        );
        assert!(
            result
                .topology_mapping
                .atoms
                .old_to_new
                .iter()
                .all(Option::is_some)
        );
        assert!(
            result
                .topology_mapping
                .bonds
                .old_to_new
                .iter()
                .all(Option::is_some)
        );
        assert!(result.edits[0].hydrogen_transfer.is_some());
    }

    #[test]
    fn every_directed_rule_applies_without_losing_source_rows() {
        let examples = [
            ("CN(=O)=O", "nitro_initial"),
            ("CC(C)=C(NC)C", "enamine_to_imine"),
            ("CC(C)=C(O)C", "enol_to_carbonyl"),
            ("CC#CO", "ynol_to_ketene"),
            ("C[N-][N+]#N", "azide_mesomer"),
        ];
        for (smiles, rule_id) in examples {
            let input = parse(smiles);
            let result = normalize_group_rules(&input, &MsReadyParams::default())
                .unwrap_or_else(|error| panic!("{rule_id} rewrite failed: {error}"));
            assert!(result.edits.iter().any(|edit| edit.rule_id == rule_id));
            assert_eq!(result.topology.atoms.len(), input.atoms.len(), "{rule_id}");
            assert_eq!(result.topology.bonds.len(), input.bonds.len(), "{rule_id}");
            assert!(
                result
                    .topology_mapping
                    .atoms
                    .old_to_new
                    .iter()
                    .all(Option::is_some)
            );
            assert!(
                result
                    .topology_mapping
                    .bonds
                    .old_to_new
                    .iter()
                    .all(Option::is_some)
            );
            assert!(result.edits.iter().all(|edit| {
                edit.formula_before == edit.formula_after
                    && edit.total_charge_before == edit.total_charge_after
            }));
            assert!(
                normalize_group_rules(&result.topology, &MsReadyParams::default())
                    .unwrap()
                    .edits
                    .is_empty(),
                "second-run {rule_id}"
            );
        }
    }

    #[test]
    fn aromatic_tautomer_like_graphs_do_not_match_nonaromatic_rules() {
        let aromatic = parse("c1ccccc1");
        for rule in RuleId::ORDER {
            assert!(
                find_candidates(&aromatic, rule).unwrap().is_empty(),
                "{}",
                rule.id()
            );
        }
    }

    #[test]
    fn explicit_and_tracked_deuterium_and_tritium_follow_the_transferred_hydrogen() {
        for (isotope, smiles) in [(2, "CC(C)=C(O[2H])C(C)C"), (3, "CC(C)=C(O[3H])C(C)C")] {
            let input = parse(smiles);
            let result =
                normalize_group_rules(&input, &MsReadyParams::default()).unwrap_or_else(|error| {
                    panic!("explicit isotope {isotope} transfer failed: {error}")
                });
            assert_eq!(
                result.edits[0].formula_before,
                result.edits[0].formula_after
            );
            assert!(result.edits[0].hydrogen_transfer.is_some());
            assert_eq!(result.topology.atoms.len(), input.atoms.len());

            let mut tracked_input = parse("CC(C)=C(O)C(C)C");
            let oxygen = tracked_input
                .atoms
                .iter()
                .position(|atom| atom.element() == cosmolkit_types::Element::O)
                .unwrap();
            tracked_input.atoms[oxygen].set_no_implicit(true);
            tracked_input.atoms[oxygen].set_tracked_isotopic_hydrogens(vec![isotope]);
            let tracked = normalize_group_rules(&tracked_input, &MsReadyParams::default())
                .unwrap_or_else(|error| {
                    panic!("tracked isotope {isotope} transfer failed: {error}")
                });
            let recipient = tracked.edits[0].hydrogen_transfer.unwrap().acceptor;
            assert_eq!(
                tracked.topology.atoms[recipient].tracked_isotopic_hydrogens(),
                &[isotope]
            );
        }
    }

    #[test]
    fn isotope_asymmetric_nitro_role_selection_is_atom_order_independent() {
        let first = parse("C[N](=O)=[18O]");
        let reordered = parse("[18O]=[N](C)=O");
        let first_result = normalize_group_rules(&first, &MsReadyParams::default()).unwrap();
        let reordered_result =
            normalize_group_rules(&reordered, &MsReadyParams::default()).unwrap();
        assert_eq!(
            canonical_graph_key(&first_result.topology).unwrap(),
            canonical_graph_key(&reordered_result.topology).unwrap()
        );
        assert_eq!(first_result.edits.len(), 1);
        assert_eq!(reordered_result.edits.len(), 1);
    }

    #[test]
    fn match_limit_and_multiple_nitro_edits_are_bounded_before_strict_validation() {
        let one_nitro = parse("CN(=O)=O");
        let limited = MsReadyParams {
            max_matches_per_sweep: 1,
            ..MsReadyParams::default()
        };
        let error = normalize_group_rules(&one_nitro, &limited).unwrap_err();
        assert_eq!(error.kind, Kind::TransformationLimit);

        let two_nitros = parse("CN(=O)=O.CN(=O)=O");
        let result = normalize_group_rules(&two_nitros, &MsReadyParams::default()).unwrap();
        assert_eq!(result.edits.len(), 2);
        assert!(
            result
                .edits
                .iter()
                .all(|edit| edit.rule_id == "nitro_initial")
        );
    }

    #[test]
    fn rule_limit_is_typed_and_does_not_return_a_truncated_graph() {
        let input = parse("CC(C)=C(O)C.CC(C)=C(O)C");
        let params = MsReadyParams {
            max_edits: 1,
            ..MsReadyParams::default()
        };
        let error = normalize_group_rules(&input, &params).unwrap_err();
        assert_eq!(error.kind, Kind::TransformationLimit);
        assert_eq!(error.stage, Stage::Normalization);
    }

    #[test]
    fn acetone_h_filled_slots_converge_across_hydrogen_depictions() {
        let expected = canonical_graph_key(&parse("CC(=O)C")).unwrap();
        for smiles in ["C=C(O)C", "[CH2]=C([OH])C", "[H]C([H])=C(O[H])C"] {
            let input = parse(smiles);
            let result = normalize_group_rules(&input, &MsReadyParams::default()).unwrap();
            assert_eq!(result.edits.len(), 1, "{smiles}");
            assert_eq!(
                crate::identify::standard_identity(&result.topology)
                    .unwrap()
                    .full_key,
                "CSCPPACGZOOCGX-UHFFFAOYSA-N",
                "real acetone identity: {smiles}"
            );
            let folded = crate::hydrogens::fold_ordinary_explicit_hydrogens(
                &input,
                &result.topology,
                &result.topology_mapping,
            )
            .unwrap();
            assert_eq!(
                canonical_graph_key(&folded.topology).unwrap(),
                expected,
                "{smiles}"
            );
            assert!(
                normalize_group_rules(&result.topology, &MsReadyParams::default())
                    .unwrap()
                    .edits
                    .is_empty()
            );
        }
        let primary = normalize_group_rules(&parse("C=CN"), &MsReadyParams::default()).unwrap();
        assert_eq!(
            primary.edits.len(),
            1,
            "remaining N-bound H fills its R slot"
        );
        assert_eq!(
            canonical_graph_key(&primary.topology).unwrap(),
            canonical_graph_key(&parse("CC=N")).unwrap()
        );
        let enamine = normalize_group_rules(&parse("C=CNC"), &MsReadyParams::default()).unwrap();
        assert_eq!(enamine.edits.len(), 1);
        assert_eq!(
            canonical_graph_key(&enamine.topology).unwrap(),
            canonical_graph_key(&parse("CC=NC")).unwrap()
        );
    }

    #[test]
    fn charged_isotope_asymmetric_nitro_converges_and_is_idempotent() {
        let expected =
            normalize_group_rules(&parse("CN(=O)=[18O]"), &MsReadyParams::default()).unwrap();
        let key = canonical_graph_key(&expected.topology).unwrap();
        for smiles in ["C[N+](=O)[18O-]", "C[N+]([O-])=[18O]", "[18O-][N+](C)=O"] {
            let result = normalize_group_rules(&parse(smiles), &MsReadyParams::default()).unwrap();
            assert_eq!(
                canonical_graph_key(&result.topology).unwrap(),
                key,
                "{smiles}"
            );
            assert!(
                normalize_group_rules(&result.topology, &MsReadyParams::default())
                    .unwrap()
                    .edits
                    .is_empty()
            );
        }
    }

    #[test]
    fn distal_ynol_radical_control_remains_valid_and_unedited() {
        let mut input = parse("C#CO");
        input.atoms[0].set_radical_electrons(1);
        crate::validate::strict_owner_valence(&input).unwrap();
        assert!(
            find_candidates(&input, RuleId::YnolToKetene)
                .unwrap()
                .is_empty()
        );
        let result = normalize_group_rules(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.topology, input);
        assert!(result.edits.is_empty());
    }

    #[test]
    fn overlapping_rules_obey_order_and_rematch_adjacent_sites_idempotently() {
        let overlap = parse("C=C(O)NC");
        assert!(
            !find_candidates(&overlap, RuleId::EnamineToImine)
                .unwrap()
                .is_empty()
        );
        assert!(
            !find_candidates(&overlap, RuleId::EnolToCarbonyl)
                .unwrap()
                .is_empty()
        );
        let result = normalize_group_rules(&overlap, &MsReadyParams::default()).unwrap();
        assert_eq!(result.edits.len(), 1);
        assert_eq!(result.edits[0].rule_id, "enamine_to_imine");
        assert_eq!(
            canonical_graph_key(&result.topology).unwrap(),
            canonical_graph_key(&parse("CC(O)=NC")).unwrap()
        );
        let first =
            normalize_group_rules(&parse("C=C(O)C=C(O)C"), &MsReadyParams::default()).unwrap();
        let permuted =
            normalize_group_rules(&parse("CC(O)=CC(O)=C"), &MsReadyParams::default()).unwrap();
        assert_eq!(first.edits.len(), 2);
        assert_eq!(
            canonical_graph_key(&first.topology).unwrap(),
            canonical_graph_key(&permuted.topology).unwrap()
        );
        for topology in [&result.topology, &first.topology] {
            assert!(
                normalize_group_rules(topology, &MsReadyParams::default())
                    .unwrap()
                    .edits
                    .is_empty()
            );
        }
    }

    #[test]
    fn real_rewrite_boundary_rejects_a_distinct_state_cycle() {
        let input = parse("C=C(O)C");
        let forward = find_candidates(&input, RuleId::EnolToCarbonyl)
            .unwrap()
            .remove(0);
        let mut reverse = forward.clone();
        for change in &mut reverse.bond_changes {
            change.order = input.bonds[change.row].order();
        }
        let transfer = reverse.hydrogen_transfer.as_mut().unwrap();
        std::mem::swap(&mut transfer.donor, &mut transfer.acceptor);
        let params = MsReadyParams::default();
        let mut budget = TransformationBudget::new(&params);
        let initial = canonical_graph_key(&input).unwrap();
        let mut attempted_edits = 0;
        // Inject a reversible test-only catalogue through the same private
        // boundary used by production. Both candidates perform real, checked
        // graph/H edits; there is no public testing API or direct set assertion.
        let error =
            normalize_group_rules_with_finder(&input, &params, &mut budget, |topology, rule, _| {
                if rule != RuleId::EnolToCarbonyl {
                    return Ok(Vec::new());
                }
                attempted_edits += 1;
                Ok(vec![if canonical_graph_key(topology)? == initial {
                    forward.clone()
                } else {
                    reverse.clone()
                }])
            })
            .unwrap_err();
        assert_eq!(attempted_edits, 2);
        assert_eq!(error.kind, Kind::TransformationCycle);
        assert_eq!(error.stage, Stage::Normalization);
        assert!(error.context.contains("enol_to_carbonyl"));
        assert_eq!(budget.edits, 1, "repeated state is rejected before commit");
    }

    #[test]
    fn canonical_state_set_mechanics_are_repeatable() {
        let topology = parse("CC");
        let key = canonical_graph_key(&topology).unwrap();
        let mut visited = BTreeSet::new();
        assert!(visited.insert(key.clone()));
        assert!(!visited.insert(key));
    }
}
