use std::collections::{BTreeSet, VecDeque};

use cosmolkit_core::{
    ValenceAssignment, ValenceModel, assign_valence_for_topology,
    assign_valence_with_options_for_topology,
};
use cosmolkit_model::{AtomId, MoleculeProperties, TopologyBlock, TopologyMapping};
use cosmolkit_types::{BondOrder, Element};

use crate::hydrogens::{
    HydrogenDisposition, HydrogenDispositionKind, ProtonCarrier, add_reservoir_proton,
    available_carriers, fold_ordinary_explicit_hydrogens, remove_reservoir_carrier,
};
use crate::rewrite::{
    TransformationBudget, canonical_graph_key, compose_mappings, formula_counts, formula_signature,
    isotope_inventory, total_charge,
};
use crate::{
    MsReadyFailure, MsReadyFailureKind as Kind, MsReadyParams, MsReadyStage as Stage,
    MsReadyWarning,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ProtonEventKind {
    InternalTransfer,
    ReservoirDeprotonation,
    ReservoirProtonation,
    LabelledReservoirBlocked,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ProtonCarrierRecord {
    InferredOrdinary,
    BracketOrdinary,
    ExplicitAtom {
        input_atom_row: usize,
        isotope: Option<u16>,
    },
    TrackedIsotope {
        isotope: u16,
    },
    Reservoir,
    None,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ProtonEvent {
    pub kind: ProtonEventKind,
    pub donor_input_row: Option<usize>,
    pub acceptor_input_row: Option<usize>,
    pub carrier: ProtonCarrierRecord,
    pub delta_hydrogens: i32,
    pub delta_charge: i32,
    pub formula_before: String,
    pub formula_after: String,
    pub total_charge_before: i32,
    pub total_charge_after: i32,
    /// Original-input-to-current mapping immediately after this transition.
    pub topology_mapping: TopologyMapping,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ChargeNormalizationRun {
    pub topology: TopologyBlock,
    pub topology_mapping: TopologyMapping,
    pub events: Vec<ProtonEvent>,
    pub warnings: Vec<MsReadyWarning>,
    pub hydrogen_dispositions: Vec<HydrogenDisposition>,
}

#[derive(Clone)]
struct ValidTransition {
    canonical: String,
    tie_key: (usize, usize, (u8, usize, u16)),
    topology: TopologyBlock,
    mapping: TopologyMapping,
    event: ProtonEvent,
    removed_input_hydrogen: Option<usize>,
}

/// Apply the frozen N+/O- proton classes, then fold ordinary explicit H via core.
pub(crate) fn normalize_proton_charges(
    input: &TopologyBlock,
    params: &MsReadyParams,
) -> Result<ChargeNormalizationRun, MsReadyFailure> {
    let mut budget = TransformationBudget::new(params);
    normalize_proton_charges_with_budget(input, params, &mut budget)
}

pub(crate) fn normalize_proton_charges_with_budget(
    input: &TopologyBlock,
    params: &MsReadyParams,
    budget: &mut TransformationBudget,
) -> Result<ChargeNormalizationRun, MsReadyFailure> {
    normalize_proton_charges_observed(input, params, budget, |_| Ok(()))
}

pub(crate) fn normalize_proton_charges_observed(
    input: &TopologyBlock,
    params: &MsReadyParams,
    budget: &mut TransformationBudget,
    mut observe: impl FnMut(&ProtonEvent) -> Result<(), MsReadyFailure>,
) -> Result<ChargeNormalizationRun, MsReadyFailure> {
    crate::validate::validate_syntax(input, &MoleculeProperties::default(), params)?;
    budget.begin_sweep();
    strict_valence(input)?;
    let mut topology = input.clone();
    let mut original_to_current = TopologyMapping::identity(input.atoms.len(), input.bonds.len());
    let mut events = Vec::new();
    let mut warnings = radical_warnings(input);
    let mut reservoir_deleted = BTreeSet::new();
    let components = component_ids(&topology);

    // BEGIN MSREADY REFERENCE §7
    // MSReady❗❗: Engineering recommendation. Charge normalization requires a local rule, an allowed hydrogen change and a valid product valence. Keep both input and product formula/charge so you can verify what changed.
    // This frozen Native policy prioritizes balanced internal transfer; the reference explicitly does not prove the historical Indigo control flow or performance.
    // END MSREADY REFERENCE
    loop {
        let candidates =
            internal_transitions(input, &topology, &original_to_current, &components, budget)?;
        let Some(selected) = select_transition(candidates) else {
            break;
        };
        budget.commit_edit().map_err(|mut error| {
            error.evidence.location.rule_id = Some(format!("{:?}", selected.event.kind));
            error.evidence.location.atom_indices = selected
                .event
                .donor_input_row
                .into_iter()
                .chain(selected.event.acceptor_input_row)
                .collect();
            error
        })?;
        if let Some(hydrogen_row) = selected.removed_input_hydrogen {
            reservoir_deleted.insert(hydrogen_row);
        }
        topology = selected.topology;
        original_to_current = selected.mapping;
        events.push(selected.event);
        observe(events.last().expect("committed proton event"))?;
    }

    // Unpaired N+ sites may exchange only ordinary H with the reservoir.
    loop {
        let candidates =
            reservoir_deprotonations(input, &topology, &original_to_current, &components, budget)?;
        let Some(selected) = select_transition(candidates) else {
            break;
        };
        budget.commit_edit().map_err(|mut error| {
            error.evidence.location.rule_id = Some(format!("{:?}", selected.event.kind));
            error.evidence.location.atom_indices = selected
                .event
                .donor_input_row
                .into_iter()
                .chain(selected.event.acceptor_input_row)
                .collect();
            error
        })?;
        if let Some(hydrogen_row) = selected.removed_input_hydrogen {
            reservoir_deleted.insert(hydrogen_row);
        }
        topology = selected.topology;
        original_to_current = selected.mapping;
        events.push(selected.event);
        observe(events.last().expect("committed proton event"))?;
    }
    let before_warning_events = events.len();
    append_labelled_reservoir_warnings(
        &topology,
        &original_to_current,
        &mut events,
        &mut warnings,
    )?;

    for event in &events[before_warning_events..] {
        observe(event)?;
    }

    // Unpaired eligible O- sites receive one ordinary reservoir proton.
    loop {
        let candidates =
            reservoir_protonations(input, &topology, &original_to_current, &components, budget)?;
        let Some(selected) = select_transition(candidates) else {
            break;
        };
        budget.commit_edit().map_err(|mut error| {
            error.evidence.location.rule_id = Some(format!("{:?}", selected.event.kind));
            error.evidence.location.atom_indices = selected
                .event
                .donor_input_row
                .into_iter()
                .chain(selected.event.acceptor_input_row)
                .collect();
            error
        })?;
        topology = selected.topology;
        original_to_current = selected.mapping;
        events.push(selected.event);
        observe(events.last().expect("committed proton event"))?;
    }

    let folded = fold_ordinary_explicit_hydrogens(input, &topology, &original_to_current)?;
    let mut hydrogen_dispositions = folded.dispositions;
    hydrogen_dispositions.extend(reservoir_deleted.into_iter().map(|input_atom_row| {
        HydrogenDisposition {
            input_atom_row,
            kind: HydrogenDispositionKind::DeletedToReservoir,
        }
    }));
    hydrogen_dispositions.sort_by_key(|disposition| disposition.input_atom_row);
    let expected_hydrogen_rows = input
        .atoms
        .iter()
        .filter(|atom| atom.atomic_number() == 1)
        .count();
    if hydrogen_dispositions.len() != expected_hydrogen_rows
        || hydrogen_dispositions
            .windows(2)
            .any(|pair| pair[0].input_atom_row == pair[1].input_atom_row)
    {
        return Err(failure(
            Kind::InvariantViolation,
            "explicit-H disposition ledger is incomplete or duplicates an input row",
        ));
    }
    let final_valence = strict_valence(&folded.topology)?;
    let final_formula = formula_counts(&folded.topology, &final_valence)?;
    let pre_fold_valence = strict_valence(&topology)?;
    let pre_fold_formula = formula_counts(&topology, &pre_fold_valence)?;
    if final_formula != pre_fold_formula
        || isotope_signature(&folded.topology) != isotope_signature(&topology)
        || total_charge(&folded.topology) != total_charge(&topology)
    {
        return Err(failure(
            Kind::InvariantViolation,
            "folding ordinary explicit H changed formula, isotope composition, or charge",
        ));
    }
    folded
        .topology_mapping
        .validate_for_counts(
            input.atoms.len(),
            folded.topology.atoms.len(),
            input.bonds.len(),
            folded.topology.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    Ok(ChargeNormalizationRun {
        topology: folded.topology,
        topology_mapping: folded.topology_mapping,
        events,
        warnings,
        hydrogen_dispositions,
    })
}

fn internal_transitions(
    original: &TopologyBlock,
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
    components: &[usize],
    budget: &mut TransformationBudget,
) -> Result<Vec<ValidTransition>, MsReadyFailure> {
    // BEGIN MSREADY REFERENCE §7
    // MSReady❗❗: A prototropic tautomer or internal proton transfer ordinarily preserves elemental/isotopic formula and total charge.
    // END MSREADY REFERENCE
    let valence = flexible_valence(topology)?;
    let donors = eligible_n_donors(topology, &valence)?;
    let acceptors = eligible_o_acceptors(topology, &valence)?;
    let mut transitions = Vec::new();
    for (donor, carriers) in donors {
        for (acceptor, _) in &acceptors {
            if components[donor] != components[*acceptor] {
                continue;
            }
            for carrier in carriers.iter().copied() {
                budget.record_candidates(1)?;
                match apply_transition(
                    original,
                    topology,
                    mapping,
                    ProtonEventKind::InternalTransfer,
                    Some(donor),
                    Some(*acceptor),
                    carrier,
                ) {
                    Ok(transition) => transitions.push(transition),
                    Err(error) if error.kind == Kind::InvalidValence => continue,
                    Err(error) => return Err(error),
                }
            }
        }
    }
    Ok(transitions)
}

fn reservoir_deprotonations(
    original: &TopologyBlock,
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
    _components: &[usize],
    budget: &mut TransformationBudget,
) -> Result<Vec<ValidTransition>, MsReadyFailure> {
    let valence = flexible_valence(topology)?;
    let donors = eligible_n_donors(topology, &valence)?;
    let mut transitions = Vec::new();
    for (donor, carriers) in donors {
        let Some(carrier) = carriers
            .iter()
            .copied()
            .find(|carrier| !carrier.is_labelled())
        else {
            continue;
        };
        budget.record_candidates(1)?;
        match apply_transition(
            original,
            topology,
            mapping,
            ProtonEventKind::ReservoirDeprotonation,
            Some(donor),
            None,
            carrier,
        ) {
            Ok(transition) => transitions.push(transition),
            Err(error) if error.kind == Kind::InvalidValence => continue,
            Err(error) => return Err(error),
        }
    }
    Ok(transitions)
}

fn reservoir_protonations(
    original: &TopologyBlock,
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
    _components: &[usize],
    budget: &mut TransformationBudget,
) -> Result<Vec<ValidTransition>, MsReadyFailure> {
    let valence = flexible_valence(topology)?;
    let acceptors = eligible_o_acceptors(topology, &valence)?;
    let mut transitions = Vec::new();
    for (acceptor, _) in acceptors {
        budget.record_candidates(1)?;
        match apply_transition(
            original,
            topology,
            mapping,
            ProtonEventKind::ReservoirProtonation,
            None,
            Some(acceptor),
            ProtonCarrier::Implicit,
        ) {
            Ok(transition) => transitions.push(transition),
            Err(error) if error.kind == Kind::InvalidValence => continue,
            Err(error) => return Err(error),
        }
    }
    Ok(transitions)
}

fn apply_transition(
    original: &TopologyBlock,
    topology: &TopologyBlock,
    original_to_current: &TopologyMapping,
    kind: ProtonEventKind,
    donor: Option<usize>,
    acceptor: Option<usize>,
    carrier: ProtonCarrier,
) -> Result<ValidTransition, MsReadyFailure> {
    let before_valence = flexible_valence(topology)?;
    let formula_before = formula_signature(topology, &before_valence)?;
    let before_counts = formula_counts(topology, &before_valence)?;
    let charge_before = total_charge(topology);
    let isotope_before = isotope_signature(topology);
    let radicals_before: Vec<_> = topology
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons())
        .collect();
    let donor_input = donor
        .map(|row| input_row_for_current(original_to_current, row, true))
        .transpose()?;
    let acceptor_input = acceptor
        .map(|row| input_row_for_current(original_to_current, row, true))
        .transpose()?;
    let carrier_record = if kind == ProtonEventKind::ReservoirProtonation {
        ProtonCarrierRecord::Reservoir
    } else {
        match carrier {
            ProtonCarrier::Implicit => ProtonCarrierRecord::InferredOrdinary,
            ProtonCarrier::Bracket => ProtonCarrierRecord::BracketOrdinary,
            ProtonCarrier::ExplicitAtom { atom, isotope, .. } => {
                ProtonCarrierRecord::ExplicitAtom {
                    input_atom_row: input_row_for_current(original_to_current, atom, true)?,
                    isotope,
                }
            }
            ProtonCarrier::TrackedIsotope { isotope, .. } => {
                ProtonCarrierRecord::TrackedIsotope { isotope }
            }
        }
    };
    let mut edit = topology
        .begin_batch_edit()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let mut removed_current_hydrogen = None;
    match kind {
        ProtonEventKind::InternalTransfer => {
            let donor = donor.expect("internal transfer has N donor");
            let acceptor = acceptor.expect("internal transfer has O acceptor");
            crate::hydrogens::transfer_carrier(&mut edit, topology, donor, acceptor, carrier)?;
            edit.atom_mut(AtomId::new(donor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
                .set_formal_charge(0);
            edit.atom_mut(AtomId::new(acceptor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
                .set_formal_charge(0);
        }
        ProtonEventKind::ReservoirDeprotonation => {
            let donor = donor.expect("reservoir deprotonation has N donor");
            removed_current_hydrogen =
                remove_reservoir_carrier(&mut edit, topology, donor, carrier)?;
            edit.atom_mut(AtomId::new(donor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
                .set_formal_charge(0);
        }
        ProtonEventKind::ReservoirProtonation => {
            let acceptor = acceptor.expect("reservoir protonation has O acceptor");
            add_reservoir_proton(&mut edit, acceptor)?;
            edit.atom_mut(AtomId::new(acceptor))
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
                .set_formal_charge(0);
        }
        ProtonEventKind::LabelledReservoirBlocked => {
            return Err(failure(
                Kind::InvariantViolation,
                "blocked labelled-H event cannot edit a topology",
            ));
        }
    }
    let mut touched = BTreeSet::new();
    touched.extend(donor);
    touched.extend(acceptor);
    if let ProtonCarrier::ExplicitAtom { atom, .. } = carrier {
        touched.insert(atom);
    }
    for row in touched {
        edit.atom_mut(AtomId::new(row))
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?
            .clear_computed_props()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    }
    let (changed, step_mapping) = edit
        .finish()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    step_mapping
        .validate_for_counts(
            topology.atoms.len(),
            changed.atoms.len(),
            topology.bonds.len(),
            changed.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    let mapping = compose_mappings(original_to_current, &step_mapping, original, &changed)?;
    let after_valence = strict_valence(&changed)?;
    let formula_after = formula_signature(&changed, &after_valence)?;
    let mut expected_counts = before_counts;
    // BEGIN MSREADY REFERENCE §7
    // MSReady❗❗: Removing a proton reduces H count by one and charge by one. Adding a proton increases both by one.
    // END MSREADY REFERENCE
    let (delta_hydrogens, delta_charge): (i32, i32) = match kind {
        ProtonEventKind::InternalTransfer => (0, 0),
        ProtonEventKind::ReservoirDeprotonation => (-1, -1),
        ProtonEventKind::ReservoirProtonation => (1, 1),
        ProtonEventKind::LabelledReservoirBlocked => (0, 0),
    };
    if delta_hydrogens < 0 {
        let hydrogens = expected_counts.get_mut(&1).ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                "reservoir loss has no formula H count",
            )
        })?;
        *hydrogens = hydrogens
            .checked_sub(delta_hydrogens.unsigned_abs())
            .ok_or_else(|| failure(Kind::InvariantViolation, "reservoir H count underflow"))?;
        if *hydrogens == 0 {
            expected_counts.remove(&1);
        }
    } else if delta_hydrogens > 0 {
        *expected_counts.entry(1).or_default() += delta_hydrogens as u32;
    }
    let actual_counts = formula_counts(&changed, &after_valence)?;
    if actual_counts != expected_counts {
        return Err(failure(
            Kind::InvariantViolation,
            format!("proton event formula mismatch: {formula_before} -> {formula_after}"),
        ));
    }
    let charge_after = total_charge(&changed);
    if charge_after != charge_before + delta_charge {
        return Err(failure(
            Kind::InvariantViolation,
            format!("proton event charge mismatch: {charge_before} -> {charge_after}"),
        ));
    }
    if isotope_signature(&changed) != isotope_before {
        return Err(failure(
            Kind::InvariantViolation,
            "proton event changed a heavy-atom, D/T, or tracked-isotope inventory",
        ));
    }
    for (old_row, old_radical) in radicals_before.iter().copied().enumerate() {
        if let Some(new_id) = step_mapping.atoms.old_to_new[old_row] {
            if changed.atoms[new_id.index()].radical_electrons() != old_radical {
                return Err(failure(
                    Kind::InvariantViolation,
                    format!("proton event changed radical electrons at source row {old_row}"),
                ));
            }
        } else if Some(old_row) != removed_current_hydrogen
            || old_radical != 0
            || topology.atoms[old_row].atomic_number() != 1
        {
            return Err(failure(
                Kind::InvariantViolation,
                format!("proton event unexpectedly deleted atom row {old_row}"),
            ));
        }
    }
    let canonical = canonical_graph_key(&changed)?;
    let input_tie_row = donor_input.or(acceptor_input).unwrap_or(usize::MAX);
    let tie_key = (
        input_tie_row,
        acceptor_input.unwrap_or(usize::MAX),
        carrier.preference(),
    );
    let event = ProtonEvent {
        kind,
        donor_input_row: donor_input,
        acceptor_input_row: acceptor_input,
        carrier: carrier_record,
        delta_hydrogens,
        delta_charge,
        formula_before,
        formula_after,
        total_charge_before: charge_before,
        total_charge_after: charge_after,
        topology_mapping: mapping.clone(),
    };
    let removed_input_hydrogen = removed_current_hydrogen
        .map(|row| input_row_for_current(original_to_current, row, true))
        .transpose()?;
    Ok(ValidTransition {
        canonical,
        tie_key,
        topology: changed,
        mapping,
        event,
        removed_input_hydrogen,
    })
}

fn append_labelled_reservoir_warnings(
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
    events: &mut Vec<ProtonEvent>,
    warnings: &mut Vec<MsReadyWarning>,
) -> Result<(), MsReadyFailure> {
    let valence = flexible_valence(topology)?;
    let donors = eligible_n_donors(topology, &valence)?;
    let mut seen = BTreeSet::new();
    for (donor, carriers) in donors {
        if carriers.iter().any(|carrier| !carrier.is_labelled()) {
            continue;
        }
        let input_row = input_row_for_current(mapping, donor, true)?;
        if !seen.insert(input_row) {
            continue;
        }
        let isotopes: Vec<_> = carriers
            .iter()
            .filter_map(|carrier| match carrier {
                ProtonCarrier::ExplicitAtom {
                    isotope: Some(isotope),
                    ..
                }
                | ProtonCarrier::TrackedIsotope { isotope, .. } => Some(*isotope),
                _ => None,
            })
            .collect();
        if isotopes.is_empty() {
            continue;
        }
        let event = ProtonEvent {
            kind: ProtonEventKind::LabelledReservoirBlocked,
            donor_input_row: Some(input_row),
            acceptor_input_row: None,
            carrier: ProtonCarrierRecord::None,
            delta_hydrogens: 0,
            delta_charge: 0,
            formula_before: formula_signature(topology, &valence)?,
            formula_after: formula_signature(topology, &valence)?,
            total_charge_before: total_charge(topology),
            total_charge_after: total_charge(topology),
            topology_mapping: mapping.clone(),
        };
        events.push(event);
        warnings.push(MsReadyWarning {
            code: "labelled_hydrogen_not_removed_from_reservoir".to_owned(),
            message: format!(
                "retained N+ at input atom row {input_row}; only isotope-labelled H {isotopes:?} is available for deprotonation"
            ),
        });
    }
    Ok(())
}

fn select_transition(mut candidates: Vec<ValidTransition>) -> Option<ValidTransition> {
    candidates.sort_by(|left, right| {
        left.canonical
            .cmp(&right.canonical)
            .then_with(|| left.tie_key.cmp(&right.tie_key))
    });
    candidates.into_iter().next()
}

fn eligible_n_donors(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<Vec<(usize, Vec<ProtonCarrier>)>, MsReadyFailure> {
    let mut donors = Vec::new();
    for (row, atom) in topology.atoms.iter().enumerate() {
        if atom.element() != Element::N
            || atom.formal_charge() != 1
            || atom.is_aromatic()
            || atom.radical_electrons() != 0
            || topology.adjacency.neighbors_of(row).iter().any(|neighbor| {
                let bond = &topology.bonds[neighbor.bond.index()];
                bond.order() != BondOrder::Single || bond.is_aromatic() || bond.query().is_some()
            })
        {
            continue;
        }
        let carriers = available_carriers(topology, valence, row)?;
        if !carriers.is_empty() {
            donors.push((row, carriers));
        }
    }
    Ok(donors)
}

fn eligible_o_acceptors(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<Vec<(usize, usize)>, MsReadyFailure> {
    let mut acceptors = Vec::new();
    for (row, atom) in topology.atoms.iter().enumerate() {
        if atom.element() != Element::O
            || atom.formal_charge() != -1
            || atom.is_aromatic()
            || atom.radical_electrons() != 0
            || crate::rules::hydrogen_count(topology, valence, row)? != 0
        {
            continue;
        }
        let neighbors = topology.adjacency.neighbors_of(row);
        if neighbors.len() != 1 {
            continue;
        }
        let neighbor = neighbors[0];
        let parent = &topology.atoms[neighbor.atom_index];
        let bond = &topology.bonds[neighbor.bond.index()];
        if !matches!(parent.element(), Element::C | Element::P | Element::S)
            || bond.order() != BondOrder::Single
            || bond.is_aromatic()
            || bond.query().is_some()
        {
            continue;
        }
        acceptors.push((row, neighbor.atom_index));
    }
    Ok(acceptors)
}

fn flexible_valence(topology: &TopologyBlock) -> Result<ValenceAssignment, MsReadyFailure> {
    assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, false)
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))
}

fn strict_valence(topology: &TopologyBlock) -> Result<ValenceAssignment, MsReadyFailure> {
    assign_valence_for_topology(topology, ValenceModel::RdkitLike)
        .map_err(|error| failure(Kind::InvalidValence, error.to_string()))
}

fn isotope_signature(topology: &TopologyBlock) -> (Vec<(u8, Option<u16>)>, Vec<u16>) {
    let (mut atom_isotopes, tracked) = isotope_inventory(topology);
    atom_isotopes.retain(|(atomic_number, isotope)| *atomic_number != 1 || isotope.is_some());
    (atom_isotopes, tracked)
}

fn input_row_for_current(
    mapping: &TopologyMapping,
    current_row: usize,
    required: bool,
) -> Result<usize, MsReadyFailure> {
    let input_row = mapping
        .atoms
        .new_to_old
        .get(current_row)
        .copied()
        .flatten()
        .map(|id| id.index());
    match (input_row, required) {
        (Some(row), _) => Ok(row),
        (None, false) => Ok(usize::MAX),
        (None, true) => Err(failure(
            Kind::InvariantViolation,
            format!("current atom row {current_row} has no original source row"),
        )),
    }
}

fn component_ids(topology: &TopologyBlock) -> Vec<usize> {
    let mut ids = vec![usize::MAX; topology.atoms.len()];
    let mut next_id = 0usize;
    for start in 0..topology.atoms.len() {
        if ids[start] != usize::MAX {
            continue;
        }
        ids[start] = next_id;
        let mut queue = VecDeque::from([start]);
        while let Some(row) = queue.pop_front() {
            for neighbor in topology.adjacency.neighbors_of(row) {
                if ids[neighbor.atom_index] == usize::MAX {
                    ids[neighbor.atom_index] = next_id;
                    queue.push_back(neighbor.atom_index);
                }
            }
        }
        next_id += 1;
    }
    ids
}

fn radical_warnings(topology: &TopologyBlock) -> Vec<MsReadyWarning> {
    // BEGIN MSREADY REFERENCE §7
    // MSReady❗❗: A new conservative profile can leave stable radicals and permanent onium chemistry unchanged or flag them.
    // END MSREADY REFERENCE
    topology
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.radical_electrons() != 0)
        .map(|(row, atom)| MsReadyWarning {
            code: "stable_radical_preserved".to_owned(),
            message: format!(
                "preserved {} radical electron(s) at input atom row {row}",
                atom.radical_electrons()
            ),
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
mod charge_tests {
    use super::*;
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

    fn parse(smiles: &str) -> TopologyBlock {
        parse_smiles(smiles, &SmilesParseParams::default())
            .unwrap_or_else(|error| panic!("cannot parse {smiles}: {error}"))
            .topology
    }

    fn hydrogen_count_from_signature(signature: &str) -> u32 {
        signature
            .split(';')
            .find_map(|element| element.strip_prefix("1:"))
            .map_or(0, |count| count.parse().unwrap())
    }

    fn heavy_formula(signature: &str) -> Vec<&str> {
        signature
            .split(';')
            .filter(|element| !element.starts_with("1:"))
            .collect()
    }

    fn mapped_charge(topology: &TopologyBlock, event: &ProtonEvent, input_row: usize) -> i8 {
        let output_row = event.topology_mapping.atoms.old_to_new[input_row]
            .expect("charge-center source row remains mapped")
            .index();
        topology.atoms[output_row].formal_charge()
    }

    #[test]
    fn internal_dt_carrier_selection_is_input_row_independent() {
        let first = parse("C[N+]([2H])([3H])CC(=O)[O-]");
        let reordered = parse("C[N+]([3H])([2H])CC(=O)[O-]");
        let left = normalize_proton_charges(&first, &MsReadyParams::default()).unwrap();
        let right = normalize_proton_charges(&reordered, &MsReadyParams::default()).unwrap();
        assert_eq!(left.events[0].kind, ProtonEventKind::InternalTransfer);
        assert_eq!(right.events[0].kind, ProtonEventKind::InternalTransfer);
        assert_eq!(
            canonical_graph_key(&left.topology).unwrap(),
            canonical_graph_key(&right.topology).unwrap()
        );
        assert_eq!(
            crate::identify::standard_identity(&left.topology)
                .unwrap()
                .full_key,
            crate::identify::standard_identity(&right.topology)
                .unwrap()
                .full_key
        );
        assert_eq!(isotope_signature(&left.topology), isotope_signature(&first));
        assert_eq!(
            isotope_signature(&right.topology),
            isotope_signature(&reordered)
        );
        assert!(
            normalize_proton_charges(&left.topology, &MsReadyParams::default())
                .unwrap()
                .events
                .is_empty()
        );
    }

    #[test]
    fn charge_edit_and_candidate_limits_are_typed_not_partial_success() {
        let input = parse("C[NH3+].CC(=O)[O-]");
        let snapshot = input.clone();
        let params = MsReadyParams {
            max_edits: 1,
            ..Default::default()
        };
        assert_eq!(
            normalize_proton_charges(&input, &params).unwrap_err().kind,
            Kind::TransformationLimit
        );
        assert_eq!(input, snapshot);
        for smiles in [
            "C[NH3+].C[NH3+]",
            "C[O-].C[O-]",
            "C[N+]([2H])([3H])CC(=O)[O-]",
        ] {
            let params = MsReadyParams {
                max_matches_per_sweep: 1,
                ..Default::default()
            };
            assert_eq!(
                normalize_proton_charges(&parse(smiles), &params)
                    .unwrap_err()
                    .kind,
                Kind::TransformationLimit,
                "{smiles}"
            );
        }
    }

    #[test]
    fn tracked_dt_internal_transfer_enumerates_both_isotopes_and_matches_explicit_products() {
        let explicit =
            normalize_proton_charges(&parse("C[N+]([2H])([3H])CC(=O)[O-]"), &Default::default())
                .unwrap();
        let expected = canonical_graph_key(&explicit.topology).unwrap();
        for isotopes in [vec![2, 3], vec![3, 2]] {
            let mut tracked = parse("C[N+]CC(=O)[O-]");
            tracked.atoms[1].set_tracked_isotopic_hydrogens(isotopes);
            let snapshot = tracked.clone();
            let result = normalize_proton_charges(&tracked, &Default::default()).unwrap();
            assert_eq!(tracked, snapshot);
            assert_eq!(result.events[0].kind, ProtonEventKind::InternalTransfer);
            assert_eq!(canonical_graph_key(&result.topology).unwrap(), expected);
            assert_eq!(
                isotope_signature(&result.topology),
                isotope_signature(&tracked)
            );
            assert_eq!(
                crate::identify::standard_identity(&result.topology)
                    .unwrap()
                    .full_key,
                crate::identify::standard_identity(&explicit.topology)
                    .unwrap()
                    .full_key
            );
        }
    }

    #[test]
    fn group_and_charge_stages_share_the_input_wide_edit_budget() {
        let input = parse("C=C(O)C.C[NH3+]");
        let params = MsReadyParams {
            max_edits: 1,
            ..Default::default()
        };
        let mut budget = TransformationBudget::new(&params);
        let group = crate::rewrite::normalize_group_rules_with_budget(&input, &params, &mut budget)
            .unwrap();
        assert_eq!(group.edits.len(), 1);
        assert_eq!(
            normalize_proton_charges_with_budget(&group.topology, &params, &mut budget)
                .unwrap_err()
                .kind,
            Kind::TransformationLimit
        );
        let params = MsReadyParams {
            max_edits: 2,
            ..Default::default()
        };
        let mut budget = TransformationBudget::new(&params);
        let group = crate::rewrite::normalize_group_rules_with_budget(&input, &params, &mut budget)
            .unwrap();
        let charge =
            normalize_proton_charges_with_budget(&group.topology, &params, &mut budget).unwrap();
        assert_eq!(charge.events.len(), 1);
        assert_eq!(total_charge(&charge.topology), 0);
        let one_site = MsReadyParams {
            max_edits: 1,
            max_matches_per_sweep: 1,
            ..Default::default()
        };
        assert_eq!(
            normalize_proton_charges(&parse("C[NH3+]"), &one_site)
                .unwrap()
                .events
                .len(),
            1
        );
        let rematch_limit = MsReadyParams {
            max_matches_per_sweep: 2,
            ..Default::default()
        };
        assert_eq!(
            normalize_proton_charges(&parse("C[NH3+].C[NH3+]"), &rematch_limit)
                .unwrap_err()
                .kind,
            Kind::TransformationLimit,
            "candidate count is cumulative across rematches"
        );
    }

    #[test]
    fn primary_secondary_and_tertiary_ammonium_exchange_ordinary_h() {
        for smiles in ["C[NH3+]", "C[NH2+]C", "C[NH+](C)C"] {
            let input = parse(smiles);
            let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
            assert_eq!(result.events.len(), 1, "{smiles}");
            assert_eq!(
                result.events[0].kind,
                ProtonEventKind::ReservoirDeprotonation
            );
            assert_eq!(result.events[0].delta_hydrogens, -1);
            assert_eq!(result.events[0].delta_charge, -1);
            assert_eq!(
                hydrogen_count_from_signature(&result.events[0].formula_before),
                hydrogen_count_from_signature(&result.events[0].formula_after) + 1
            );
            assert_eq!(
                heavy_formula(&result.events[0].formula_before),
                heavy_formula(&result.events[0].formula_after)
            );
            assert_eq!(
                mapped_charge(
                    &result.topology,
                    &result.events[0],
                    result.events[0].donor_input_row.unwrap()
                ),
                0
            );
            assert_eq!(
                result.events[0].total_charge_after,
                result.events[0].total_charge_before - 1
            );
            assert_eq!(
                total_charge(&result.topology),
                result.events[0].total_charge_after
            );
        }
    }

    #[test]
    fn carboxylate_and_alkoxide_exchange_a_reservoir_proton() {
        for smiles in ["CC(=O)[O-]", "C[O-]"] {
            let input = parse(smiles);
            let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
            assert_eq!(result.events.len(), 1, "{smiles}");
            assert_eq!(result.events[0].kind, ProtonEventKind::ReservoirProtonation);
            assert_eq!(result.events[0].delta_hydrogens, 1);
            assert_eq!(result.events[0].delta_charge, 1);
            assert_eq!(
                hydrogen_count_from_signature(&result.events[0].formula_after),
                hydrogen_count_from_signature(&result.events[0].formula_before) + 1
            );
            assert_eq!(
                heavy_formula(&result.events[0].formula_before),
                heavy_formula(&result.events[0].formula_after)
            );
            assert_eq!(
                mapped_charge(
                    &result.topology,
                    &result.events[0],
                    result.events[0].acceptor_input_row.unwrap()
                ),
                0
            );
        }
    }

    #[test]
    fn glycine_uses_balanced_internal_transfer_before_reservoir_events() {
        let input = parse("[NH3+]CC(=O)[O-]");
        let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.events.len(), 1);
        assert_eq!(result.events[0].kind, ProtonEventKind::InternalTransfer);
        assert_eq!(result.events[0].delta_hydrogens, 0);
        assert_eq!(result.events[0].delta_charge, 0);
        assert_eq!(
            result.events[0].formula_before,
            result.events[0].formula_after
        );
        assert_eq!(
            mapped_charge(
                &result.topology,
                &result.events[0],
                result.events[0].donor_input_row.unwrap()
            ),
            0
        );
        assert_eq!(
            mapped_charge(
                &result.topology,
                &result.events[0],
                result.events[0].acceptor_input_row.unwrap()
            ),
            0
        );
        assert_eq!(total_charge(&result.topology), 0);
    }

    #[test]
    fn permanent_and_aromatic_ammonium_nitro_and_azide_charges_are_preserved() {
        for smiles in [
            "C[N+](C)(C)C",
            "[nH+]1ccccc1",
            "C[N+](=O)[O-]",
            "C[N-][N+]#N",
        ] {
            let input = parse(smiles);
            let result = normalize_proton_charges(&input, &MsReadyParams::default())
                .unwrap_or_else(|error| panic!("{smiles}: {error}"));
            assert!(result.events.is_empty(), "{smiles}");
            assert_eq!(result.topology.atoms.len(), input.atoms.len());
            assert_eq!(total_charge(&result.topology), total_charge(&input));
        }
    }

    #[test]
    fn multiple_sites_pair_within_component_then_exchange_residual_protons() {
        let input = parse("C[NH+](C)CC(CC[NH3+])C(=O)[O-]");
        let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.events.len(), 2);
        assert_eq!(result.events[0].kind, ProtonEventKind::InternalTransfer);
        assert_eq!(
            result.events[1].kind,
            ProtonEventKind::ReservoirDeprotonation
        );
        assert_eq!(result.events[1].delta_hydrogens, -1);
        assert_eq!(result.events[1].delta_charge, -1);
        assert_eq!(total_charge(&result.topology), 0);
    }

    #[test]
    fn explicit_hydrogen_rows_are_disposed_and_mapped_through_folding() {
        let input = parse("[N+]([H])([H])([H])CC(=O)[O-]");
        let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.events.len(), 1);
        assert_eq!(result.events[0].kind, ProtonEventKind::InternalTransfer);
        assert_eq!(result.hydrogen_dispositions.len(), 3);
        assert!(
            result
                .hydrogen_dispositions
                .iter()
                .all(|disposition| matches!(
                    disposition.kind,
                    HydrogenDispositionKind::FoldedInto { .. }
                ))
        );
        assert!(
            result
                .topology_mapping
                .atoms
                .old_to_new
                .iter()
                .all(|mapped| mapped.is_none()
                    || mapped.is_some_and(|id| id.index() < result.topology.atoms.len()))
        );
        for (bond_row, bond) in input.bonds.iter().enumerate() {
            let touches_hydrogen = input.atoms[bond.begin().index()].atomic_number() == 1
                || input.atoms[bond.end().index()].atomic_number() == 1;
            assert_eq!(
                result.topology_mapping.bonds.old_to_new[bond_row].is_none(),
                touches_hydrogen
            );
        }
    }

    #[test]
    fn explicit_ordinary_reservoir_hydrogen_is_deleted_and_mapped() {
        let input = parse("C[N+]([H])([H])C");
        let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.events.len(), 1);
        assert_eq!(
            result.events[0].kind,
            ProtonEventKind::ReservoirDeprotonation
        );
        assert!(
            result
                .hydrogen_dispositions
                .iter()
                .any(|disposition| matches!(
                    disposition.kind,
                    HydrogenDispositionKind::DeletedToReservoir
                ))
        );
        assert!(
            result
                .hydrogen_dispositions
                .iter()
                .any(|disposition| matches!(
                    disposition.kind,
                    HydrogenDispositionKind::FoldedInto { .. }
                ))
        );
        for (bond_row, bond) in input.bonds.iter().enumerate() {
            let touches_hydrogen = input.atoms[bond.begin().index()].atomic_number() == 1
                || input.atoms[bond.end().index()].atomic_number() == 1;
            assert_eq!(
                result.topology_mapping.bonds.old_to_new[bond_row].is_none(),
                touches_hydrogen
            );
        }
    }

    #[test]
    fn explicit_dt_transfers_internally_but_reservoir_exchange_is_blocked() {
        for (isotope, smiles) in [
            (2, "C[N+]([2H])(C)CC(=O)[O-]"),
            (3, "C[N+]([3H])(C)CC(=O)[O-]"),
        ] {
            let input = parse(smiles);
            let result = normalize_proton_charges(&input, &MsReadyParams::default())
                .unwrap_or_else(|error| panic!("isotope {isotope}: {error}"));
            assert_eq!(result.events.len(), 1);
            assert_eq!(result.events[0].kind, ProtonEventKind::InternalTransfer);
            assert!(matches!(
                result.events[0].carrier,
                ProtonCarrierRecord::ExplicitAtom { isotope: Some(value), .. } if value == isotope
            ));
            assert!(
                result
                    .topology
                    .atoms
                    .iter()
                    .any(|atom| { atom.atomic_number() == 1 && atom.isotope() == Some(isotope) })
            );
        }

        let mut tracked = parse("C[NH+](C)CC(=O)[O-]");
        tracked.atoms[1].set_explicit_hydrogens(0);
        tracked.atoms[1].set_tracked_isotopic_hydrogens(vec![2]);
        let tracked_result = normalize_proton_charges(&tracked, &MsReadyParams::default())
            .unwrap_or_else(|error| panic!("tracked D: {error}"));
        assert_eq!(tracked_result.events.len(), 1);
        assert_eq!(
            tracked_result.events[0].kind,
            ProtonEventKind::InternalTransfer
        );
        assert_eq!(
            tracked_result.events[0].carrier,
            ProtonCarrierRecord::TrackedIsotope { isotope: 2 }
        );
        assert!(
            tracked_result
                .topology
                .atoms
                .iter()
                .any(|atom| { atom.tracked_isotopic_hydrogens().contains(&2) })
        );

        let labelled_only = parse("C[N+]([2H])(C)C");
        let blocked = normalize_proton_charges(&labelled_only, &MsReadyParams::default()).unwrap();
        assert_eq!(total_charge(&blocked.topology), 1);
        assert!(
            blocked
                .events
                .iter()
                .any(|event| event.kind == ProtonEventKind::LabelledReservoirBlocked)
        );
        assert!(
            blocked
                .warnings
                .iter()
                .any(|warning| warning.code == "labelled_hydrogen_not_removed_from_reservoir")
        );
    }

    #[test]
    fn nitrate_comparison_templates_do_not_become_proton_events() {
        for smiles in ["[OH][N+](=O)[O-]", "[O-][N+](=O)[O-]"] {
            let input = parse(smiles);
            let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
            assert!(result.events.is_empty(), "{smiles}");
            assert_eq!(total_charge(&result.topology), total_charge(&input));
        }
    }

    #[test]
    fn radicals_are_preserved_and_flagged_without_charge_repair() {
        let mut input = parse("[CH3].CC(=O)[O-]");
        input.atoms[0].set_radical_electrons(1);
        let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(
            result.topology.atoms[0].radical_electrons(),
            input.atoms[0].radical_electrons()
        );
        assert!(
            result
                .warnings
                .iter()
                .any(|warning| warning.code == "stable_radical_preserved")
        );
        assert_eq!(total_charge(&result.topology), 0);
    }

    #[test]
    fn disconnected_ammonium_and_carboxylate_do_not_transfer_between_components() {
        let input = parse("C[NH3+].CC(=O)[O-]");
        let result = normalize_proton_charges(&input, &MsReadyParams::default()).unwrap();
        assert_eq!(result.events.len(), 2);
        assert_eq!(
            result.events[0].kind,
            ProtonEventKind::ReservoirDeprotonation
        );
        assert_eq!(result.events[1].kind, ProtonEventKind::ReservoirProtonation);
        assert_eq!(total_charge(&result.topology), 0);
    }
}
