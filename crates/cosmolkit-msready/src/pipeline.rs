use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, CoordinateBlock, MoleculeProperties, TopologyBlock, TopologyMapping,
};
use cosmolkit_smiles::{SmilesParseParams, SmilesRecord, parse_smiles};

use crate::charge::{ProtonCarrierRecord, ProtonEvent, ProtonEventKind};
use crate::finalize::{FinalizedComponent, finalize_component, writer_params};
use crate::rewrite::{RuleEdit, TransformationBudget, canonical_graph_key, compose_mappings};
use crate::*;

pub const RULE_SHA256: &str = "385df67d5c2fcff4c91c9c1e2a08939531bb79afa8f2f69a3a64a2d89ddd0610";
pub const EXCLUSION_SHA256: &str =
    "5b11ff6fbffb5ee0169e73c9428d85344598eac33a76b4e96b8aad06e93c4322";
pub const EXCLUSION_VARIANTS_SHA256: &str =
    "61184edcae96bf9e7267b26088c417d8134f181ffbb9eda5b68f5d3ef96cac47";
pub const NITRATE_TEMPLATES_SHA256: &str =
    "a73adc8cc780c2ddf8cb94c09ff7cb7c5c043fb0a6e2e450830f27284cf1bbb2";

/// Standardize a detached definite graph without mutating its carriers or properties.
pub fn standardize(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &MsReadyParams,
) -> Result<MsReadyResult, MsReadyFailure> {
    run(
        topology,
        properties,
        params,
        input_link(params, None),
        false,
    )
}

/// Parse with the SMILES owner without generic chemistry cleanup or H removal.
/// Actual input text and caller archival metadata are distinct provenance fields.
pub fn standardize_smiles(
    text: &str,
    params: &MsReadyParams,
) -> Result<MsReadyResult, MsReadyFailure> {
    let link = input_link(params, Some(text));
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        skip_cleanup: true,
        parse_name: false,
        ..Default::default()
    };
    let record = parse_smiles(text, &parser).map_err(|error| {
        let mut failure = MsReadyFailure::new(
            MsReadyFailureKind::InvalidInput,
            MsReadyStage::Input,
            error.to_string(),
        );
        failure.evidence.input_link = link.clone();
        failure.evidence.location.owner_category =
            Some(format!("cosmolkit_smiles::SmilesParseError::{error:?}"));
        failure
    })?;
    run(&record.topology, &record.properties, params, link, true)
}

fn input_link(params: &MsReadyParams, source: Option<&str>) -> MsReadyInputLink {
    MsReadyInputLink {
        record_id: params.record_id.clone(),
        original_smiles: params.original_smiles.clone(),
        source_smiles: source.map(str::to_owned),
    }
}

#[derive(Default)]
struct Progress {
    exclusions: Vec<MsReadyExclusion>,
    events: Vec<MsReadyEvent>,
    warnings: Vec<MsReadyWarning>,
    location: MsReadyLocation,
    input_components: Vec<usize>,
}

struct Fragment {
    topology: TopologyBlock,
    mapping: TopologyMapping,
}

fn run(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &MsReadyParams,
    link: MsReadyInputLink,
    initialize_radicals: bool,
) -> Result<MsReadyResult, MsReadyFailure> {
    // Approved concept reference: MSREADY_SMILES_implementation_reference.md;
    // SHA-256 0065003553dbc57f77632a4430bed30e3e189803b828f3dde08331443677900d.
    // BEGIN MSREADY REFERENCE §11
    // MSReady❗❗: result.merge_identical_products_keep_all_origin_links()
    // END MSREADY REFERENCE
    // Native behavior: frozen exact exclusions/repairs, whole-profile fixed point,
    // shared input-wide budget, isotope/stereo policy and complete atom accounting.
    // This is not execution equivalence to the historical one-pass workflow.
    // Native complexity: ordered O(P log P) identity aggregation; detached owner
    // stage clones, canonical ranking, graph matching and identity retain their
    // owner costs. No historical performance-equivalence claim.
    let original = MsReadyOriginal {
        topology: topology.clone(),
        properties: properties.clone(),
    };
    let mut progress = Progress::default();
    let result = (|| {
        validate_syntax(topology, properties, params)?;
        let raw_parts = fragments(
            topology,
            topology,
            &TopologyMapping::identity(topology.atoms.len(), topology.bonds.len()),
        )?;
        let mut input_components = vec![0; topology.atoms.len()];
        for (index, part) in raw_parts.iter().enumerate() {
            for id in part.mapping.atoms.new_to_old.iter().flatten() {
                input_components[id.index()] = index;
            }
        }
        progress.input_components = input_components.clone();
        // Exact raw recognition precedes any retained text radical initialization
        // and provides completed comparison evidence if a later stage fails.
        let mut initial_matches = Vec::new();
        for (index, part) in raw_parts.iter().enumerate() {
            let matches = record_exclusions(
                topology,
                part,
                index,
                0,
                &input_components,
                &BTreeMap::new(),
            )?;
            initial_matches.push(!matches.is_empty());
            progress.exclusions.extend(matches);
        }
        let mut working = topology.clone();
        let mut budget = TransformationBudget::new(params);
        budget.begin_sweep();
        if initialize_radicals {
            // Observation is a bounded preflight, not another chemical sweep:
            // reserve the shared disconnection allowance across components but
            // charge actual candidates only when the initialized graph matches.
            let mut observation_remaining = budget.remaining_matches();
            for (index, part) in raw_parts.iter().enumerate() {
                if initial_matches[index] {
                    continue;
                }
                progress.location = location_for(part, index, &input_components);
                // Supervisor-approved observation routing for the exact bonded
                // (+1,-1) representation repair. Owner radical inference must
                // observe the validated spectator state, not Li+ with a pending
                // covalent depiction. No chemical projection event is emitted.
                // All C/O/H/radical/isotope attributes remain untouched, and the
                // initialized analysis graph reruns the actual frozen guards.
                let candidates = crate::disconnect::alkali_oxygen_candidates_limited(
                    &part.topology,
                    observation_remaining,
                )
                .map_err(|error| {
                    let mut error = remap_error(error, &part.mapping);
                    error.stage = MsReadyStage::Disconnection;
                    error
                })?;
                observation_remaining -= candidates.len();
                let mut observation = part
                    .topology
                    .begin_batch_edit()
                    .map_err(|error| invariant(&error.to_string()))?;
                for candidate in candidates {
                    if candidate.metal_charge_delta == 0 {
                        observation
                            .remove_bond(cosmolkit_model::BondId::new(candidate.bond_index))
                            .map_err(|error| invariant(&error.to_string()))?;
                    }
                }
                let (observation, _) = observation
                    .finish()
                    .map_err(|error| invariant(&error.to_string()))?;
                let assignment =
                    cosmolkit_core::assign_radicals(&observation).map_err(|error| {
                        let mut failure = MsReadyFailure::new(
                            MsReadyFailureKind::OwnerError,
                            MsReadyStage::Input,
                            error.to_string(),
                        );
                        failure.evidence.location.owner_category =
                            Some(format!("cosmolkit_core::RadicalError::{error:?}"));
                        failure
                    })?;
                for (row, atom) in part.topology.atoms.iter().enumerate() {
                    if atom.radical_electrons() != 0 || assignment.radical_electrons[row] == 0 {
                        continue;
                    }
                    let input_row = mapped_atom(&part.mapping, row)?;
                    working.atoms[input_row]
                        .set_radical_electrons(assignment.radical_electrons[row]);
                    let mut event = empty_event(
                        MsReadyEventKind::ParserRadicalInitialization,
                        MsReadyStage::Input,
                        "owner_parser_radical_initialization",
                        part.mapping.clone(),
                    );
                    event.location = progress.location.clone();
                    event.location.rule_id = Some(event.rule_id.clone());
                    event.location.atom_indices = vec![input_row];
                    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
                        &part.topology,
                        cosmolkit_core::ValenceModel::RdkitLike,
                        false,
                    )
                    .map_err(|error| {
                        MsReadyFailure::new(
                            MsReadyFailureKind::OwnerError,
                            MsReadyStage::Input,
                            error.to_string(),
                        )
                    })?;
                    event.formula_before =
                        crate::rewrite::formula_signature(&part.topology, &valence)?;
                    event.formula_after = event.formula_before.clone();
                    event.total_charge_before = crate::rewrite::total_charge(&part.topology);
                    event.total_charge_after = event.total_charge_before;
                    event.radical_electrons_before = Some(0);
                    event.radical_electrons_after = Some(assignment.radical_electrons[row]);
                    progress.events.push(event);
                    progress.warnings.push(MsReadyWarning { code: "parser_radical_initialized".into(), message: format!("owner inferred {} radical electron(s) at original input atom {input_row}; raw input snapshot preserved", assignment.radical_electrons[row]) });
                }
            }
        }
        progress.location = MsReadyLocation::default();
        let disconnected = crate::disconnect::disconnect_all_alkali_oxygen_observed(
            &working,
            Some(&mut budget),
            |_, edit| {
                let mut event = empty_event(
                    MsReadyEventKind::AlkaliDisconnection,
                    MsReadyStage::Disconnection,
                    "alkali_oxygen_disconnection",
                    edit.topology_mapping.clone(),
                );
                event.location = MsReadyLocation {
                    input_component: Some(input_components[edit.metal_atom]),
                    rule_id: Some(event.rule_id.clone()),
                    atom_indices: vec![edit.metal_atom, edit.oxygen_atom],
                    bond_indices: vec![edit.bond_index],
                    ..Default::default()
                };
                event.formula_before = edit.formula_before.clone();
                event.formula_after = edit.formula_after.clone();
                event.total_charge_before = edit.total_charge_before;
                event.total_charge_after = edit.total_charge_after;
                event.atom_charge_deltas = vec![
                    (edit.metal_atom, edit.metal_charge_delta),
                    (edit.oxygen_atom, edit.oxygen_charge_delta),
                ];
                event.hydrogen_counts_before = edit
                    .hydrogen_counts_by_input_atom
                    .iter()
                    .copied()
                    .enumerate()
                    .collect();
                progress.events.push(event);
                Ok(())
            },
        )
        .map_err(|mut error| {
            error.stage = MsReadyStage::Disconnection;
            error
        })?;
        let parts = fragments(
            topology,
            &disconnected.topology,
            &disconnected.topology_mapping,
        )?;
        let mut post_disconnection_exclusions = Vec::new();
        let mut retained = Vec::new();
        // Compare ALL components before strict retained-domain failures; organic
        // counterions are retained unless an exact catalogue graph matches.
        for (index, part) in parts.into_iter().enumerate() {
            let matches = record_exclusions(
                topology,
                &part,
                index,
                0,
                &input_components,
                &BTreeMap::new(),
            )?;
            if matches.is_empty() {
                retained.push((index, part));
            } else {
                post_disconnection_exclusions.extend(matches);
            }
        }
        progress.exclusions = post_disconnection_exclusions;
        let mut products = BTreeMap::<(String, String), MsReadyProduct>::new();
        let mut h_dispositions = BTreeMap::<usize, MsReadyHydrogenDispositionKind>::new();
        for (component_index, part) in retained {
            progress.location = location_for(&part, component_index, &input_components);
            let mut current = part.topology;
            let mut mapping = part.mapping;
            // Strict checks may use only the already-adopted guarded neutral
            // nitro validation copy; the actual later edit is always logged.
            let validation = crate::rewrite::nitro_validation_graph(&current)?;
            crate::validate::strict_owner_valence(&validation)
                .map_err(|e| remap_error(e, &mapping))?;
            retained_metal_states(&validation).map_err(|e| remap_error(e, &mapping))?;
            let mut seen = BTreeSet::from([canonical_graph_key(&current)?]);
            loop {
                let before = canonical_graph_key(&current)?;
                let group_input_mapping = mapping.clone();
                let groups = crate::rewrite::normalize_group_rules_observed(
                    &current,
                    params,
                    &mut budget,
                    |edit| {
                        let event =
                            group_event(topology, &group_input_mapping, edit, &progress.location)?;
                        progress.events.push(event);
                        Ok(())
                    },
                )
                .map_err(|e| remap_error(e, &mapping))?;
                mapping = compose_mappings(
                    &mapping,
                    &groups.topology_mapping,
                    topology,
                    &groups.topology,
                )?;
                current = groups.topology;
                crate::validate::strict_owner_valence(&current)
                    .map_err(|e| remap_error(e, &mapping))?;
                let charge_input_mapping = mapping.clone();
                let charge = crate::charge::normalize_proton_charges_observed(
                    &current,
                    params,
                    &mut budget,
                    |edit| {
                        progress.events.push(proton_event(
                            topology,
                            &charge_input_mapping,
                            edit,
                            &progress.location,
                        )?);
                        Ok(())
                    },
                )
                .map_err(|e| remap_error(e, &mapping))?;
                remember_hydrogens(&charge.hydrogen_dispositions, &mapping, &mut h_dispositions)?;
                for warning in charge.warnings {
                    let contextual = MsReadyWarning {
                        code: warning.code,
                        message: format!(
                            "input component {:?}, post-disconnection component {component_index}, original atoms {:?}: {}",
                            progress.location.input_component,
                            progress.location.atom_indices,
                            warning.message
                        ),
                    };
                    if !progress.warnings.contains(&contextual) {
                        progress.warnings.push(contextual);
                    }
                }
                mapping = compose_mappings(
                    &mapping,
                    &charge.topology_mapping,
                    topology,
                    &charge.topology,
                )?;
                current = charge.topology;
                let after = canonical_graph_key(&current)?;
                if !visit_profile_state(&before, after, &mut seen)? {
                    break;
                }
            }
            // Connectivity is re-acquired through the owner, even though the
            // frozen rules currently preserve heavy-component connectivity.
            for (final_index, final_part) in fragments(topology, &current, &mapping)?
                .into_iter()
                .enumerate()
            {
                let matches = record_exclusions(
                    topology,
                    &final_part,
                    component_index,
                    final_index,
                    &input_components,
                    &h_dispositions,
                )?;
                if !matches.is_empty() {
                    progress.exclusions.extend(matches);
                    continue;
                }
                let record = SmilesRecord {
                    topology: final_part.topology,
                    coordinates: CoordinateBlock::default(),
                    properties: MoleculeProperties::default(),
                };
                let finalized = finalize_component(topology, &record, &final_part.mapping)
                    .map_err(|e| remap_error(e, &final_part.mapping))?;
                remember_hydrogens(
                    &finalized.hydrogen_dispositions,
                    &TopologyMapping::identity(topology.atoms.len(), topology.bonds.len()),
                    &mut h_dispositions,
                )?;
                for (sub_index, sub_part) in fragments(
                    topology,
                    &finalized.record.topology,
                    &finalized.original_to_final,
                )?
                .into_iter()
                .enumerate()
                {
                    let matches = record_exclusions(
                        topology,
                        &sub_part,
                        component_index,
                        final_index + sub_index,
                        &input_components,
                        &h_dispositions,
                    )?;
                    if !matches.is_empty() {
                        progress.exclusions.extend(matches);
                        continue;
                    }
                    let final_value =
                        if sub_part.topology.atoms.len() == finalized.record.topology.atoms.len() {
                            finalized.clone()
                        } else {
                            finalize_component(
                                topology,
                                &SmilesRecord {
                                    topology: sub_part.topology,
                                    coordinates: CoordinateBlock::default(),
                                    properties: MoleculeProperties::default(),
                                },
                                &sub_part.mapping,
                            )?
                        };
                    let origin = origin_for(
                        topology,
                        &final_value.original_to_final,
                        component_index,
                        final_index + sub_index,
                        &input_components,
                        &h_dispositions,
                    )?;
                    insert_product(topology, &mut products, final_value, origin)?;
                }
            }
        }
        let products: Vec<_> = products.into_values().collect();
        let dispositions =
            account_atoms(topology, &products, &progress.exclusions, &h_dispositions)?;
        if products.is_empty() && progress.exclusions.is_empty() {
            return Err(invariant("no product without an explicit exclusion reason"));
        }
        let status = if products.is_empty() {
            MsReadyStatus::NoEligibleComponent
        } else {
            MsReadyStatus::Products(products)
        };
        Ok(MsReadyResult {
            original: original.clone(),
            input_link: link.clone(),
            run_info: run_info(),
            status,
            exclusions: progress.exclusions.clone(),
            events: progress.events.clone(),
            warnings: progress.warnings.clone(),
            atom_dispositions: dispositions,
        })
    })();
    result.map_err(|mut error: MsReadyFailure| {
        error.original = Some(original);
        error.evidence.input_link = link;
        let location = &mut error.evidence.location;
        if location.component_index.is_none() {
            location.component_index = progress.location.component_index;
        }
        if location.input_component.is_none() {
            location.input_component = location
                .atom_indices
                .first()
                .and_then(|row| progress.input_components.get(*row).copied())
                .or(progress.location.input_component);
        }
        if location.atom_indices.is_empty() {
            location.atom_indices = progress.location.atom_indices;
        }
        error.evidence.exclusions = progress.exclusions;
        error.evidence.events = progress.events;
        error.evidence.warnings = progress.warnings;
        error
    })
}

pub fn run_info() -> MsReadyRunInfo {
    MsReadyRunInfo {
        profile_id: PROFILE_ID,
        profile_version: PROFILE_VERSION,
        behavior: BEHAVIOR,
        rule_sha256: RULE_SHA256,
        exclusion_sha256: EXCLUSION_SHA256,
        exclusion_variants_sha256: EXCLUSION_VARIANTS_SHA256,
        nitrate_templates_sha256: NITRATE_TEMPLATES_SHA256,
        inchi_engine_version: crate::identify::INCHI_ENGINE_VERSION,
        writer_params: writer_params(),
    }
}

fn fragments(
    original: &TopologyBlock,
    current: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Result<Vec<Fragment>, MsReadyFailure> {
    crate::components::split_mapped_components(
        current,
        &CoordinateBlock::default(),
        &MoleculeProperties::default(),
    )?
    .into_iter()
    .map(|part| {
        let composed = compose_mappings(mapping, &part.topology_mapping, original, &part.topology)?;
        Ok(Fragment {
            topology: part.topology,
            mapping: composed,
        })
    })
    .collect()
}
fn mapped_atom(mapping: &TopologyMapping, row: usize) -> Result<usize, MsReadyFailure> {
    mapping
        .atoms
        .new_to_old
        .get(row)
        .and_then(|id| id.map(|id| id.index()))
        .ok_or_else(|| invariant("required original atom mapping absent"))
}
fn remap_error(mut error: MsReadyFailure, mapping: &TopologyMapping) -> MsReadyFailure {
    error.evidence.location.atom_indices = error
        .evidence
        .location
        .atom_indices
        .iter()
        .filter_map(|row| {
            mapping
                .atoms
                .new_to_old
                .get(*row)
                .and_then(|id| id.map(|id| id.index()))
        })
        .collect();
    error.evidence.location.bond_indices = error
        .evidence
        .location
        .bond_indices
        .iter()
        .filter_map(|row| {
            mapping
                .bonds
                .new_to_old
                .get(*row)
                .and_then(|id| id.map(|id| id.index()))
        })
        .collect();
    error
}
fn location_for(part: &Fragment, index: usize, components: &[usize]) -> MsReadyLocation {
    let atoms: Vec<_> = part
        .mapping
        .atoms
        .new_to_old
        .iter()
        .flatten()
        .map(|id| id.index())
        .collect();
    MsReadyLocation {
        input_component: atoms.first().map(|row| components[*row]),
        component_index: Some(index),
        atom_indices: atoms,
        bond_indices: part
            .mapping
            .bonds
            .new_to_old
            .iter()
            .flatten()
            .map(|id| id.index())
            .collect(),
        ..Default::default()
    }
}
fn origin_for(
    original: &TopologyBlock,
    mapping: &TopologyMapping,
    index: usize,
    final_index: usize,
    components: &[usize],
    hydrogens: &BTreeMap<usize, MsReadyHydrogenDispositionKind>,
) -> Result<MsReadyOrigin, MsReadyFailure> {
    let mut atom_indices: Vec<_> = mapping
        .atoms
        .new_to_old
        .iter()
        .flatten()
        .map(|id| id.index())
        .collect();
    let mut dispositions = Vec::new();
    for (&row, &kind) in hydrogens {
        if let MsReadyHydrogenDispositionKind::FoldedInto { input_parent_row } = kind {
            if mapping.atoms.old_to_new[input_parent_row].is_some() {
                atom_indices.push(row);
                dispositions.push(MsReadyHydrogenDisposition {
                    input_atom_row: row,
                    kind,
                });
            }
        }
    }
    for &row in &atom_indices {
        if original.atoms[row].atomic_number() == 1 {
            if let Some(id) = mapping.atoms.old_to_new[row] {
                let kind = original.atoms[row].isotope().map_or(
                    MsReadyHydrogenDispositionKind::Retained {
                        output_atom_row: id.index(),
                    },
                    |isotope| MsReadyHydrogenDispositionKind::RetainedIsotope {
                        output_atom_row: id.index(),
                        isotope,
                    },
                );
                dispositions.push(MsReadyHydrogenDisposition {
                    input_atom_row: row,
                    kind,
                });
            }
        }
    }
    atom_indices.sort_unstable();
    atom_indices.dedup();
    dispositions.sort_by_key(|d| d.input_atom_row);
    let input_component = atom_indices
        .first()
        .map(|row| components[*row])
        .ok_or_else(|| invariant("product/exclusion has no original component origin"))?;
    Ok(MsReadyOrigin {
        input_component,
        component_index: index,
        final_component_index: final_index,
        atom_indices,
        bond_indices: mapping
            .bonds
            .new_to_old
            .iter()
            .flatten()
            .map(|id| id.index())
            .collect(),
        topology_mapping: mapping.clone(),
        hydrogen_dispositions: dispositions,
    })
}
fn record_exclusions(
    original: &TopologyBlock,
    part: &Fragment,
    index: usize,
    final_index: usize,
    components: &[usize],
    hydrogens: &BTreeMap<usize, MsReadyHydrogenDispositionKind>,
) -> Result<Vec<MsReadyExclusion>, MsReadyFailure> {
    let matches = match_exact_exclusions(&part.topology)?;
    if matches.is_empty() {
        return Ok(Vec::new());
    }
    let origin = origin_for(
        original,
        &part.mapping,
        index,
        final_index,
        components,
        hydrogens,
    )?;
    matches
        .into_iter()
        .map(|matched| {
            let atom_mappings = matched
                .atom_mappings
                .iter()
                .map(|rows| {
                    rows.iter()
                        .map(|row| mapped_atom(&part.mapping, *row))
                        .collect()
                })
                .collect::<Result<_, _>>()?;
            let folded_hydrogen_mappings = matched
                .folded_hydrogen_mappings
                .iter()
                .map(|(h, p)| {
                    Ok((
                        mapped_atom(&part.mapping, *h)?,
                        mapped_atom(&part.mapping, *p)?,
                    ))
                })
                .collect::<Result<_, MsReadyFailure>>()?;
            Ok(MsReadyExclusion {
                source_id: matched.source_id,
                source_sha256: EXCLUSION_SHA256,
                source_field: matched.source_field,
                template_id: matched.template_id,
                adaptation_id: matched.adaptation_id,
                comparison_only: matched.comparison_only,
                matched_formula: matched
                    .matched_formula
                    .ok_or_else(|| invariant("matched exclusion lacks formula metadata"))?,
                matched_formal_charge: matched.matched_formal_charge,
                origin: origin.clone(),
                atom_mappings,
                folded_hydrogen_mappings,
                reason: "exact_whole_component_exclusion",
            })
        })
        .collect()
}
fn remember_hydrogens(
    dispositions: &[MsReadyHydrogenDisposition],
    stage_input: &TopologyMapping,
    output: &mut BTreeMap<usize, MsReadyHydrogenDispositionKind>,
) -> Result<(), MsReadyFailure> {
    for disposition in dispositions {
        let row = mapped_atom(stage_input, disposition.input_atom_row)?;
        let kind = match disposition.kind {
            MsReadyHydrogenDispositionKind::FoldedInto { input_parent_row } => {
                MsReadyHydrogenDispositionKind::FoldedInto {
                    input_parent_row: mapped_atom(stage_input, input_parent_row)?,
                }
            }
            MsReadyHydrogenDispositionKind::DeletedToReservoir => {
                MsReadyHydrogenDispositionKind::DeletedToReservoir
            }
            _ => continue,
        };
        if let Some(previous) = output.insert(row, kind) {
            if previous != kind {
                return Err(invariant(
                    "input hydrogen received conflicting dispositions",
                ));
            }
        }
    }
    Ok(())
}
fn empty_event(
    kind: MsReadyEventKind,
    stage: MsReadyStage,
    rule: &str,
    mapping: TopologyMapping,
) -> MsReadyEvent {
    MsReadyEvent {
        kind,
        stage,
        rule_id: rule.into(),
        appearance_order: None,
        location: MsReadyLocation::default(),
        donor_input_row: None,
        acceptor_input_row: None,
        carrier: MsReadyProtonCarrier::None,
        delta_hydrogens: 0,
        delta_charge: 0,
        atom_charge_deltas: Vec::new(),
        hydrogen_counts_before: Vec::new(),
        radical_electrons_before: None,
        radical_electrons_after: None,
        formula_before: String::new(),
        formula_after: String::new(),
        total_charge_before: 0,
        total_charge_after: 0,
        topology_mapping: mapping,
    }
}
fn group_event(
    original: &TopologyBlock,
    input: &TopologyMapping,
    edit: &RuleEdit,
    location: &MsReadyLocation,
) -> Result<MsReadyEvent, MsReadyFailure> {
    let mapping = event_mapping(original, input, &edit.topology_mapping)?;
    let mut event = empty_event(
        MsReadyEventKind::GroupRewrite,
        MsReadyStage::Normalization,
        edit.rule_id,
        mapping,
    );
    event.location = location.clone();
    event.location.rule_id = Some(edit.rule_id.into());
    event.location.atom_indices = edit
        .matched_atom_rows
        .iter()
        .map(|row| mapped_atom(input, *row))
        .collect::<Result<_, _>>()?;
    event.location.bond_indices = edit
        .matched_bond_rows
        .iter()
        .map(|row| {
            input.bonds.new_to_old[*row]
                .map(|id| id.index())
                .ok_or_else(|| invariant("group edit lost original bond"))
        })
        .collect::<Result<_, _>>()?;
    event.appearance_order = Some(edit.appearance_order);
    event.formula_before = edit.formula_before.clone();
    event.formula_after = edit.formula_after.clone();
    event.total_charge_before = edit.total_charge_before;
    event.total_charge_after = edit.total_charge_after;
    if let Some(transfer) = edit.hydrogen_transfer {
        event.donor_input_row = Some(mapped_atom(input, transfer.donor)?);
        event.acceptor_input_row = Some(mapped_atom(input, transfer.acceptor)?);
        event.carrier = match transfer.carrier {
            crate::rules::HydrogenCarrier::Implicit => MsReadyProtonCarrier::InferredOrdinary,
            crate::rules::HydrogenCarrier::Bracket => MsReadyProtonCarrier::BracketOrdinary,
            crate::rules::HydrogenCarrier::ExplicitAtom { atom, isotope, .. } => {
                MsReadyProtonCarrier::ExplicitAtom {
                    input_atom_row: mapped_atom(input, atom)?,
                    isotope,
                }
            }
            crate::rules::HydrogenCarrier::TrackedIsotope { isotope, .. } => {
                MsReadyProtonCarrier::TrackedIsotope { isotope }
            }
        };
    }
    Ok(event)
}
fn proton_event(
    original: &TopologyBlock,
    input: &TopologyMapping,
    edit: &ProtonEvent,
    location: &MsReadyLocation,
) -> Result<MsReadyEvent, MsReadyFailure> {
    let kind = match edit.kind {
        ProtonEventKind::InternalTransfer => MsReadyEventKind::InternalProtonTransfer,
        ProtonEventKind::ReservoirDeprotonation => MsReadyEventKind::ReservoirDeprotonation,
        ProtonEventKind::ReservoirProtonation => MsReadyEventKind::ReservoirProtonation,
        ProtonEventKind::LabelledReservoirBlocked => MsReadyEventKind::LabelledReservoirBlocked,
    };
    let mapping = event_mapping(original, input, &edit.topology_mapping)?;
    let mut event = empty_event(
        kind,
        MsReadyStage::Normalization,
        &format!("{:?}", edit.kind),
        mapping,
    );
    event.location = location.clone();
    event.location.rule_id = Some(event.rule_id.clone());
    event.donor_input_row = edit
        .donor_input_row
        .map(|row| mapped_atom(input, row))
        .transpose()?;
    event.acceptor_input_row = edit
        .acceptor_input_row
        .map(|row| mapped_atom(input, row))
        .transpose()?;
    event.location.atom_indices = event
        .donor_input_row
        .into_iter()
        .chain(event.acceptor_input_row)
        .collect();
    event.location.bond_indices.clear();
    event.carrier = match edit.carrier {
        ProtonCarrierRecord::InferredOrdinary => MsReadyProtonCarrier::InferredOrdinary,
        ProtonCarrierRecord::BracketOrdinary => MsReadyProtonCarrier::BracketOrdinary,
        ProtonCarrierRecord::ExplicitAtom {
            input_atom_row,
            isotope,
        } => MsReadyProtonCarrier::ExplicitAtom {
            input_atom_row: mapped_atom(input, input_atom_row)?,
            isotope,
        },
        ProtonCarrierRecord::TrackedIsotope { isotope } => {
            MsReadyProtonCarrier::TrackedIsotope { isotope }
        }
        ProtonCarrierRecord::Reservoir => MsReadyProtonCarrier::Reservoir,
        ProtonCarrierRecord::None => MsReadyProtonCarrier::None,
    };
    event.delta_hydrogens = edit.delta_hydrogens;
    event.delta_charge = edit.delta_charge;
    if !matches!(edit.kind, ProtonEventKind::LabelledReservoirBlocked) {
        if let Some(row) = event.donor_input_row {
            event.atom_charge_deltas.push((row, -1));
        }
        if let Some(row) = event.acceptor_input_row {
            event.atom_charge_deltas.push((row, 1));
        }
    }
    event.formula_before = edit.formula_before.clone();
    event.formula_after = edit.formula_after.clone();
    event.total_charge_before = edit.total_charge_before;
    event.total_charge_after = edit.total_charge_after;
    Ok(event)
}
fn event_mapping(
    original: &TopologyBlock,
    input: &TopologyMapping,
    event: &TopologyMapping,
) -> Result<TopologyMapping, MsReadyFailure> {
    crate::rewrite::compose_mapping_counts(
        input,
        event,
        original.atoms.len(),
        original.bonds.len(),
        event.atoms.new_to_old.len(),
        event.bonds.new_to_old.len(),
    )
}

fn insert_product(
    original: &TopologyBlock,
    products: &mut BTreeMap<(String, String), MsReadyProduct>,
    finalized: FinalizedComponent,
    mut origin: MsReadyOrigin,
) -> Result<(), MsReadyFailure> {
    let key = (
        finalized.identity.full_key.clone(),
        finalized.canonical_smiles.clone(),
    );
    if let Some(product) = products.get_mut(&key) {
        let mappings = crate::exclusions::exact_search_mappings(
            &finalized.record.topology,
            &product.topology,
        )?;
        let rows = mappings
            .first()
            .ok_or_else(|| invariant("same exact product identity lacks graph isomorphism"))?;
        let mut alignment =
            TopologyMapping::identity(rows.len(), finalized.record.topology.bonds.len());
        alignment.atoms.old_to_new = rows.iter().map(|row| Some(AtomId::new(*row))).collect();
        for (old, &new) in rows.iter().enumerate() {
            alignment.atoms.new_to_old[new] = Some(AtomId::new(old));
        }
        for bond in &finalized.record.topology.bonds {
            let begin = AtomId::new(rows[bond.begin().index()]);
            let end = AtomId::new(rows[bond.end().index()]);
            let target = product
                .topology
                .bonds
                .iter()
                .find(|b| {
                    (b.begin() == begin && b.end() == end) || (b.begin() == end && b.end() == begin)
                })
                .ok_or_else(|| invariant("isomorphic product bond missing"))?;
            alignment.bonds.old_to_new[bond.id().index()] = Some(target.id());
            alignment.bonds.new_to_old[target.id().index()] = Some(bond.id());
        }
        origin.topology_mapping = compose_mappings(
            &origin.topology_mapping,
            &alignment,
            original,
            &product.topology,
        )?;
        for disposition in &mut origin.hydrogen_dispositions {
            match &mut disposition.kind {
                MsReadyHydrogenDispositionKind::Retained { output_atom_row }
                | MsReadyHydrogenDispositionKind::RetainedIsotope {
                    output_atom_row, ..
                } => *output_atom_row = rows[*output_atom_row],
                _ => {}
            }
        }
        product.origins.push(origin);
        product.multiplicity = product.origins.len();
    } else {
        let identity = finalized.identity;
        products.insert(
            key,
            MsReadyProduct {
                topology: finalized.record.topology,
                canonical_smiles: finalized.canonical_smiles,
                inchi: identity.inchi,
                inchi_key: identity.full_key,
                identifier: MsReadyIdentifierDiagnostics {
                    engine_version: identity.engine_version,
                    return_code: identity.return_code,
                    message: identity.message,
                    log: identity.log,
                    diagnostics: identity.diagnostics,
                },
                formula: finalized.composition.formula,
                isotope_composition: finalized
                    .composition
                    .counts
                    .into_iter()
                    .map(
                        |((atomic_number, isotope), count)| MsReadyCompositionEntry {
                            atomic_number,
                            isotope,
                            count,
                        },
                    )
                    .collect(),
                formal_charge: finalized.composition.formal_charge,
                neutral_atom_mass: finalized.composition.neutral_atom_mass,
                origins: vec![origin],
                multiplicity: 1,
            },
        );
    }
    Ok(())
}
fn account_atoms(
    original: &TopologyBlock,
    products: &[MsReadyProduct],
    exclusions: &[MsReadyExclusion],
    hydrogens: &BTreeMap<usize, MsReadyHydrogenDispositionKind>,
) -> Result<Vec<MsReadyAtomDisposition>, MsReadyFailure> {
    let mut ledger = BTreeMap::new();
    for (product_index, product) in products.iter().enumerate() {
        for (origin_index, origin) in product.origins.iter().enumerate() {
            for &row in &origin.atom_indices {
                let kind = if let Some(id) = origin.topology_mapping.atoms.old_to_new[row] {
                    MsReadyAtomDispositionKind::Retained {
                        product_index,
                        origin_index,
                        output_atom_row: id.index(),
                    }
                } else if let Some(MsReadyHydrogenDispositionKind::FoldedInto {
                    input_parent_row,
                }) = hydrogens.get(&row)
                {
                    MsReadyAtomDispositionKind::FoldedHydrogen {
                        product_index,
                        origin_index,
                        input_parent_row: *input_parent_row,
                    }
                } else {
                    return Err(invariant("origin atom lacks output or H-fold disposition"));
                };
                if ledger.insert(row, kind).is_some() {
                    return Err(invariant("input atom belongs to multiple retained origins"));
                }
            }
        }
    }
    for (index, exclusion) in exclusions.iter().enumerate() {
        for &row in &exclusion.origin.atom_indices {
            match ledger
                .entry(row)
                .or_insert_with(|| MsReadyAtomDispositionKind::Excluded {
                    exclusion_indices: Vec::new(),
                }) {
                MsReadyAtomDispositionKind::Excluded { exclusion_indices } => {
                    exclusion_indices.push(index)
                }
                _ => return Err(invariant("input atom both retained and excluded")),
            }
        }
    }
    for (&row, kind) in hydrogens {
        if *kind == MsReadyHydrogenDispositionKind::DeletedToReservoir {
            if ledger
                .insert(row, MsReadyAtomDispositionKind::DeletedToReservoir)
                .is_some()
            {
                return Err(invariant("reservoir H has another disposition"));
            }
        }
    }
    if ledger.len() != original.atoms.len() || ledger.keys().copied().ne(0..original.atoms.len()) {
        return Err(invariant("input atom disposition ledger is incomplete"));
    }
    Ok(ledger
        .into_iter()
        .map(|(input_atom_row, kind)| MsReadyAtomDisposition {
            input_atom_row,
            kind,
        })
        .collect())
}
fn invariant(message: &str) -> MsReadyFailure {
    MsReadyFailure::new(
        MsReadyFailureKind::InvariantViolation,
        MsReadyStage::Accounting,
        message,
    )
}

fn visit_profile_state(
    before: &str,
    after: String,
    seen: &mut BTreeSet<String>,
) -> Result<bool, MsReadyFailure> {
    if after == before {
        return Ok(false);
    }
    if !seen.insert(after) {
        return Err(MsReadyFailure::new(
            MsReadyFailureKind::TransformationCycle,
            MsReadyStage::Normalization,
            "whole-profile group/proton cycle revisits a distinct state",
        ));
    }
    Ok(true)
}

#[cfg(test)]
mod pipeline_tests {
    use super::*;
    #[test]
    fn whole_profile_boundary_uses_real_normalized_graph_states_and_ignores_no_ops() {
        let enol = parse_smiles("C=C(O)C", &SmilesParseParams::default()).unwrap();
        let ketone = parse_smiles("CC(=O)C", &SmilesParseParams::default()).unwrap();
        let enol = canonical_graph_key(&enol.topology).unwrap();
        let ketone = canonical_graph_key(&ketone.topology).unwrap();
        assert_ne!(enol, ketone);
        let mut seen = BTreeSet::from([enol.clone()]);
        assert!(!visit_profile_state(&enol, enol.clone(), &mut seen).unwrap());
        assert!(visit_profile_state(&enol, ketone.clone(), &mut seen).unwrap());
        assert!(!visit_profile_state(&ketone, ketone.clone(), &mut seen).unwrap());
        let failure = visit_profile_state(&ketone, enol, &mut seen).unwrap_err();
        assert_eq!(failure.kind, MsReadyFailureKind::TransformationCycle);
        assert_eq!(failure.stage, MsReadyStage::Normalization);
    }
}
