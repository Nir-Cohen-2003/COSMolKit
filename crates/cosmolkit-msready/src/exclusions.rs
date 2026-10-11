use std::sync::OnceLock;

use crate::{MsReadyFailure, MsReadyFailureKind, MsReadyStage};
use cosmolkit_model::{
    AtomId, AtomSpec, BondId, BondSpec, CoordinateBlock, MoleculeProperties, QueryAtom, QueryBond,
    QueryGraph, TopologyBlock,
};
use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

const SOURCE_ROWS: &str = include_str!("data_exclusions_source.json");
const VARIANTS: &str = include_str!("data_exclusion_variants.json");
const NITRATE_TEMPLATES: &str = include_str!("data_nitrate_templates.json");

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExclusionMatch {
    pub source_id: String,
    pub source_field: String,
    pub template_id: Option<String>,
    pub adaptation_id: Option<String>,
    pub comparison_only: bool,
    /// Template atom rows mapped to the original component's atom rows.
    pub atom_mappings: Vec<Vec<usize>>,
    /// Ordinary explicit-H atom row to its original parent atom row.
    pub folded_hydrogen_mappings: Vec<(usize, usize)>,
    pub matched_formal_charge: i32,
    pub matched_formula: Option<String>,
}

#[derive(Clone)]
struct Candidate {
    source_id: String,
    source_field: String,
    template_id: Option<String>,
    adaptation_id: Option<String>,
    comparison_only: bool,
    topology: TopologyBlock,
}

fn catalogue() -> Result<&'static [Candidate], MsReadyFailure> {
    static CATALOGUE: OnceLock<Result<Vec<Candidate>, MsReadyFailure>> = OnceLock::new();
    CATALOGUE
        .get_or_init(build_catalogue)
        .as_deref()
        .map_err(Clone::clone)
}

fn build_catalogue() -> Result<Vec<Candidate>, MsReadyFailure> {
    let mut candidates = Vec::new();
    let mut source_fields = std::collections::BTreeMap::new();
    let rows: serde_json::Value =
        serde_json::from_str(SOURCE_ROWS).map_err(|e| failure(e.to_string()))?;
    let rows = rows
        .as_array()
        .ok_or_else(|| failure("exclusion data is not an array"))?;
    for row in rows {
        let id = row
            .get("ID")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("source row lacks ID"))?;
        let (field, smiles) = selected_source_structure(row)
            .ok_or_else(|| failure(format!("source row {id} lacks selected SMILES field")))?;
        source_fields.insert(id.to_owned(), field);
        let topology = parse_packaged(smiles, &format!("source row {id}"))?;
        let raw_nitrate = id == "83";
        let raw_pf6 = id == "58";
        candidates.push(Candidate {
            source_id: id.to_owned(),
            source_field: field.to_owned(),
            template_id: if raw_nitrate {
                Some("83_raw_h2no3".to_owned())
            } else if raw_pf6 {
                Some("58_raw_pf6".to_owned())
            } else {
                None
            },
            adaptation_id: None,
            comparison_only: raw_nitrate || raw_pf6,
            topology,
        });
    }

    let variants: serde_json::Value =
        serde_json::from_str(VARIANTS).map_err(|e| failure(e.to_string()))?;
    if variants.get("schema_version").and_then(|v| v.as_str()) != Some("1.0.0")
        || variants.get("profile_id").and_then(|v| v.as_str()) != Some("cosmolkit_msready_v1")
    {
        return Err(failure(
            "exclusion variant catalogue profile/version mismatch",
        ));
    }
    for variant in variants
        .get("variants")
        .and_then(|v| v.as_array())
        .ok_or_else(|| failure("exclusion variants missing"))?
    {
        let ids = variant
            .get("source_ids")
            .and_then(|v| v.as_array())
            .ok_or_else(|| failure("variant source_ids missing"))?;
        let smiles = variant
            .get("smiles")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("variant SMILES missing"))?;
        let template_id = variant
            .get("variant_id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("variant id missing"))?;
        let adaptation_id = variant
            .get("adaptation_id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("variant adaptation id missing"))?;
        let source_field = variant
            .get("source_field")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("variant source field missing"))?;
        if !matches!(source_field, "Structure(SMILES)" | "CANONICAL_SMILES") {
            return Err(failure(format!(
                "variant {template_id} has unapproved source field {source_field}"
            )));
        }
        let topology = parse_packaged(smiles, &format!("variant {template_id}"))?;
        for id in ids {
            let id = id
                .as_str()
                .ok_or_else(|| failure("variant source id is not text"))?;
            if source_fields.get(id).copied() != Some(source_field) {
                return Err(failure(format!(
                    "variant {template_id} source field disagrees with source row {id}"
                )));
            }
            candidates.push(Candidate {
                source_id: id.to_owned(),
                source_field: source_field.to_owned(),
                template_id: Some(template_id.to_owned()),
                adaptation_id: Some(adaptation_id.to_owned()),
                comparison_only: variant
                    .get("comparison_only")
                    .and_then(|v| v.as_bool())
                    .unwrap_or(false),
                topology: topology.clone(),
            });
        }
    }

    let nitrate: serde_json::Value =
        serde_json::from_str(NITRATE_TEMPLATES).map_err(|e| failure(e.to_string()))?;
    for template in nitrate
        .get("templates")
        .and_then(|v| v.as_array())
        .ok_or_else(|| failure("nitrate templates missing"))?
    {
        let template_id = template
            .get("template_id")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("nitrate template id missing"))?;
        if template_id == "83_raw_h2no3" {
            continue; // The source row already contributes the verbatim raw template.
        }
        let smiles = template
            .get("smiles")
            .and_then(|v| v.as_str())
            .ok_or_else(|| failure("nitrate template SMILES missing"))?;
        let topology = parse_packaged(smiles, &format!("nitrate template {template_id}"))?;
        let adaptation_id = template
            .get("adaptation_id")
            .and_then(|v| v.as_str())
            .or_else(|| nitrate.get("adaptation_id").and_then(|v| v.as_str()));
        candidates.push(Candidate {
            source_id: template
                .get("source_id")
                .and_then(|v| v.as_str())
                .unwrap_or("83")
                .to_owned(),
            source_field: nitrate
                .get("source_field")
                .and_then(|v| v.as_str())
                .unwrap_or("CANONICAL_SMILES")
                .to_owned(),
            template_id: Some(template_id.to_owned()),
            adaptation_id: adaptation_id.map(str::to_owned),
            comparison_only: template
                .get("comparison_only")
                .and_then(|v| v.as_bool())
                .unwrap_or(false),
            topology,
        });
    }
    Ok(candidates)
}

fn selected_source_structure(row: &serde_json::Value) -> Option<(&'static str, &str)> {
    row.get("Structure(SMILES)")
        .and_then(|v| v.as_str())
        .map(|smiles| ("Structure(SMILES)", smiles))
        .or_else(|| {
            row.get("CANONICAL_SMILES")
                .and_then(|v| v.as_str())
                .map(|smiles| ("CANONICAL_SMILES", smiles))
        })
}

fn parse_packaged(smiles: &str, label: &str) -> Result<TopologyBlock, MsReadyFailure> {
    parse_smiles(smiles, &SmilesParseParams::default())
        .map(|record| record.topology)
        .map_err(|e| failure(format!("invalid packaged {label}: {e}")))
}

/// Match an entire component against the frozen source rows and finite comparison variants.
/// The input graph is borrowed and never normalized or edited.
pub fn match_exact_exclusions(
    component: &TopologyBlock,
) -> Result<Vec<ExclusionMatch>, MsReadyFailure> {
    let folded_component = fold_ordinary_hydrogens(component)?;
    let mut matches = Vec::new();
    for candidate in catalogue()? {
        let folded_candidate = fold_ordinary_hydrogens(&candidate.topology)?;
        let atom_mappings =
            exact_search_mappings(&folded_candidate.topology, &folded_component.topology)?;
        if atom_mappings.is_empty() {
            continue;
        }
        let atom_mappings = atom_mappings
            .iter()
            .map(|mapping| {
                mapping
                    .iter()
                    .map(|&row| folded_component.old_rows[row])
                    .collect()
            })
            .collect();
        let mut folded_hydrogen_mappings: Vec<_> = folded_component.folded_hydrogens.clone();
        folded_hydrogen_mappings.sort_unstable();
        folded_hydrogen_mappings.dedup();
        matches.push(ExclusionMatch {
            source_id: candidate.source_id.clone(),
            source_field: candidate.source_field.clone(),
            template_id: candidate.template_id.clone(),
            adaptation_id: candidate.adaptation_id.clone(),
            comparison_only: candidate.comparison_only,
            atom_mappings,
            folded_hydrogen_mappings,
            matched_formal_charge: component
                .atoms
                .iter()
                .map(|atom| i32::from(atom.formal_charge()))
                .sum(),
            // Exact exclusion-only raw states may have unusual strict valence.
            // Metadata uses an explicit non-strict owner assignment; errors are
            // propagated, never silently replaced by a missing formula.
            matched_formula: Some({
                let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                    component,
                    cosmolkit_core::ValenceModel::RdkitLike,
                    false,
                )
                .map_err(|error| failure(format!("matched exclusion formula valence: {error}")))?;
                cosmolkit_descriptors::molecular_formula_with_valence(
                    component,
                    true,
                    false,
                    Some(&assignment),
                )
                .map_err(|error| failure(format!("matched exclusion formula: {error}")))?
            }),
        });
    }
    Ok(matches)
}

/// The raw-valence exception and exclusion comparison share the same template
/// and complete atom/bond/H/isotope predicates, including ordinary-H folding.
pub(crate) fn is_exact_raw_nitrate(topology: &TopologyBlock) -> Result<bool, MsReadyFailure> {
    let raw = catalogue()?
        .iter()
        .find(|entry| entry.template_id.as_deref() == Some("83_raw_h2no3"))
        .ok_or_else(|| failure("raw nitrate template missing from catalogue"))?;
    let target = fold_ordinary_hydrogens(topology)?;
    let template = fold_ordinary_hydrogens(&raw.topology)?;
    Ok(!exact_search_mappings(&template.topology, &target.topology)?.is_empty())
}

pub fn source_exclusion_count() -> Result<usize, MsReadyFailure> {
    let rows: serde_json::Value =
        serde_json::from_str(SOURCE_ROWS).map_err(|e| failure(e.to_string()))?;
    Ok(rows
        .as_array()
        .ok_or_else(|| failure("exclusion data is not an array"))?
        .len())
}

struct FoldedTopology {
    topology: TopologyBlock,
    /// Folded atom row -> original graph atom row.
    old_rows: Vec<usize>,
    /// Original ordinary explicit-H row -> original parent row.
    folded_hydrogens: Vec<(usize, usize)>,
}

fn fold_ordinary_hydrogens(topology: &TopologyBlock) -> Result<FoldedTopology, MsReadyFailure> {
    let ordinary_hydrogens: Vec<_> = topology
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(row, atom)| {
            (atom.element() == cosmolkit_types::Element::H && atom.isotope().is_none())
                .then_some(row)
        })
        .collect();
    if ordinary_hydrogens.is_empty() {
        return Ok(FoldedTopology {
            topology: topology.clone(),
            old_rows: (0..topology.atoms.len()).collect(),
            folded_hydrogens: Vec::new(),
        });
    }
    let mut params = cosmolkit_core::RemoveHsParams::default();
    params.remove_isotopes = false;
    params.remove_and_track_isotopes = false;
    params.remove_with_wedged_bond = false;
    params.update_explicit_count = true;
    params.sanitize = false;
    let result = cosmolkit_core::remove_hydrogens_with_params(
        topology.clone(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &params,
    )
    .map_err(|e| failure(format!("core ordinary-H folding failed: {e}")))?;
    let old_rows: Vec<_> = result
        .mapping
        .atoms
        .new_to_old
        .iter()
        .map(|id| {
            id.map(|id| id.index())
                .ok_or_else(|| failure("core hydrogen mapping lacks original atom row"))
        })
        .collect::<Result<_, _>>()?;
    let mut folded_hydrogens = Vec::new();
    for hydrogen in ordinary_hydrogens {
        if result
            .mapping
            .atoms
            .old_to_new
            .get(hydrogen)
            .is_some_and(Option::is_none)
        {
            let neighbors = topology.adjacency.neighbors_of(hydrogen);
            if neighbors.len() == 1
                && topology.atoms[neighbors[0].atom_index].element() != cosmolkit_types::Element::H
            {
                folded_hydrogens.push((hydrogen, neighbors[0].atom_index));
            }
        }
    }
    Ok(FoldedTopology {
        topology: result.topology,
        old_rows,
        folded_hydrogens,
    })
}

/// The SEARCH owner supplies whole-graph isomorphism/mappings; exact source-owned
/// atom and bond fields are then checked against every complete bijection.
pub(crate) fn exact_search_mappings(
    query_topology: &TopologyBlock,
    target: &TopologyBlock,
) -> Result<Vec<Vec<usize>>, MsReadyFailure> {
    if query_topology.atoms.len() != target.atoms.len()
        || query_topology.bonds.len() != target.bonds.len()
    {
        return Ok(Vec::new());
    }
    let atoms = query_topology
        .atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| QueryAtom::new(AtomId::new(i), AtomSpec::new(atom.element())))
        .collect();
    let bonds = query_topology
        .bonds
        .iter()
        .enumerate()
        .map(|(i, bond)| {
            let spec = BondSpec::new(bond.begin(), bond.end(), bond.order())
                .with_aromatic(bond.is_aromatic());
            QueryBond::new(BondId::new(i), spec)
        })
        .collect();
    let query = QueryGraph::from_parts(atoms, bonds, [], vec![], vec![], vec![])
        .map_err(|e| failure(format!("could not build exact exclusion query: {e}")))?;
    let search_results = cosmolkit_search::match_query(&query, target)
        .map_err(|e| failure(format!("owner whole-graph search failed: {e}")))?;
    let query_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        query_topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .ok();
    let target_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        target,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .ok();
    let mut mappings = Vec::new();
    for result in search_results {
        let mapping = result.atom_mapping;
        if mapping.len() != query_topology.atoms.len()
            || result.bond_mapping.len() != query_topology.bonds.len()
        {
            continue;
        }
        let atom_match = query_topology.atoms.iter().enumerate().all(|(left, a)| {
            let right = mapping[left];
            let Some(b) = target.atoms.get(right) else {
                return false;
            };
            a.formal_charge() == b.formal_charge()
                && a.isotope() == b.isotope()
                && a.radical_electrons() == b.radical_electrons()
                && ordinary_hydrogen_count_matches(
                    i32::from(a.explicit_hydrogens()),
                    query_valence
                        .as_ref()
                        .and_then(|v| v.implicit_hydrogens.get(left))
                        .copied(),
                    i32::from(b.explicit_hydrogens()),
                    target_valence
                        .as_ref()
                        .and_then(|v| v.implicit_hydrogens.get(right))
                        .copied(),
                )
                && a.tracked_isotopic_hydrogens() == b.tracked_isotopic_hydrogens()
                && a.is_aromatic() == b.is_aromatic()
        });
        if !atom_match {
            continue;
        }
        let bond_match = query_topology.bonds.iter().enumerate().all(|(left, a)| {
            let right = result.bond_mapping[left];
            let Some(b) = target.bonds.get(right) else {
                return false;
            };
            let begin = mapping[a.begin().index()];
            let end = mapping[a.end().index()];
            let endpoints = (b.begin().index() == begin && b.end().index() == end)
                || (b.begin().index() == end && b.end().index() == begin);
            endpoints && a.order() == b.order() && a.is_aromatic() == b.is_aromatic()
        });
        if bond_match {
            mappings.push(mapping);
        }
    }
    mappings.sort();
    mappings.dedup();
    Ok(mappings)
}

fn ordinary_hydrogen_count_matches(
    left_explicit: i32,
    left_implicit: Option<i32>,
    right_explicit: i32,
    right_implicit: Option<i32>,
) -> bool {
    match (left_implicit, right_implicit) {
        (Some(left), Some(right)) => left_explicit + left == right_explicit + right,
        (None, None) => left_explicit == right_explicit,
        _ => false,
    }
}

fn failure(context: impl Into<String>) -> MsReadyFailure {
    MsReadyFailure {
        evidence: Default::default(),
        kind: MsReadyFailureKind::OwnerError,
        stage: MsReadyStage::Validation,
        context: context.into(),
        original: None,
    }
}

#[cfg(test)]
mod exclusion_tests {
    use super::*;

    #[test]
    fn every_pinned_source_row_matches_as_a_whole_component() {
        let rows: serde_json::Value = serde_json::from_str(SOURCE_ROWS).unwrap();
        for row in rows.as_array().unwrap() {
            let id = row["ID"].as_str().unwrap();
            let (_, smiles) = selected_source_structure(row).unwrap();
            let graph = parse_smiles(smiles, &SmilesParseParams::default())
                .unwrap_or_else(|e| panic!("source {id}: {e}"))
                .topology;
            let matched = match_exact_exclusions(&graph).unwrap();
            assert!(
                matched.iter().any(|matched| matched.source_id == id),
                "source {id}"
            );
            if id == "58" {
                let pf6 = matched
                    .iter()
                    .find(|matched| matched.template_id.as_deref() == Some("58_raw_pf6"))
                    .unwrap();
                assert!(pf6.comparison_only);
                assert_eq!(pf6.source_field, "CANONICAL_SMILES");
            }
        }
    }

    #[test]
    fn every_finite_variant_is_packaged_and_self_matching() {
        let data: serde_json::Value = serde_json::from_str(VARIANTS).unwrap();
        for variant in data["variants"].as_array().unwrap() {
            let id = variant["variant_id"].as_str().unwrap();
            let smiles = variant["smiles"].as_str().unwrap();
            let graph = parse_smiles(smiles, &SmilesParseParams::default())
                .unwrap_or_else(|e| panic!("variant {id}: {e}"))
                .topology;
            let matched = match_exact_exclusions(&graph).unwrap();
            let matched = matched
                .iter()
                .find(|matched| matched.template_id.as_deref() == Some(id))
                .unwrap_or_else(|| panic!("variant {id}"));
            if id == "58_hexafluorophosphate_anion" {
                assert!(matched.comparison_only);
            }
        }
    }

    #[test]
    fn nitrate_templates_remain_comparison_only_and_exact() {
        assert_eq!(source_exclusion_count().unwrap(), 32);
        let cases = [
            ("ON(=O)O", "83_raw_h2no3", "H2NO3", 0),
            ("[OH][N+](=O)[O-]", "83_nitric_acid", "HNO3", 0),
            ("[O-][N+](=O)[O-]", "83_nitrate", "NO3-", -1),
        ];
        for (smiles, template, formula, charge) in cases {
            let graph = parse_smiles(smiles, &SmilesParseParams::default())
                .unwrap()
                .topology;
            let matched = match_exact_exclusions(&graph)
                .unwrap()
                .into_iter()
                .find(|matched| matched.template_id.as_deref() == Some(template))
                .unwrap();
            assert!(matched.comparison_only);
            assert_eq!(matched.source_id, "83");
            assert_eq!(matched.source_field, "CANONICAL_SMILES");
            assert_eq!(matched.matched_formal_charge, charge);
            assert_eq!(matched.matched_formula.as_deref(), Some(formula));
            assert!(!matched.atom_mappings.is_empty());
            if template == "83_raw_h2no3" {
                assert!(matched.adaptation_id.is_none());
            } else {
                assert_eq!(
                    matched.adaptation_id.as_deref(),
                    Some("nitrate_comparison_charge_adaptation_v1")
                );
            }
        }
        for smiles in ["O=N(O)O", "[H]ON(=O)O", "[H]O[N+](=O)[O-]"] {
            let graph = parse_smiles(smiles, &SmilesParseParams::default())
                .unwrap()
                .topology;
            assert!(
                match_exact_exclusions(&graph)
                    .unwrap()
                    .iter()
                    .any(|m| m.comparison_only),
                "{smiles}"
            );
        }
        let mixed = parse_smiles("CC.O[N](=O)O", &SmilesParseParams::default()).unwrap();
        let original = mixed.topology.clone();
        let parts = crate::components::split_mapped_components(
            &mixed.topology,
            &mixed.coordinates,
            &mixed.properties,
        )
        .unwrap();
        assert!(parts.iter().any(|part| {
            match_exact_exclusions(&part.topology)
                .unwrap()
                .iter()
                .any(|m| m.template_id.as_deref() == Some("83_raw_h2no3"))
        }));
        assert_eq!(mixed.topology, original);
    }

    #[test]
    fn symmetry_order_and_ordinary_explicit_h_are_supported_without_isotope_folding() {
        for smiles in ["O=C(O)O", "OC(=O)O", "[H]OC(=O)O"] {
            let graph = parse_smiles(smiles, &SmilesParseParams::default())
                .unwrap()
                .topology;
            let result = match_exact_exclusions(&graph).unwrap();
            assert!(
                result.iter().any(|matched| matched.source_id == "29"),
                "{smiles}: {result:?}"
            );
            if smiles.starts_with("[H]") {
                assert_eq!(result[0].folded_hydrogen_mappings.len(), 1);
                assert!(result[0].atom_mappings[0].iter().all(|row| *row != 0));
            }
        }
        let isotope = parse_smiles("[2H]OC(=O)O", &SmilesParseParams::default())
            .unwrap()
            .topology;
        assert!(match_exact_exclusions(&isotope).unwrap().is_empty());
    }

    #[test]
    fn duplicate_source_ids_and_core_and_isotope_negatives_are_preserved() {
        let carbonate = parse_smiles("O=C(O)O", &SmilesParseParams::default())
            .unwrap()
            .topology;
        let carbonate_matches = match_exact_exclusions(&carbonate).unwrap();
        let ids: Vec<_> = carbonate_matches
            .iter()
            .map(|m| m.source_id.as_str())
            .collect();
        assert!(ids.contains(&"29") && ids.contains(&"155"));
        for smiles in [
            "[18O-][N+](=O)[O-]",
            "O[N+](=O)[18O-]",
            "[18O]N(=O)O",
            "O=N[O-]",
            "CO[N+](=O)[O-]",
            "C[N+](=O)[O-]",
            "ON(=O)[O-]",
            "O=C(O)CO",
            "CCCl",
            "CP(=O)(O)O",
            "CS(=O)(=O)O",
        ] {
            let graph = parse_smiles(smiles, &SmilesParseParams::default())
                .unwrap()
                .topology;
            assert!(
                match_exact_exclusions(&graph).unwrap().is_empty(),
                "{smiles}"
            );
        }
    }
}
