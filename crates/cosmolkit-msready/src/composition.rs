use std::collections::BTreeMap;

use cosmolkit_core::{ValenceAssignment, total_hydrogen_count_from_validated};
use cosmolkit_descriptors::{molecular_formula_with_valence, neutral_atom_exact_mass_with_valence};
use cosmolkit_model::TopologyBlock;

use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};

/// Counts include inferred/bracket ordinary H and explicit isotope-bearing atoms.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ProductComposition {
    pub counts: BTreeMap<(u8, Option<u16>), u32>,
    pub formula: String,
    pub formal_charge: i32,
    pub neutral_atom_mass: f64,
}

/// Requires materialized tracked isotope H and the final owner assignment.
pub(crate) fn product_composition(
    topology: &TopologyBlock,
    assignment: &ValenceAssignment,
) -> Result<ProductComposition, MsReadyFailure> {
    // BEGIN MSREADY REFERENCE §10
    // MSReady❗❗: Store the final elemental formula, isotope composition, formal charge and neutral-atom formula mass separately.
    // END MSREADY REFERENCE
    // Native behavior: isotope-aware formula and neutral-atom mass delegate to
    // descriptors; no electron constant/mass table or adduct calculation here.
    // Native complexity: one O(V) count pass with O(log K) ordered-map updates;
    // descriptor kernels borrow the same assignment, without recalculating it.
    let mut counts = BTreeMap::<(u8, Option<u16>), u32>::new();
    let mut formal_charge = 0_i32;
    let mut add = |key, amount: u32| -> Result<(), MsReadyFailure> {
        let count = counts.entry(key).or_default();
        *count = count.checked_add(amount).ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                "product composition count overflow",
            )
        })?;
        Ok(())
    };
    for atom in &topology.atoms {
        if !atom.tracked_isotopic_hydrogens().is_empty() {
            return Err(failure(
                Kind::InvariantViolation,
                "composition requires materialized isotope H",
            ));
        }
        add((atom.atomic_number(), atom.isotope()), 1)?;
        let hydrogens = total_hydrogen_count_from_validated(topology, assignment, atom.id(), false)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        if hydrogens != 0 {
            add((1, None), hydrogens)?;
        }
        formal_charge = formal_charge
            .checked_add(i32::from(atom.formal_charge()))
            .ok_or_else(|| failure(Kind::InvariantViolation, "product formal charge overflow"))?;
    }
    let formula = molecular_formula_with_valence(topology, true, false, Some(assignment))
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let neutral_atom_mass = neutral_atom_exact_mass_with_valence(topology, Some(assignment))
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    if !neutral_atom_mass.is_finite() {
        return Err(failure(Kind::InvariantViolation, "non-finite product mass"));
    }
    Ok(ProductComposition {
        counts,
        formula,
        formal_charge,
        neutral_atom_mass,
    })
}

fn failure(kind: Kind, context: impl Into<String>) -> MsReadyFailure {
    MsReadyFailure {
        evidence: Default::default(),
        kind,
        stage: Stage::Finalization,
        context: context.into(),
        original: None,
    }
}
