use cosmolkit_core::{SanitizeOperations, SanitizeParams, sanitize_topology};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock, TopologyMapping};
use cosmolkit_smiles::{SmilesRecord, SmilesWriteParams, write_smiles_with_params};
use cosmolkit_types::{BondDirection, BondStereo, ChiralTag};

use crate::composition::{ProductComposition, product_composition};
use crate::hydrogens::{
    HydrogenDisposition, fold_ordinary_explicit_hydrogens, materialize_tracked_isotopic_hydrogens,
};
use crate::identify::{StandardIdentity, standard_identity};
use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};

// Persisted/raw as well as computed carriers must be removed. These are the
// source-named stereo/perception/drawing properties used by core, SMILES and IO.
const SEARCH_CARRIER_PROPERTIES: &[&str] = &[
    "_CIPCode",
    "_CIPRank",
    "_CIPNeighborOrder",
    "_ChiralAtomRank",
    "_ChiralityPossible",
    "_NonExplicit3DChirality",
    "_UnknownStereo",
    "_ringStereoAtoms",
    "_ringStereochemCand",
    "_brokenChirality",
    "_StereochemDone",
    "_potentialStereo",
    "_needsDetectAtomStereo",
    "_needsDetectBondStereo",
    "_MolFileChiralFlag",
    "_MolFileStereoCare",
    "_MolFileBondStereo",
    "_MolFileBondCfg",
    "molStereoCare",
    "molAtomMapNumber",
    "molFileAlias",
    "molFileValue",
    "_MolFile_PXA",
    "_MolFileInfo",
    "_MolFileInfoLine",
    "_MolFileComments",
    "_chiralPermutation",
    "_STEREO",
    "_atropcw",
    "_atropccw",
    "_stereoGroup",
    "stereoGroupAtom",
];

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct FinalizedComponent {
    pub record: SmilesRecord,
    pub canonical_smiles: String,
    pub identity: StandardIdentity,
    pub composition: ProductComposition,
    pub original_to_final: TopologyMapping,
    pub hydrogen_dispositions: Vec<HydrogenDisposition>,
}

/// Finalize one already-normalized retained component; no exclusion/proton rules.
pub(crate) fn finalize_component(
    original: &TopologyBlock,
    current: &SmilesRecord,
    original_to_current: &TopologyMapping,
) -> Result<FinalizedComponent, MsReadyFailure> {
    // BEGIN MSREADY REFERENCE §9
    // MSReady❗❗: Clear tetrahedral, double-bond, allene/cumulene and supported enhanced or axial stereo annotations.
    // MSReady❗❗: Clear bond directions and molecule-level stereo groups, and prevent coordinates from causing stereo to be perceived again during serialization.
    // END MSREADY REFERENCE
    // Native behavior: physically clear derived carriers; isotope preservation
    // overrides archived isotope-removal settings. Owner H mapping, selected
    // aromaticity/H/valence, SMILES, descriptors and real Standard InChI are
    // reused, never a generic cleanup or identifier graph roundtrip.
    // Native complexity: fixed linear carrier passes and detached owner clones;
    // aromatic/ring, canonical-writer and identifier costs remain owner-defined.
    // No historical execution/performance equivalence is asserted.
    crate::validate_syntax(
        &current.topology,
        &current.properties,
        &crate::MsReadyParams::default(),
    )?;
    let mut graph = current.topology.clone();
    let mut properties = current.properties.clone();
    clear_search_carriers(&mut graph, &mut properties)?;
    let before_graph = materialize_tracked_isotopic_hydrogens(&graph)?;
    let before_valence = crate::validate::strict_owner_valence(&before_graph)?;
    let before = product_composition(&before_graph, &before_valence)?;

    let folded = fold_ordinary_explicit_hydrogens(original, &graph, original_to_current)?;
    graph = materialize_tracked_isotopic_hydrogens(&folded.topology)?;
    // Materialization appends only new isotope H atoms/bonds, with no source
    // atom identity. Original explicit D/T and every heavy atom retain their map.
    let mut original_to_final = folded.topology_mapping;
    original_to_final
        .atoms
        .new_to_old
        .resize(graph.atoms.len(), None);
    original_to_final
        .bonds
        .new_to_old
        .resize(graph.bonds.len(), None);
    original_to_final
        .validate_for_counts(
            original.atoms.len(),
            graph.atoms.len(),
            original.bonds.len(),
            graph.bonds.len(),
        )
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    let operations = SanitizeOperations::PROPERTIES
        | SanitizeOperations::SYMM_RINGS
        | SanitizeOperations::KEKULIZE
        | SanitizeOperations::SET_AROMATICITY
        | SanitizeOperations::SET_CONJUGATION
        | SanitizeOperations::SET_HYBRIDIZATION
        | SanitizeOperations::ADJUST_HS;
    let sanitized = sanitize_topology(&graph, &SanitizeParams { operations }).map_err(|error| {
        let kind = if matches!(error, cosmolkit_core::SanitizeError::Properties { .. }) {
            Kind::InvalidValence
        } else {
            Kind::OwnerError
        };
        failure(kind, error.to_string())
    })?;
    graph = sanitized.topology;
    clear_search_carriers(&mut graph, &mut properties)?;
    // No radical perception stage: existing electron counts are immutable.
    for (before_atom, after_atom) in folded.topology.atoms.iter().zip(&graph.atoms) {
        if before_atom.radical_electrons() != after_atom.radical_electrons() {
            return Err(failure(
                Kind::InvariantViolation,
                "finalization changed radical electrons",
            ));
        }
    }
    graph
        .validate()
        .map_err(|error| failure(Kind::InvariantViolation, error.to_string()))?;
    crate::retained_metal_states(&graph)?;
    let assignment = crate::validate::strict_owner_valence(&graph)?;
    let composition = product_composition(&graph, &assignment)?;
    if before.counts != composition.counts
        || before.formal_charge != composition.formal_charge
        || (before.neutral_atom_mass - composition.neutral_atom_mass).abs() > 1e-9
    {
        return Err(failure(
            Kind::InvariantViolation,
            "finalization changed isotope composition, charge or mass",
        ));
    }
    let record = SmilesRecord {
        topology: graph,
        coordinates: CoordinateBlock::default(),
        properties,
    };
    let notation = write_smiles_with_params(&record, &writer_params())
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let canonical_smiles = String::from_utf8(notation.as_bytes().to_vec())
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    let identity = standard_identity(&record.topology)?;
    if canonical_smiles.contains(['@', '/', '\\'])
        || identity
            .inchi
            .split('/')
            .any(|layer| matches!(layer.as_bytes().first(), Some(b'b' | b't' | b'm' | b's')))
    {
        return Err(failure(
            Kind::InvariantViolation,
            "stereo reintroduced by writer/identifier",
        ));
    }
    Ok(FinalizedComponent {
        record,
        canonical_smiles,
        identity,
        composition,
        original_to_final,
        hydrogen_dispositions: folded.dispositions,
    })
}

pub(crate) fn writer_params() -> SmilesWriteParams {
    SmilesWriteParams {
        do_isomeric_smiles: true,
        do_kekule: false,
        canonical: true,
        clean_stereo: false,
        rooted_at_atom: None,
        all_bonds_explicit: false,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: true,
    }
}

fn clear_search_carriers(
    graph: &mut TopologyBlock,
    properties: &mut MoleculeProperties,
) -> Result<(), MsReadyFailure> {
    for atom in &mut graph.atoms {
        atom.set_chiral_tag(ChiralTag::Unspecified);
        atom.set_chiral_permutation(None);
        atom.set_unknown_stereo(false);
        atom.set_mol_parity(None);
        atom.set_mol_inversion_flag(None);
        atom.set_atom_map(None);
        atom.clear_computed_props()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        for key in SEARCH_CARRIER_PROPERTIES {
            atom.clear_prop(key)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        }
    }
    for bond in &mut graph.bonds {
        bond.set_stereo(BondStereo::None)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        bond.set_stereo_atoms(None);
        bond.set_direction(BondDirection::None);
        bond.set_unknown_stereo(false);
        bond.clear_computed_props()
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        for key in SEARCH_CARRIER_PROPERTIES {
            bond.clear_prop(key)
                .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
        }
    }
    graph.stereo_groups.clear();
    properties
        .clear_computed_props()
        .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    for key in SEARCH_CARRIER_PROPERTIES {
        properties
            .clear_prop(key)
            .map_err(|error| failure(Kind::OwnerError, error.to_string()))?;
    }
    Ok(())
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

#[cfg(test)]
mod finalize_tests {
    use super::*;
    use cosmolkit_model::{AtomId, Conformer2D, Conformer3D, StereoGroup, StereoGroupKind};
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};
    use cosmolkit_types::BondOrder;

    fn parse(text: &str) -> SmilesRecord {
        parse_smiles(text, &SmilesParseParams::default()).unwrap()
    }
    fn finalize(record: &SmilesRecord) -> FinalizedComponent {
        let original = record.clone();
        let product = finalize_component(
            &record.topology,
            record,
            &TopologyMapping::identity(record.topology.atoms.len(), record.topology.bonds.len()),
        )
        .unwrap();
        assert_eq!(
            *record, original,
            "original topology/coordinates/properties mutated"
        );
        assert_clean(&product);
        product
    }
    fn assert_clean(product: &FinalizedComponent) {
        let graph = &product.record.topology;
        assert_eq!(product.record.coordinates, CoordinateBlock::default());
        assert!(graph.stereo_groups.is_empty());
        for atom in &graph.atoms {
            assert_eq!(atom.chiral_tag(), ChiralTag::Unspecified);
            assert_eq!(atom.chiral_permutation(), None);
            assert!(!atom.unknown_stereo());
            assert_eq!(atom.mol_parity(), None);
            assert_eq!(atom.mol_inversion_flag(), None);
            assert_eq!(atom.atom_map(), None);
            assert!(atom.tracked_isotopic_hydrogens().is_empty());
            for key in SEARCH_CARRIER_PROPERTIES {
                assert!(atom.prop(key).is_none(), "atom {key}");
            }
        }
        for bond in &graph.bonds {
            assert_eq!(bond.stereo(), BondStereo::None);
            assert!(bond.stereo_atom_references().is_empty());
            assert_eq!(bond.direction(), BondDirection::None);
            assert!(!bond.unknown_stereo());
            for key in SEARCH_CARRIER_PROPERTIES {
                assert!(bond.prop(key).is_none(), "bond {key}");
            }
        }
        for key in SEARCH_CARRIER_PROPERTIES {
            assert!(
                product.record.properties.prop(key).is_none(),
                "molecule {key}"
            );
        }
        assert!(!product.canonical_smiles.contains(['@', '/', '\\']));
        assert!(product.identity.inchi.starts_with("InChI=1S/"));
        assert_eq!(product.identity.full_key.len(), 27);
        assert!(
            !product
                .identity
                .inchi
                .split('/')
                .any(|layer| matches!(layer.as_bytes().first(), Some(b'b' | b't' | b'm' | b's')))
        );
        let reparsed = parse(&product.canonical_smiles);
        let assignment = crate::validate::strict_owner_valence(&reparsed.topology).unwrap();
        let composition = product_composition(&reparsed.topology, &assignment).unwrap();
        assert_eq!(composition.counts, product.composition.counts);
        assert_eq!(composition.formula, product.composition.formula);
        assert_eq!(composition.formal_charge, product.composition.formal_charge);
        assert!(
            (composition.neutral_atom_mass - product.composition.neutral_atom_mass).abs() < 1e-9
        );
        assert_eq!(
            standard_identity(&reparsed.topology).unwrap().full_key,
            product.identity.full_key
        );
    }
    fn assert_same(left: &FinalizedComponent, right: &FinalizedComponent) {
        assert_eq!(left.canonical_smiles, right.canonical_smiles);
        assert_eq!(left.identity.inchi, right.identity.inchi);
        assert_eq!(left.identity.full_key, right.identity.full_key);
        assert_eq!(left.composition, right.composition);
    }

    #[test]
    fn tetrahedral_enantiomers_have_one_real_stereo_free_identity() {
        let left = parse("C[C@H](O)C(=O)O");
        assert!(
            left.topology
                .atoms
                .iter()
                .any(|a| a.chiral_tag() != ChiralTag::Unspecified)
        );
        let left = finalize(&left);
        let right = finalize(&parse("C[C@@H](O)C(=O)O"));
        assert_same(&left, &right);
        assert_eq!(left.identity.full_key, "JVTAAEKCZFNVCJ-UHFFFAOYSA-N");
        assert_eq!(left.composition.formula, "C3H6O3");
        assert_eq!(left.record.topology.atoms.len(), 6);
    }

    #[test]
    fn ez_directions_and_reference_atoms_are_physically_removed() {
        let left = parse("C/C=C/C");
        assert!(
            left.topology
                .bonds
                .iter()
                .any(|b| b.direction() != BondDirection::None || b.stereo() != BondStereo::None)
        );
        let left = finalize(&left);
        assert_same(&left, &finalize(&parse("C/C=C\\C")));
        assert_eq!(left.composition.formula, "C4H8");
        assert_eq!(
            left.record
                .topology
                .bonds
                .iter()
                .filter(|b| b.order() == BondOrder::Double)
                .count(),
            1
        );
    }

    #[test]
    fn enhanced_groups_maps_drawing_and_raw_perception_carriers_are_cleared() {
        let mut record = parse("C[C@H](O)C(=O)O");
        for kind in [
            StereoGroupKind::Absolute,
            StereoGroupKind::Or,
            StereoGroupKind::And,
        ] {
            record.topology.stereo_groups.push(StereoGroup::new(
                kind,
                vec![AtomId::new(1)],
                Vec::new(),
            ));
        }
        record.topology.atoms[1].set_atom_map(Some(42));
        record.topology.atoms[1].set_mol_parity(Some(1));
        record.topology.atoms[1].set_mol_inversion_flag(Some(2));
        record.topology.atoms[1].set_unknown_stereo(true);
        record.topology.bonds[0].set_direction(BondDirection::BeginWedge);
        record.topology.bonds[0].set_unknown_stereo(true);
        // Raw/non-computed carriers are deliberately malformed: removal, not
        // perception or fallback, must prevent their re-use by the owners.
        for key in SEARCH_CARRIER_PROPERTIES {
            record.topology.atoms[1]
                .set_prop(*key, "raw drawing/stereo")
                .unwrap();
            record.topology.bonds[0]
                .set_prop(*key, "raw drawing/stereo")
                .unwrap();
            record
                .properties
                .set_prop(*key, "raw drawing/stereo")
                .unwrap();
        }
        record.topology.atoms[1]
            .set_prop("application_note", "preserved")
            .unwrap();
        record
            .properties
            .set_prop("record_note", "preserved")
            .unwrap();
        let result = finalize(&record);
        assert_same(&result, &finalize(&parse("CC(O)C(=O)O")));
        assert!(
            result.record.topology.atoms[1]
                .prop("application_note")
                .is_some()
        );
        assert!(result.record.properties.prop("record_note").is_some());
    }

    #[test]
    fn axial_and_allene_permutations_do_not_survive_in_graph_or_identifiers() {
        for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
            let mut record = parse("Cc1cccc(C)c1-c1c(C)cccc1C");
            let bond = record
                .topology
                .bonds
                .iter()
                .position(|b| {
                    !b.is_aromatic()
                        && record.topology.atoms[b.begin().index()].is_aromatic()
                        && record.topology.atoms[b.end().index()].is_aromatic()
                })
                .unwrap();
            record.topology.bonds[bond].set_stereo(stereo).unwrap();
            record.topology.bonds[bond]
                .set_source_stereo_atom_references(vec![AtomId::new(1), AtomId::new(8)]);
            let product = finalize(&record);
            assert_same(&product, &finalize(&parse("Cc1cccc(C)c1-c1c(C)cccc1C")));
        }
        let mut allene = parse("CC=C=CC");
        allene.topology.atoms[2].set_chiral_tag(ChiralTag::Allene);
        allene.topology.atoms[2].set_chiral_permutation(Some(1));
        allene.topology.atoms[2]
            .set_prop("_chiralPermutation", 1_u32)
            .unwrap();
        assert_same(&finalize(&allene), &finalize(&parse("CC=C=CC")));
    }

    #[test]
    fn actual_coordinate_induced_chirality_is_removed_without_touching_original_conformers() {
        let mut record = parse("FC(Cl)(Br)I");
        record.coordinates.conformers_3d.push(Conformer3D::new(
            0,
            vec![
                [1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0],
                [-1.0, -1.0, 1.0],
                [-1.0, 1.0, -1.0],
                [1.0, -1.0, -1.0],
            ],
            true,
        ));
        record
            .coordinates
            .conformers_2d
            .push(Conformer2D::new(1, vec![[1.0, 2.0]; 5]));
        let valence = crate::validate::strict_owner_valence(&record.topology).unwrap();
        record.topology = cosmolkit_core::assign_chiral_tags_from_structure(
            &record.topology,
            &record.coordinates,
            &valence,
            &cosmolkit_core::StructureTagParams::default(),
        )
        .unwrap()
        .topology;
        assert_ne!(
            record.topology.atoms[1].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert!(
            record.topology.atoms[1]
                .prop("_NonExplicit3DChirality")
                .is_some()
        );
        assert_same(&finalize(&record), &finalize(&parse("FC(Cl)(Br)I")));
    }

    #[test]
    fn ordinary_h_fold_has_owner_mapping_and_preserves_explicit_dt() {
        let record = parse("[H]OC([2H])([3H])C");
        let h_row = record
            .topology
            .atoms
            .iter()
            .position(|a| a.atomic_number() == 1 && a.isotope().is_none())
            .unwrap();
        let product = finalize(&record);
        assert_eq!(
            product.record.topology.atoms.len(),
            record.topology.atoms.len() - 1
        );
        assert_eq!(product.original_to_final.atoms.old_to_new[h_row], None);
        assert!(
            product
                .hydrogen_dispositions
                .iter()
                .any(|h| h.input_atom_row == h_row
                    && matches!(
                        h.kind,
                        crate::hydrogens::HydrogenDispositionKind::FoldedInto { .. }
                    ))
        );
        for (row, atom) in record
            .topology
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number() == 1 && a.isotope().is_some())
        {
            let mapped = product.original_to_final.atoms.old_to_new[row].unwrap();
            assert_eq!(
                product.record.topology.atoms[mapped.index()].isotope(),
                atom.isotope()
            );
        }
        assert_eq!(product.composition.counts[&(1, Some(2))], 1);
        assert_eq!(product.composition.counts[&(1, Some(3))], 1);
        assert!(product.canonical_smiles.contains("[2H]"));
        assert!(product.canonical_smiles.contains("[3H]"));
        assert_eq!(product.composition.counts[&(1, None)], 4);
    }

    #[test]
    fn tracked_dt_materialization_matches_explicit_labels_and_locations() {
        for (text, row, explicit) in [
            ("CCO", 2, "CCO[2H]"),
            ("CC[OH]", 1, "CC([2H])O"),
            ("CC[O]", 2, "CCO[2H]"),
        ] {
            let mut tracked = parse(text);
            tracked.topology.atoms[row].set_tracked_isotopic_hydrogens(vec![2]);
            // The no-implicit [O] carrier is a bound labelled hydroxyl, not an
            // independent radical after materialization.
            if text == "CC[O]" {
                tracked.topology.atoms[row].set_radical_electrons(0);
            }
            let result = finalize(&tracked);
            assert_same(&result, &finalize(&parse(explicit)));
            assert_eq!(
                result.record.topology.atoms.len(),
                tracked.topology.atoms.len() + 1
            );
            assert_eq!(
                result.original_to_final.atoms.new_to_old.last(),
                Some(&None)
            );
            tracked.topology.atoms[row].set_tracked_isotopic_hydrogens(vec![3]);
            let tritium = finalize(&tracked);
            assert_ne!(result.identity.full_key, tritium.identity.full_key);
            assert!(tritium.composition.neutral_atom_mass > result.composition.neutral_atom_mass);
        }
        let mut tracked = parse("CCO");
        tracked.topology.atoms[1].set_tracked_isotopic_hydrogens(vec![2, 3]);
        assert_same(&finalize(&tracked), &finalize(&parse("CC([2H])([3H])O")));
    }

    #[test]
    fn isotope_induced_tetrahedral_stereo_clears_without_erasing_labels() {
        let left = finalize(&parse("[2H][C@H](O)C"));
        let right = finalize(&parse("[2H][C@@H](O)C"));
        assert_same(&left, &right);
        assert_ne!(
            left.identity.full_key,
            finalize(&parse("CCO")).identity.full_key
        );
        assert!(left.identity.inchi.contains("/i"));
        let carbon = finalize(&parse("C[13CH2]O"));
        assert!(carbon.canonical_smiles.contains("13C"));
        assert!(carbon.composition.formula.contains("[13C]"));
        assert_eq!(carbon.composition.counts[&(6, Some(13))], 1);
        assert_ne!(
            carbon.identity.full_key,
            finalize(&parse("CCO")).identity.full_key
        );
        assert!((carbon.composition.neutral_atom_mass - 47.045219652).abs() < 1e-8);
    }

    #[test]
    fn aromatic_and_kekule_graphs_converge_including_nh() {
        for (aromatic, kekule) in [("c1ccccc1", "C1=CC=CC=C1"), ("c1cc[nH]c1", "C1=CNC=C1")] {
            let left = finalize(&parse(aromatic));
            let right = finalize(&parse(kekule));
            assert_same(&left, &right);
            assert!(right.record.topology.atoms.iter().all(|a| a.is_aromatic()));
            assert!(right.record.topology.bonds.iter().all(|b| b.is_aromatic()));
        }
    }

    #[test]
    fn ions_and_isotope_labels_have_neutral_atom_not_ionic_mass() {
        for (text, formula, mass) in [
            ("CCO", "C2H6O", 46.041864812),
            ("[NH4+]", "H4N+", 18.034374128),
            ("[OH-]", "HO-", 17.002739652),
            ("[O-]C(=O)C(=O)[O-]", "C2O4-2", 87.97965848),
        ] {
            let product = finalize(&parse(text));
            assert_eq!(product.composition.formula, formula);
            assert!(
                (product.composition.neutral_atom_mass - mass).abs() < 1e-8,
                "{text}: {:?}",
                product.composition
            );
            let assignment =
                crate::validate::strict_owner_valence(&product.record.topology).unwrap();
            let ionic = cosmolkit_descriptors::exact_molecular_weight_with_valence(
                &product.record.topology,
                false,
                Some(&assignment),
            )
            .unwrap();
            if product.composition.formal_charge == 0 {
                assert_eq!(ionic, product.composition.neutral_atom_mass);
            } else {
                assert_ne!(ionic, product.composition.neutral_atom_mass);
            }
        }
        let ion = finalize(&parse("[23Na+]"));
        assert_eq!(ion.composition.counts.len(), 1);
        assert_eq!(ion.composition.counts[&(11, Some(23))], 1);
        assert_eq!(ion.composition.formal_charge, 1);
        assert!(ion.canonical_smiles.contains("23Na"));
        assert!((ion.composition.neutral_atom_mass - 22.98976928).abs() < 1e-8);
    }

    #[test]
    fn retained_metal_single_bonds_and_radicals_survive_finalization_and_real_identity() {
        for text in [
            "C[Hg+]",
            "C[Hg]C",
            "C[Sn](C)(C)C",
            "C[202Hg]C",
            "C[120Sn](C)(C)C",
            "[13CH3][Hg+]",
            "[CH2]CC[Hg]C",
            "[CH3]",
        ] {
            let record = parse(text);
            let product = finalize(&record);
            assert_eq!(
                product.record.topology.atoms.len(),
                record.topology.atoms.len()
            );
            assert_eq!(
                product.record.topology.bonds.len(),
                record.topology.bonds.len()
            );
            for (before, after) in record
                .topology
                .atoms
                .iter()
                .zip(&product.record.topology.atoms)
            {
                assert_eq!(before.element(), after.element());
                assert_eq!(before.isotope(), after.isotope());
                assert_eq!(before.formal_charge(), after.formal_charge());
                assert_eq!(before.radical_electrons(), after.radical_electrons());
            }
            for bond in &record.topology.bonds {
                let after = &product.record.topology.bonds[bond.id().index()];
                assert_eq!(
                    (after.begin(), after.end(), after.order()),
                    (bond.begin(), bond.end(), bond.order())
                );
            }
            if text.contains("Hg") || text.contains("Sn") {
                assert!(product.record.topology.bonds.iter().any(|b| {
                    [b.begin(), b.end()].iter().any(|a| {
                        matches!(
                            product.record.topology.atoms[a.index()].atomic_number(),
                            50 | 80
                        )
                    }) && b.order() == BondOrder::Single
                }));
            }
        }
    }

    #[test]
    fn finalization_is_idempotent_and_permutation_independent_with_mapped_h() {
        let record = parse("[H]O[13CH2]C");
        let product = finalize(&record);
        assert_same(&product, &finalize(&product.record));
        let (reordered, _) = record
            .topology
            .reordered_atoms(&[
                AtomId::new(3),
                AtomId::new(2),
                AtomId::new(1),
                AtomId::new(0),
            ])
            .unwrap();
        assert_same(
            &product,
            &finalize(&SmilesRecord {
                topology: reordered,
                ..record.clone()
            }),
        );
        let mut mapped = parse("[H]OC");
        mapped.topology.atoms[0].set_atom_map(Some(99));
        let result = finalize(&mapped);
        assert_eq!(result.record.topology.atoms.len(), 2);
        assert_eq!(result.original_to_final.atoms.old_to_new[0], None);
        assert_eq!(result.original_to_final.bonds.old_to_new[0], None);
        result
            .original_to_final
            .validate_for_counts(3, 2, 2, 1)
            .unwrap();
    }
}
