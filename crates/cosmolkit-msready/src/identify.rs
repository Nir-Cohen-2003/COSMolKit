use cosmolkit_core::{
    KekulizeParams, ValenceModel, assign_valence_for_topology,
    assign_valence_with_options_for_topology, element_info, kekulize,
};
use cosmolkit_inchi::{
    InchiAtom, InchiBond, InchiBondType, InchiChiralTag, InchiMolecule, InchiToolkitError,
    MolToInchiToolkit, inchi_to_inchi_key, mol_to_inchi,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element, TopologyBlock,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StandardIdentity {
    pub inchi: String,
    pub full_key: String,
    pub engine_version: &'static str,
    pub return_code: i32,
    pub message: Vec<u8>,
    pub log: Vec<u8>,
    pub diagnostics: Vec<String>,
}

pub(crate) const INCHI_ENGINE_VERSION: &str = "1.07.5";

fn toolkit_error(message: impl Into<String>) -> InchiToolkitError {
    InchiToolkitError {
        kind: "MSReadyToolkitError",
        message: message.into(),
    }
}

fn model_graph(molecule: &InchiMolecule) -> Result<TopologyBlock, InchiToolkitError> {
    let atoms = molecule
        .atoms()
        .iter()
        .map(|a| {
            let element = u8::try_from(a.atomic_number)
                .ok()
                .and_then(Element::from_atomic_number)
                .ok_or_else(|| {
                    toolkit_error(format!("invalid atomic number {}", a.atomic_number))
                })?;
            let mut spec = AtomSpec::new(element)
                .with_formal_charge(
                    i8::try_from(a.formal_charge)
                        .map_err(|_| toolkit_error("formal charge outside model range"))?,
                )
                .with_explicit_hydrogens(
                    u8::try_from(a.num_explicit_hydrogens)
                        .map_err(|_| toolkit_error("explicit H count outside model range"))?,
                )
                .with_aromatic(a.is_aromatic)
                .with_radical_electrons(
                    u8::try_from(a.num_radical_electrons)
                        .map_err(|_| toolkit_error("radical count outside model range"))?,
                )
                .with_no_implicit(a.no_implicit);
            if a.isotope != 0 {
                spec = spec.with_isotope(
                    u16::try_from(a.isotope)
                        .map_err(|_| toolkit_error("isotope outside model range"))?,
                );
            }
            Ok(Atom::from_spec(AtomId::new(0), spec))
        })
        .collect::<Result<Vec<_>, InchiToolkitError>>()?;
    let bonds = molecule
        .bonds()
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let order = match b.bond_type {
                InchiBondType::Single => BondOrder::Single,
                InchiBondType::Double => BondOrder::Double,
                InchiBondType::Triple => BondOrder::Triple,
                InchiBondType::Aromatic => BondOrder::Aromatic,
                other => return Err(toolkit_error(format!("unsupported adapter bond {other:?}"))),
            };
            let spec = BondSpec::new(
                AtomId::new(b.begin_atom_index() as usize),
                AtomId::new(b.end_atom_index() as usize),
                order,
            )
            .with_aromatic(b.is_aromatic);
            Ok(Bond::from_spec(BondId::new(i), spec))
        })
        .collect::<Result<Vec<_>, InchiToolkitError>>()?;
    let atoms = atoms
        .into_iter()
        .enumerate()
        .map(|(i, a)| a.with_id(AtomId::new(i)))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .map_err(|e| toolkit_error(e.to_string()))
}

fn adapter_graph(
    topology: &TopologyBlock,
    molecule: &mut InchiMolecule,
) -> Result<(), InchiToolkitError> {
    let mut atoms = Vec::with_capacity(topology.atoms.len());
    for atom in &topology.atoms {
        let mut adapter = molecule.atoms()[atom.id().index()].clone();
        adapter.is_aromatic = atom.is_aromatic();
        adapter.formal_charge = i32::from(atom.formal_charge());
        adapter.num_explicit_hydrogens = u32::from(atom.explicit_hydrogens());
        adapter.num_radical_electrons = u32::from(atom.radical_electrons());
        atoms.push(adapter);
    }
    let bonds = topology
        .bonds
        .iter()
        .map(|bond| {
            let mut adapter = InchiBond::new(
                bond.begin().index() as u32,
                bond.end().index() as u32,
                match bond.order() {
                    BondOrder::Single => InchiBondType::Single,
                    BondOrder::Double => InchiBondType::Double,
                    BondOrder::Triple => InchiBondType::Triple,
                    BondOrder::Aromatic => InchiBondType::Aromatic,
                    _ => InchiBondType::Unspecified,
                },
            );
            adapter.is_aromatic = bond.is_aromatic();
            adapter
        })
        .collect();
    molecule
        .replace_graph(atoms, bonds, Vec::new())
        .map_err(|e| toolkit_error(format!("adapter graph synchronization: {e:?}")))
}

fn topology_to_inchi_molecule(
    topology: &TopologyBlock,
) -> Result<InchiMolecule, InchiToolkitError> {
    topology
        .validate()
        .map_err(|e| toolkit_error(e.to_string()))?;
    let atoms = topology
        .atoms
        .iter()
        .map(|atom| InchiAtom {
            atomic_number: i32::from(atom.atomic_number()),
            formal_charge: i32::from(atom.formal_charge()),
            num_explicit_hydrogens: u32::from(atom.explicit_hydrogens()),
            is_aromatic: atom.is_aromatic(),
            isotope: u32::from(atom.isotope().unwrap_or(0)),
            num_radical_electrons: u32::from(atom.radical_electrons()),
            no_implicit: atom.no_implicit(),
            chiral_tag: InchiChiralTag::Unspecified,
            cip_rank: None,
            cached_explicit_valence: None,
            cached_implicit_valence: None,
        })
        .collect();
    let bonds = topology
        .bonds
        .iter()
        .map(|bond| {
            let mut adapter = InchiBond::new(
                bond.begin().index() as u32,
                bond.end().index() as u32,
                match bond.order() {
                    BondOrder::Single => InchiBondType::Single,
                    BondOrder::Double => InchiBondType::Double,
                    BondOrder::Triple => InchiBondType::Triple,
                    BondOrder::Aromatic => InchiBondType::Aromatic,
                    unsupported => {
                        return Err(toolkit_error(format!(
                            "unsupported bond for identity: {unsupported:?}"
                        )));
                    }
                },
            );
            // Core kekulization inspects the aromatic flag as well as order.
            // This transient adapter must transport both independent carriers.
            adapter.is_aromatic = bond.is_aromatic();
            Ok(adapter)
        })
        .collect::<Result<Vec<_>, _>>()?;
    InchiMolecule::try_from_graph(atoms, bonds, Vec::new())
        .map_err(|e| toolkit_error(format!("adapter topology conversion: {e:?}")))
}

pub(crate) fn standard_identity(
    topology: &TopologyBlock,
) -> Result<StandardIdentity, crate::MsReadyFailure> {
    use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};
    let fail = |kind, context: String| MsReadyFailure {
        evidence: Default::default(),
        kind,
        stage: Stage::Identifier,
        context,
        original: Some(crate::MsReadyOriginal {
            topology: topology.clone(),
            properties: cosmolkit_model::MoleculeProperties::default(),
        }),
    };
    let comparison = crate::hydrogens::materialize_tracked_isotopic_hydrogens(topology)?;
    let mut adapter = topology_to_inchi_molecule(&comparison)
        .map_err(|e| fail(Kind::UnsupportedIdentifier, e.message))?;
    let generated = mol_to_inchi(&mut OwnerToolkit, &adapter, None)
        .map_err(|e| fail(Kind::UnsupportedIdentifier, e.to_string()))?;
    let inchi = String::from_utf8(generated.inchi.clone())
        .map_err(|e| fail(Kind::InvariantViolation, e.to_string()))?;
    if !inchi.starts_with("InChI=1S/") {
        return Err(fail(
            Kind::InvariantViolation,
            format!("non-Standard InChI prefix: {inchi}"),
        ));
    }
    let keyed = inchi_to_inchi_key(&generated.inchi)
        .map_err(|e| fail(Kind::UnsupportedIdentifier, e.to_string()))?;
    let full_key =
        String::from_utf8(keyed.key).map_err(|e| fail(Kind::InvariantViolation, e.to_string()))?;
    if full_key.len() != 27 || !full_key.is_ascii() {
        return Err(fail(
            Kind::InvariantViolation,
            format!("unexpected full InChIKey {full_key}"),
        ));
    }
    // Adapter state is transient and never replaces the analysis topology.
    let _ = &mut adapter;
    Ok(StandardIdentity {
        inchi,
        full_key,
        engine_version: INCHI_ENGINE_VERSION,
        return_code: generated.return_values.return_code,
        message: generated.return_values.message,
        log: generated.return_values.log,
        diagnostics: generated
            .diagnostics
            .iter()
            .map(|d| d.message.clone())
            .collect(),
    })
}

struct OwnerToolkit;
impl MolToInchiToolkit for OwnerToolkit {
    fn needs_update_property_cache(
        &mut self,
        molecule: &InchiMolecule,
    ) -> Result<bool, InchiToolkitError> {
        Ok(molecule
            .atoms()
            .iter()
            .any(|a| a.cached_explicit_valence.is_none() || a.cached_implicit_valence.is_none()))
    }
    fn update_property_cache(
        &mut self,
        molecule: &mut InchiMolecule,
        strict: bool,
    ) -> Result<(), InchiToolkitError> {
        let topology = model_graph(molecule)?;
        let valence =
            assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, strict)
                .map_err(|e| toolkit_error(e.to_string()))?;
        for (i, atom) in molecule.atom_properties_mut().iter_mut().enumerate() {
            atom.cached_explicit_valence = valence.explicit_valence.get(i).copied();
            atom.cached_implicit_valence = valence.implicit_hydrogens.get(i).copied();
        }
        Ok(())
    }
    fn kekulize(
        &mut self,
        molecule: &mut InchiMolecule,
        mark_atoms_bonds: bool,
    ) -> Result<(), InchiToolkitError> {
        let topology = model_graph(molecule)?;
        let assignment = kekulize(
            &topology,
            &KekulizeParams {
                mark_atoms_bonds,
                ..KekulizeParams::default()
            },
        )
        .map_err(|e| toolkit_error(e.to_string()))?;
        adapter_graph(&assignment.topology, molecule)
    }
    fn element_symbol(&mut self, atomic_number: i32) -> Result<Vec<u8>, InchiToolkitError> {
        let z = u8::try_from(atomic_number)
            .map_err(|_| toolkit_error("atomic number outside model range"))?;
        let element =
            Element::from_atomic_number(z).ok_or_else(|| toolkit_error("unknown atomic number"))?;
        Ok(element.symbol().as_bytes().to_vec())
    }
    fn atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
        let z = u8::try_from(atomic_number)
            .map_err(|_| toolkit_error("atomic number outside model range"))?;
        let element =
            Element::from_atomic_number(z).ok_or_else(|| toolkit_error("unknown atomic number"))?;
        Ok(element_info(element).atomic_weight)
    }
    fn total_num_hydrogens(
        &mut self,
        molecule: &InchiMolecule,
        atom_index: u32,
    ) -> Result<u32, InchiToolkitError> {
        let topology = model_graph(molecule)?;
        let i =
            usize::try_from(atom_index).map_err(|_| toolkit_error("atom index out of range"))?;
        let atom = topology
            .atoms
            .get(i)
            .ok_or_else(|| toolkit_error("atom index out of range"))?;
        let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
            .map_err(|e| toolkit_error(e.to_string()))?;
        let explicit_h_neighbors = topology
            .adjacency
            .neighbors_of(i)
            .iter()
            .filter(|n| topology.atoms[n.atom_index].element() == Element::H)
            .count();
        u32::try_from(
            atom.explicit_hydrogens() as usize
                + explicit_h_neighbors
                + valence.implicit_hydrogens[i].max(0) as usize,
        )
        .map_err(|_| toolkit_error("hydrogen count overflow"))
    }
    fn calc_implicit_valence(
        &mut self,
        molecule: &mut InchiMolecule,
        atom_index: u32,
    ) -> Result<i32, InchiToolkitError> {
        self.update_property_cache(molecule, true)?;
        let i =
            usize::try_from(atom_index).map_err(|_| toolkit_error("atom index out of range"))?;
        molecule
            .atoms()
            .get(i)
            .and_then(|a| a.cached_implicit_valence)
            .ok_or_else(|| toolkit_error("implicit valence missing"))
    }
    fn total_degree(
        &mut self,
        molecule: &InchiMolecule,
        atom_index: u32,
    ) -> Result<u32, InchiToolkitError> {
        let i =
            usize::try_from(atom_index).map_err(|_| toolkit_error("atom index out of range"))?;
        u32::try_from(
            molecule
                .bonds()
                .iter()
                .filter(|b| b.begin_atom_index() as usize == i || b.end_atom_index() as usize == i)
                .count(),
        )
        .map_err(|_| toolkit_error("degree overflow"))
    }
}

#[cfg(test)]
mod identity_tests {
    use super::*;
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles_complete_source};
    fn parsed(smiles: &str) -> TopologyBlock {
        parse_smiles_complete_source(smiles, &SmilesParseParams::default())
            .unwrap()
            .topology
    }
    #[test]
    fn acetone_diagnostic_key_and_unnormalized_enol_are_distinct() {
        let acetone = standard_identity(&parsed("CC(=O)C")).unwrap();
        let enol = standard_identity(&parsed("C=C(O)C")).unwrap();
        assert_eq!(acetone.full_key, "CSCPPACGZOOCGX-UHFFFAOYSA-N");
        assert_ne!(acetone.full_key, enol.full_key);
        assert_ne!(acetone.inchi, enol.inchi);
    }

    #[test]
    fn aromatic_adapter_flags_reach_owner_kekulization_and_real_standard_identity() {
        let aromatic = parsed("c1ccccc1");
        let snapshot = aromatic.clone();
        let adapter = topology_to_inchi_molecule(&aromatic).unwrap();
        assert!(adapter.bonds().iter().all(|bond| bond.is_aromatic));
        let roundtrip = model_graph(&adapter).unwrap();
        for (actual, source) in roundtrip.bonds.iter().zip(&aromatic.bonds) {
            assert_eq!(
                (
                    actual.order(),
                    actual.is_aromatic(),
                    actual.begin(),
                    actual.end()
                ),
                (
                    source.order(),
                    source.is_aromatic(),
                    source.begin(),
                    source.end()
                )
            );
        }
        let identity = standard_identity(&aromatic).unwrap();
        assert_eq!(identity.full_key, "UHOVQNZJYSORNB-UHFFFAOYSA-N");
        let kekule = standard_identity(&parsed("C1=CC=CC=C1")).unwrap();
        assert_eq!(identity.full_key, kekule.full_key);
        assert_eq!(identity.inchi, kekule.inchi);
        assert_eq!(aromatic, snapshot);
    }

    #[test]
    fn tracked_isotopes_match_explicit_identity_and_ranking_without_mutation() {
        let mut keys = Vec::new();
        let mut rankings = Vec::new();
        for isotope in [2, 3] {
            for (parent, explicit) in [
                (2, format!("CCO[{isotope}H]")),
                (0, format!("[{isotope}H]CCO")),
            ] {
                let mut tracked = parsed("CCO");
                tracked.atoms[parent].set_tracked_isotopic_hydrogens(vec![isotope]);
                let snapshot = tracked.clone();
                let identity = standard_identity(&tracked).unwrap();
                let ranking = crate::rewrite::canonical_graph_key(&tracked).unwrap();
                let explicit = parsed(&explicit);
                assert_eq!(
                    identity.full_key,
                    standard_identity(&explicit).unwrap().full_key
                );
                assert_eq!(
                    ranking,
                    crate::rewrite::canonical_graph_key(&explicit).unwrap()
                );
                assert_eq!(tracked, snapshot);
                keys.push(identity.full_key);
                rankings.push(ranking);
            }
        }
        for left in 0..keys.len() {
            for right in left + 1..keys.len() {
                assert_ne!(
                    keys[left], keys[right],
                    "D/T and label positions must remain distinct"
                );
                assert_ne!(rankings[left], rankings[right]);
            }
        }
        // Bracket/no-implicit tracked representation used by proton-moving rules.
        let mut tracked = parsed("CO");
        tracked.atoms[1].set_no_implicit(true);
        tracked.atoms[1].set_tracked_isotopic_hydrogens(vec![2]);
        assert_eq!(
            standard_identity(&tracked).unwrap().full_key,
            standard_identity(&parsed("CO[2H]")).unwrap().full_key
        );
    }

    #[test]
    fn actual_standard_engine_identity_proof_and_diagnostics() {
        for smiles in [
            "CC(=O)C",
            "C=C(O)C",
            "CN",
            "C[NH3+]",
            "C[N+](C)(C)C",
            "[13CH4]",
            "[2H]O[3H]",
            "[CH3]",
            "C[Hg+]",
            "C[Hg]C",
            "C[Sn](C)(C)C",
            "[201Hg]C",
            "C[201Hg+]",
            "[201Sn](C)(C)(C)C",
            "[13CH3][Na+]",
            "[6Li+]",
            "[Cl-]",
        ] {
            let graph = parsed(smiles);
            let result = standard_identity(&graph).unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert!(
                result.inchi.starts_with("InChI=1S/"),
                "{smiles}: {}",
                result.inchi
            );
            assert_eq!(result.full_key.len(), 27, "{smiles}: {}", result.full_key);
            assert_eq!(result.engine_version, "1.07.5");
            assert_eq!(
                graph,
                parsed(smiles),
                "identity conversion changed caller topology: {smiles}"
            );
            eprintln!(
                "{smiles}\t{}\t{}\tengine={} code={} message={:?} log={:?} diagnostics={:?}",
                result.full_key,
                result.inchi,
                result.engine_version,
                result.return_code,
                result.message,
                result.log,
                result.diagnostics
            );
        }
    }
}
