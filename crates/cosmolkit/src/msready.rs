//! Shared-borrow bridge to the detached Native MSReady owner.
use crate::{Molecule, MsReadyFailure, MsReadyParams, MsReadyResult};

impl Molecule {
    /// Derives stereo-free, isotope-preserving MSReady products with default limits.
    ///
    /// Returns detached values, never live molecules, and leaves all input blocks
    /// and runtime caches unchanged. This finite Native profile does not claim
    /// historical KNIME/Indigo/EPA equivalence. Requires opt-in `cap-msready`.
    pub fn with_msready(&self) -> Result<MsReadyResult, MsReadyFailure> {
        self.with_msready_with_params(&MsReadyParams::default())
    }

    /// Derives MSReady products with archival metadata and explicit resource limits.
    ///
    /// Chemistry is frozen by the profile. `MsReadyParams::max_edits` counts
    /// committed edits across all stages/components of the original input.
    /// `MsReadyParams::max_matches_per_sweep` is not an input-wide candidate cap:
    /// one group sweep pools all six appearances and their rematches; a separate
    /// proton sweep pools internal/reservoir phases and their rematches. The
    /// initial disconnection scan is its own sweep. New component/group/proton
    /// sweeps reset candidate counts, never the accumulated edit count.
    ///
    /// A failure retains original input and completed evidence; it does not
    /// return a partial successful collection.
    pub fn with_msready_with_params(
        &self,
        params: &MsReadyParams,
    ) -> Result<MsReadyResult, MsReadyFailure> {
        cosmolkit_msready::standardize(self.topology(), self.properties(), params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn shared_borrow_preserves_all_four_runtime_blocks_on_success_and_failure() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().unwrap();
        let peer = molecule.clone();
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();
        let original_topology = molecule.topology().clone();
        let original_properties = molecule.properties().clone();
        molecule.with_msready().unwrap();
        molecule
            .with_msready_with_params(&MsReadyParams {
                max_edits: 0,
                ..Default::default()
            })
            .unwrap_err();
        for value in [&molecule, &peer] {
            assert!(Arc::ptr_eq(&topology, &value.topology_arc_runtime()));
            assert!(Arc::ptr_eq(&coordinates, &value.coordinates_arc_runtime()));
            assert!(Arc::ptr_eq(&properties, &value.properties_arc_runtime()));
            assert!(Arc::ptr_eq(&cache, &value.derived_cache_arc_runtime()));
            assert_eq!(value.topology(), &original_topology);
            assert_eq!(value.properties(), &original_properties);
        }
    }
}
