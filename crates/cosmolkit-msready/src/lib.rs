//! Detached, isotope-preserving and stereo-free Native MSReady standardization.
//! This finite published-concept profile does not claim historical workflow parity.
mod charge;
mod components;
mod composition;
mod disconnect;
mod exclusions;
mod finalize;
mod hydrogens;
mod identify;
mod pipeline;
mod profile;
mod rewrite;
mod rules;
mod types;
mod validate;

#[cfg(test)]
use disconnect::disconnect_all_alkali_oxygen;
pub(crate) use exclusions::match_exact_exclusions;
pub use pipeline::{
    EXCLUSION_SHA256, EXCLUSION_VARIANTS_SHA256, NITRATE_TEMPLATES_SHA256, RULE_SHA256, run_info,
    standardize, standardize_smiles,
};
pub use profile::{BEHAVIOR, NITRATE_ADAPTATION_ID, PROFILE_ID, PROFILE_VERSION};
pub use types::*;
#[cfg(test)]
use validate::validate_input;
pub(crate) use validate::{retained_metal_states, validate_syntax};

#[cfg(test)]
mod profile_tests {
    use super::*;
    #[test]
    fn frozen_profile_and_nitrate_comparison_adaptation_are_explicit() {
        assert_eq!(PROFILE_ID, "cosmolkit_msready_v1");
        assert_eq!(PROFILE_VERSION, "1.0.0");
        assert_eq!(BEHAVIOR, "Native");
        assert_eq!(
            NITRATE_ADAPTATION_ID,
            "nitrate_comparison_charge_adaptation_v1"
        );
        assert_eq!(MsReadyParams::default().max_edits, 1024);
        assert_eq!(MsReadyParams::default().max_matches_per_sweep, 4096);
    }
}
