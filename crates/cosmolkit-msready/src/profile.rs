use crate::{MsReadyFailure, MsReadyFailureKind, MsReadyStage};

pub const PROFILE_ID: &str = "cosmolkit_msready_v1";
pub const PROFILE_VERSION: &str = "1.0.0";
pub const BEHAVIOR: &str = "Native";
pub const NITRATE_ADAPTATION_ID: &str = "nitrate_comparison_charge_adaptation_v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MetalPredicate {
    HgCationMonovalent,
    HgNeutralDivalent,
    SnNeutralTetravalent,
    AlkaliOxygenDisconnection,
}

pub(crate) fn unsupported_at(row: usize, context: impl Into<String>) -> MsReadyFailure {
    unsupported(context).with_atoms(vec![row])
}

pub fn unsupported(context: impl Into<String>) -> MsReadyFailure {
    MsReadyFailure {
        evidence: Default::default(),
        kind: MsReadyFailureKind::UnsupportedChemistry,
        stage: MsReadyStage::Validation,
        context: context.into(),
        original: None,
    }
}
