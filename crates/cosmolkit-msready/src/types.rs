use cosmolkit_model::{MoleculeProperties, TopologyBlock, TopologyMapping};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MsReadyStage {
    Input,
    Validation,
    Disconnection,
    Exclusion,
    Normalization,
    Finalization,
    Identifier,
    Accounting,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MsReadyFailureKind {
    InvalidInput,
    InvalidValence,
    UnsupportedChemistry,
    UnsupportedIdentifier,
    InvariantViolation,
    OwnerError,
    TransformationLimit,
    TransformationCycle,
}

/// All atom/bond rows refer to the immutable original input, not fragment rows.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MsReadyLocation {
    pub input_component: Option<usize>,
    pub component_index: Option<usize>,
    pub rule_id: Option<String>,
    pub atom_indices: Vec<usize>,
    pub bond_indices: Vec<usize>,
    pub owner_category: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MsReadyInputLink {
    pub record_id: Option<String>,
    /// Optional caller archival metadata, not the actual parsed text.
    pub original_smiles: Option<String>,
    pub source_smiles: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MsReadyFailure {
    pub kind: MsReadyFailureKind,
    pub stage: MsReadyStage,
    pub context: String,
    pub original: Option<MsReadyOriginal>,
    pub evidence: MsReadyFailureEvidence,
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MsReadyFailureEvidence {
    pub input_link: MsReadyInputLink,
    pub location: MsReadyLocation,
    pub exclusions: Vec<MsReadyExclusion>,
    pub events: Vec<MsReadyEvent>,
    pub warnings: Vec<MsReadyWarning>,
}
impl MsReadyFailure {
    pub(crate) fn new(
        kind: MsReadyFailureKind,
        stage: MsReadyStage,
        context: impl Into<String>,
    ) -> Self {
        Self {
            kind,
            stage,
            context: context.into(),
            original: None,
            evidence: MsReadyFailureEvidence::default(),
        }
    }
    pub(crate) fn with_atoms(mut self, rows: Vec<usize>) -> Self {
        self.evidence.location.atom_indices = rows;
        self
    }
    pub(crate) fn with_bonds(mut self, rows: Vec<usize>) -> Self {
        self.evidence.location.bond_indices = rows;
        self
    }
}
impl std::fmt::Display for MsReadyFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?} at {:?}: {}", self.kind, self.stage, self.context)
    }
}
impl std::error::Error for MsReadyFailure {}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyParams {
    pub record_id: Option<String>,
    pub original_smiles: Option<String>,
    pub max_edits: usize,
    pub max_matches_per_sweep: usize,
}
impl Default for MsReadyParams {
    fn default() -> Self {
        Self {
            record_id: None,
            original_smiles: None,
            max_edits: 1024,
            max_matches_per_sweep: 4096,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct MsReadyOriginal {
    pub topology: TopologyBlock,
    pub properties: MoleculeProperties,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyRunInfo {
    pub profile_id: &'static str,
    pub profile_version: &'static str,
    pub behavior: &'static str,
    pub rule_sha256: &'static str,
    pub exclusion_sha256: &'static str,
    pub exclusion_variants_sha256: &'static str,
    pub nitrate_templates_sha256: &'static str,
    pub inchi_engine_version: &'static str,
    pub writer_params: cosmolkit_smiles::SmilesWriteParams,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MsReadyResult {
    pub original: MsReadyOriginal,
    pub input_link: MsReadyInputLink,
    pub run_info: MsReadyRunInfo,
    pub status: MsReadyStatus,
    pub exclusions: Vec<MsReadyExclusion>,
    pub events: Vec<MsReadyEvent>,
    pub warnings: Vec<MsReadyWarning>,
    pub atom_dispositions: Vec<MsReadyAtomDisposition>,
}
#[derive(Clone, Debug, PartialEq)]
pub enum MsReadyStatus {
    Products(Vec<MsReadyProduct>),
    NoEligibleComponent,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MsReadyProduct {
    pub topology: TopologyBlock,
    pub canonical_smiles: String,
    pub inchi: String,
    pub inchi_key: String,
    pub identifier: MsReadyIdentifierDiagnostics,
    pub formula: String,
    pub isotope_composition: Vec<MsReadyCompositionEntry>,
    pub formal_charge: i32,
    pub neutral_atom_mass: f64,
    pub origins: Vec<MsReadyOrigin>,
    pub multiplicity: usize,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyIdentifierDiagnostics {
    pub engine_version: &'static str,
    pub return_code: i32,
    pub message: Vec<u8>,
    pub log: Vec<u8>,
    pub diagnostics: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyCompositionEntry {
    pub atomic_number: u8,
    pub isotope: Option<u16>,
    pub count: u32,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyOrigin {
    pub input_component: usize,
    pub component_index: usize,
    pub final_component_index: usize,
    pub atom_indices: Vec<usize>,
    pub bond_indices: Vec<usize>,
    pub topology_mapping: TopologyMapping,
    pub hydrogen_dispositions: Vec<MsReadyHydrogenDisposition>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyExclusion {
    pub source_id: String,
    pub source_sha256: &'static str,
    pub source_field: String,
    pub template_id: Option<String>,
    pub adaptation_id: Option<String>,
    pub comparison_only: bool,
    pub matched_formula: String,
    pub matched_formal_charge: i32,
    pub origin: MsReadyOrigin,
    pub atom_mappings: Vec<Vec<usize>>,
    pub folded_hydrogen_mappings: Vec<(usize, usize)>,
    pub reason: &'static str,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MsReadyEventKind {
    ParserRadicalInitialization,
    AlkaliDisconnection,
    GroupRewrite,
    InternalProtonTransfer,
    ReservoirDeprotonation,
    ReservoirProtonation,
    LabelledReservoirBlocked,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MsReadyProtonCarrier {
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
pub struct MsReadyEvent {
    pub kind: MsReadyEventKind,
    pub stage: MsReadyStage,
    pub rule_id: String,
    pub appearance_order: Option<usize>,
    pub location: MsReadyLocation,
    pub donor_input_row: Option<usize>,
    pub acceptor_input_row: Option<usize>,
    pub carrier: MsReadyProtonCarrier,
    pub delta_hydrogens: i32,
    pub delta_charge: i32,
    pub atom_charge_deltas: Vec<(usize, i8)>,
    pub hydrogen_counts_before: Vec<(usize, u32)>,
    pub radical_electrons_before: Option<u8>,
    pub radical_electrons_after: Option<u8>,
    pub formula_before: String,
    pub formula_after: String,
    pub total_charge_before: i32,
    pub total_charge_after: i32,
    pub topology_mapping: TopologyMapping,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyWarning {
    pub code: String,
    pub message: String,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MsReadyHydrogenDispositionKind {
    FoldedInto {
        input_parent_row: usize,
    },
    Retained {
        output_atom_row: usize,
    },
    RetainedIsotope {
        output_atom_row: usize,
        isotope: u16,
    },
    DeletedToReservoir,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MsReadyHydrogenDisposition {
    pub input_atom_row: usize,
    pub kind: MsReadyHydrogenDispositionKind,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MsReadyAtomDispositionKind {
    Retained {
        product_index: usize,
        origin_index: usize,
        output_atom_row: usize,
    },
    FoldedHydrogen {
        product_index: usize,
        origin_index: usize,
        input_parent_row: usize,
    },
    DeletedToReservoir,
    Excluded {
        exclusion_indices: Vec<usize>,
    },
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsReadyAtomDisposition {
    pub input_atom_row: usize,
    pub kind: MsReadyAtomDispositionKind,
}
