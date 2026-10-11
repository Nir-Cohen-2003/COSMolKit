use cosmolkit_core::{ValenceAssignment, ValenceModel, assign_valence_with_options_for_topology};
use cosmolkit_model::{BondOrder, TopologyBlock};
use cosmolkit_types::Element;

use crate::{MsReadyFailure, MsReadyFailureKind as Kind, MsReadyStage as Stage};

pub(crate) const RULE_DATA: &str = include_str!("data/rules.json");

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum RuleId {
    NitroInitial,
    EnamineToImine,
    EnolToCarbonyl,
    YnolToKetene,
    AzideMesomer,
    NitroFinal,
}

impl RuleId {
    pub(crate) const ORDER: [Self; 6] = [
        Self::NitroInitial,
        Self::EnamineToImine,
        Self::EnolToCarbonyl,
        Self::YnolToKetene,
        Self::AzideMesomer,
        Self::NitroFinal,
    ];

    pub(crate) const fn id(self) -> &'static str {
        match self {
            Self::NitroInitial => "nitro_initial",
            Self::EnamineToImine => "enamine_to_imine",
            Self::EnolToCarbonyl => "enol_to_carbonyl",
            Self::YnolToKetene => "ynol_to_ketene",
            Self::AzideMesomer => "azide_mesomer",
            Self::NitroFinal => "nitro_final",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum HydrogenCarrier {
    Implicit,
    Bracket,
    ExplicitAtom {
        atom: usize,
        bond: usize,
        isotope: Option<u16>,
    },
    TrackedIsotope {
        index: usize,
        isotope: u16,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct HydrogenTransfer {
    pub donor: usize,
    pub acceptor: usize,
    pub carrier: HydrogenCarrier,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BondChange {
    pub row: usize,
    pub order: BondOrder,
    pub endpoints: Option<(usize, usize)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ChargeChange {
    pub row: usize,
    pub charge: i8,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RuleCandidate {
    pub rule: RuleId,
    /// Core atoms in stable source-role order, followed by matched R atoms.
    pub atom_rows: Vec<usize>,
    /// Source bonds corresponding to the template's core bonds.
    pub bond_rows: Vec<usize>,
    pub bond_changes: Vec<BondChange>,
    pub charge_changes: Vec<ChargeChange>,
    pub hydrogen_transfer: Option<HydrogenTransfer>,
}

pub(crate) fn validate_rule_data() -> Result<(), MsReadyFailure> {
    let document: serde_json::Value = serde_json::from_str(RULE_DATA).map_err(|error| {
        failure(
            Kind::InvariantViolation,
            format!("cannot decode packaged rule data: {error}"),
        )
    })?;
    if document
        .get("schema_version")
        .and_then(|value| value.as_str())
        != Some("1.0.0")
        || document.get("profile_id").and_then(|value| value.as_str()) != Some(crate::PROFILE_ID)
    {
        return Err(failure(
            Kind::InvariantViolation,
            "rule data profile or schema version does not match the Native profile",
        ));
    }
    let appearances = document
        .get("appearances")
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| failure(Kind::InvariantViolation, "rule appearances are missing"))?;
    if appearances.len() != RuleId::ORDER.len() {
        return Err(failure(
            Kind::InvariantViolation,
            format!(
                "rule data has {} appearances, expected 6",
                appearances.len()
            ),
        ));
    }
    for (index, (entry, rule)) in appearances.iter().zip(RuleId::ORDER).enumerate() {
        if entry.get("rule_id").and_then(|value| value.as_str()) != Some(rule.id())
            || entry.get("order").and_then(serde_json::Value::as_u64) != Some(index as u64 + 1)
        {
            return Err(failure(
                Kind::InvariantViolation,
                format!(
                    "packaged rule appearance {} is out of frozen order",
                    rule.id()
                ),
            ));
        }
        let source_hash = entry
            .get("source_sha256")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        if source_hash.len() != 64 || !source_hash.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            return Err(failure(
                Kind::InvariantViolation,
                format!("rule {} has an invalid source hash", rule.id()),
            ));
        }
    }
    Ok(())
}

pub(crate) fn find_candidates(
    topology: &TopologyBlock,
    rule: RuleId,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    find_candidates_limited(topology, rule, usize::MAX)
}

pub(crate) fn find_candidates_limited(
    topology: &TopologyBlock,
    rule: RuleId,
    match_limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    topology
        .validate()
        .map_err(|error| failure(Kind::InvalidInput, error.to_string()))?;
    match rule {
        RuleId::NitroInitial => find_nitro_initial(topology, match_limit),
        RuleId::EnamineToImine => find_enamine_to_imine(topology, match_limit),
        RuleId::EnolToCarbonyl => find_enol_to_carbonyl(topology, match_limit),
        RuleId::YnolToKetene => find_ynol_to_ketene(topology, match_limit),
        RuleId::AzideMesomer => find_azide_mesomer(topology, match_limit),
        RuleId::NitroFinal => find_nitro_final(topology, match_limit),
    }
}

fn find_nitro_initial(
    topology: &TopologyBlock,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: workflow_configuration/Tauto_Neut (#3652)/N mesos (#3406)/Table Creator (#2133)/settings.xml; SHA-256 120847c464cf0ffd4899c088ec8a3a808d53b51e374cdbc32f8a70c9c5038a6e
    // MSReady❗❗: source-backed Native interpretation; the archive is static and has no executable Indigo oracle. Canonical tie selection serializes every concrete candidate; relative source matcher complexity cannot be established from the archive.
    // <entry key="2" type="xstring" value="[*]N(=O)=O&gt;&gt;[O-][N+]([*])=O"/>
    // END MSREADY MRV SOURCE
    find_nitro(topology, RuleId::NitroInitial, limit)
}

fn find_nitro_final(
    topology: &TopologyBlock,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: workflow_configuration/Tauto_Neut (#3652)/N mesos (#3430)/Table Creator (#2158)/settings.xml; SHA-256 52e77e915133920e6d7ac35968b5a69d493df16e50a95499801068fb924dc528
    // MSReady❗❗: this is a distinct final-stage appearance with the same raw reaction; source execution is unavailable, and the per-candidate graph clone/canonical write cost has no measured source baseline.
    // <entry key="2" type="xstring" value="[*]N(=O)=O&gt;&gt;[O-][N+]([*])=O"/>
    // END MSREADY MRV SOURCE
    find_nitro(topology, RuleId::NitroFinal, limit)
}

fn find_nitro(
    topology: &TopologyBlock,
    rule: RuleId,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    let mut candidates = Vec::new();
    for (n_row, nitrogen) in topology.atoms.iter().enumerate() {
        if nitrogen.element() != Element::N
            || !matches!(nitrogen.formal_charge(), 0 | 1)
            || nitrogen.is_aromatic()
            || nitrogen.radical_electrons() != 0
            || nitrogen.explicit_hydrogens() != 0
            || !nitrogen.tracked_isotopic_hydrogens().is_empty()
        {
            continue;
        }
        let mut oxygens = Vec::with_capacity(2);
        let mut substituent = None;
        let mut core_bonds = Vec::with_capacity(3);
        let mut valid = true;
        for neighbor in topology.adjacency.neighbors_of(n_row) {
            let bond = &topology.bonds[neighbor.bond.index()];
            let atom = &topology.atoms[neighbor.atom_index];
            if bond.query().is_some() || bond.is_aromatic() {
                valid = false;
                break;
            }
            match (atom.element(), bond.order()) {
                (Element::O, BondOrder::Single | BondOrder::Double)
                    if ((bond.order() == BondOrder::Double && atom.formal_charge() == 0)
                        || (nitrogen.formal_charge() == 1
                            && bond.order() == BondOrder::Single
                            && atom.formal_charge() == -1))
                        && atom.explicit_hydrogens() == 0
                        && atom.tracked_isotopic_hydrogens().is_empty()
                        && atom.radical_electrons() == 0
                        && !atom.is_aromatic()
                        && topology
                            .adjacency
                            .neighbors_of(neighbor.atom_index)
                            .iter()
                            .all(|other| other.atom_index == n_row) =>
                {
                    oxygens.push(neighbor.atom_index);
                    core_bonds.push(neighbor.bond.index());
                }
                (element, BondOrder::Single) if element != Element::H && substituent.is_none() => {
                    substituent = Some(neighbor.atom_index);
                    core_bonds.push(neighbor.bond.index());
                }
                _ => {
                    valid = false;
                    break;
                }
            }
        }
        if !valid || oxygens.len() != 2 || substituent.is_none() || core_bonds.len() != 3 {
            continue;
        }
        let double_count = core_bonds
            .iter()
            .filter(|&&row| topology.bonds[row].order() == BondOrder::Double)
            .count();
        if double_count != if nitrogen.formal_charge() == 0 { 2 } else { 1 } {
            continue;
        }
        // Both assignments are considered; isotopic asymmetry is settled by
        // the canonical product graph in rewrite.rs, never by atom row.
        for (single_oxygen, double_oxygen) in [(oxygens[0], oxygens[1]), (oxygens[1], oxygens[0])] {
            let single_bond = topology
                .adjacency
                .neighbors_of(n_row)
                .iter()
                .find(|neighbor| neighbor.atom_index == single_oxygen)
                .expect("nitro oxygen was found in N adjacency")
                .bond
                .index();
            let double_bond = topology
                .adjacency
                .neighbors_of(n_row)
                .iter()
                .find(|neighbor| neighbor.atom_index == double_oxygen)
                .expect("matched O")
                .bond
                .index();
            let mut bonds = core_bonds.clone();
            bonds.sort_unstable();
            if candidates.len() >= limit {
                return Err(failure(
                    Kind::TransformationLimit,
                    format!("rule candidate count exceeds limit {limit}"),
                ));
            }
            candidates.push(RuleCandidate {
                rule,
                atom_rows: vec![
                    n_row,
                    substituent.expect("checked R"),
                    single_oxygen,
                    double_oxygen,
                ],
                bond_rows: bonds,
                bond_changes: vec![
                    BondChange {
                        row: single_bond,
                        order: BondOrder::Single,
                        endpoints: None,
                    },
                    BondChange {
                        row: double_bond,
                        order: BondOrder::Double,
                        endpoints: None,
                    },
                ],
                charge_changes: vec![
                    ChargeChange {
                        row: n_row,
                        charge: 1,
                    },
                    ChargeChange {
                        row: single_oxygen,
                        charge: -1,
                    },
                    ChargeChange {
                        row: double_oxygen,
                        charge: 0,
                    },
                ],
                hydrogen_transfer: None,
            });
        }
    }
    Ok(candidates)
}

fn find_enamine_to_imine(
    topology: &TopologyBlock,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: mrv/Tauto_Neut (#3652)__Metanode (#3383)__MarvinSketch (#3381).mrv; SHA-256 6b9e65c06cfd0dcf5bd391f64f63c9dfc0205a0087fec32bb90b551825a9724f
    // MSReady❗❗: verbatim R-group/bond correspondence is anchored below; Native concretization is not an executed workflow proof. The matcher scans adjacency and the rewrite engine clones/writes per candidate; relative source complexity is unresolved.
    // <atomArray><atom id="a1" elementType="C" x2="-7.700000000000001" y2="1.5950666666666669"></atom><atom id="a2" elementType="N" x2="-6.104933333333334" y2="1.7050133333333335"></atom><atom id="a3" elementType="R" rgroupRef="1" x2="-8.559973333333335" y2="2.8726133333333337" ligandOrder="a1"></atom><atom id="a4" elementType="C" x2="-8.376480000000000" y2="0.21149333333333337"></atom><atom id="a5" elementType="R" rgroupRef="2" x2="-9.912746666666667" y2="0.10565333333333335" ligandOrder="a4"></atom><atom id="a6" elementType="R" rgroupRef="3" x2="-5.244960000000001" y2="0.4274666666666667" ligandOrder="a2"></atom><atom id="a7" elementType="R" rgroupRef="4" x2="-7.5166933333333334" y2="-1.0660533333333335" ligandOrder="a4"></atom><atom id="a8" elementType="H" x2="-5.428505685081719" y2="3.0885037151867056"></atom></atomArray>
    // <bondArray><bond id="b1" atomRefs2="a1 a2" order="1"></bond><bond id="b2" atomRefs2="a1 a3" order="1"></bond><bond id="b3" atomRefs2="a1 a4" order="2"></bond><bond id="b4" atomRefs2="a4 a5" order="1"></bond><bond id="b5" atomRefs2="a2 a6" order="1"></bond><bond id="b6" atomRefs2="a4 a7" order="1"></bond><bond id="b7" atomRefs2="a2 a8" order="1"></bond></bondArray>
    // <bondArray><bond id="b1" atomRefs2="a1 a2" order="2"></bond><bond id="b2" atomRefs2="a2 a3" order="1"></bond><bond id="b3" atomRefs2="a2 a4" order="1"></bond><bond id="b4" atomRefs2="a4 a5" order="1"></bond><bond id="b5" atomRefs2="a1 a6" order="1"></bond><bond id="b6" atomRefs2="a4 a7" order="1"></bond></bondArray>
    // END MSREADY MRV SOURCE
    let valence = assign_for_match(topology)?;
    let mut candidates = Vec::new();
    for (n_row, nitrogen) in topology.atoms.iter().enumerate() {
        if nitrogen.element() != Element::N
            || nitrogen.formal_charge() != 0
            || nitrogen.is_aromatic()
            || nitrogen.radical_electrons() != 0
        {
            continue;
        }
        for n_neighbor in topology.adjacency.neighbors_of(n_row) {
            let alpha = n_neighbor.atom_index;
            if topology.atoms[alpha].element() != Element::C
                || topology.atoms[alpha].formal_charge() != 0
                || topology.atoms[alpha].is_aromatic()
                || topology.atoms[alpha].radical_electrons() != 0
                || topology.bonds[n_neighbor.bond.index()].order() != BondOrder::Single
                || topology.bonds[n_neighbor.bond.index()].query().is_some()
            {
                continue;
            }
            let n_heavy: Vec<_> = topology
                .adjacency
                .neighbors_of(n_row)
                .iter()
                .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() != 1)
                .collect();
            if n_heavy.iter().any(|neighbor| neighbor.atom_index != alpha)
                && n_heavy
                    .iter()
                    .filter(|neighbor| neighbor.atom_index != alpha)
                    .count()
                    != 1
            {
                continue;
            }
            if !n_heavy.iter().any(|neighbor| neighbor.atom_index == alpha) {
                continue;
            }
            let n_h = hydrogen_count(topology, &valence, n_row)?;
            let n_r_branches = n_heavy
                .iter()
                .filter(|neighbor| neighbor.atom_index != alpha)
                .count();
            // The source N has both the donated H and an R slot. A second
            // owner-assigned H may fill that R slot, but the donated H cannot.
            if n_h == 0 || n_r_branches + (n_h.saturating_sub(1) as usize) < 1 {
                continue;
            }
            let mut n_bonds = Vec::new();
            let mut attachment_ok = true;
            for neighbor in topology.adjacency.neighbors_of(n_row) {
                let bond = &topology.bonds[neighbor.bond.index()];
                if bond.order() != BondOrder::Single || bond.query().is_some() || bond.is_aromatic()
                {
                    attachment_ok = false;
                    break;
                }
                n_bonds.push(neighbor.bond.index());
            }
            if !attachment_ok {
                continue;
            }
            for c_neighbor in topology.adjacency.neighbors_of(alpha) {
                if c_neighbor.atom_index == n_row
                    || topology.atoms[c_neighbor.atom_index].element() != Element::C
                    || topology.atoms[c_neighbor.atom_index].formal_charge() != 0
                    || topology.atoms[c_neighbor.atom_index].is_aromatic()
                    || topology.atoms[c_neighbor.atom_index].radical_electrons() != 0
                    || topology.bonds[c_neighbor.bond.index()].order() != BondOrder::Double
                    || topology.bonds[c_neighbor.bond.index()].query().is_some()
                    || topology.bonds[c_neighbor.bond.index()].is_aromatic()
                {
                    continue;
                }
                let beta = c_neighbor.atom_index;
                if !single_attachments(topology, alpha, &[n_row, beta])
                    || !single_attachments(topology, beta, &[alpha])
                    || attachment_slot_count(topology, &valence, alpha, &[n_row, beta])? < 1
                    || attachment_slot_count(topology, &valence, beta, &[alpha])? < 2
                {
                    continue;
                }
                let h_bond = topology
                    .adjacency
                    .neighbors_of(n_row)
                    .iter()
                    .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 1)
                    .map(|neighbor| neighbor.bond.index())
                    .collect::<Vec<_>>();
                let mut core_bonds = vec![n_neighbor.bond.index(), c_neighbor.bond.index()];
                core_bonds.extend(n_bonds.iter().copied());
                core_bonds.extend(h_bond);
                core_bonds.sort_unstable();
                core_bonds.dedup();
                for transfer in hydrogen_transfers(topology, &valence, n_row, beta)? {
                    if candidates.len() >= limit {
                        return Err(failure(
                            Kind::TransformationLimit,
                            format!("rule candidate count exceeds limit {limit}"),
                        ));
                    }
                    candidates.push(RuleCandidate {
                        rule: RuleId::EnamineToImine,
                        atom_rows: vec![n_row, alpha, beta],
                        bond_rows: core_bonds.clone(),
                        bond_changes: vec![
                            BondChange {
                                row: n_neighbor.bond.index(),
                                order: BondOrder::Double,
                                endpoints: None,
                            },
                            BondChange {
                                row: c_neighbor.bond.index(),
                                order: BondOrder::Single,
                                endpoints: None,
                            },
                        ],
                        charge_changes: Vec::new(),
                        hydrogen_transfer: Some(transfer),
                    });
                }
            }
        }
    }
    Ok(candidates)
}

fn find_enol_to_carbonyl(
    topology: &TopologyBlock,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2253).mrv; SHA-256 26c8b16905fc41235cb9a63b4fffcab237b5952774f0d17144107d6b84ab28a9
    // MSReady❗❗: mapped C/O atoms and each labeled R group remain attached to their original atoms. The archive is static; the concrete-H-slot adaptation and per-candidate canonicalization have no executed upstream/complexity baseline.
    // <atomArray><atom id="a1" elementType="R" rgroupRef="1" x2="-0.5513627814580033" y2="2.843655795071509" ligandOrder="a2"></atom><atom id="a2" elementType="C" x2="-0.5513676717940985" y2="1.303655795079274"></atom><atom id="a3" elementType="R" rgroupRef="3" x2="0.7823090048591634" y2="-1.006382321190216" ligandOrder="a4"></atom><atom id="a4" elementType="C" x2="0.7823090048591652" y2="0.5336515599278648"></atom><atom id="a5" elementType="O" x2="2.1159930169963985" y2="1.3036430896017572"></atom><atom id="a6" elementType="R" rgroupRef="2" x2="-1.8850492387834574" y2="0.533660030238448" ligandOrder="a2"></atom><atom id="a7" elementType="H" x2="3.4496672484614432" y2="0.533634619306707"></atom></atomArray>
    // <bondArray><bond id="b1" atomRefs2="a2 a1" order="1"></bond><bond id="b2" atomRefs2="a4 a2" order="2"></bond><bond id="b3" atomRefs2="a3 a4" order="1"></bond><bond id="b4" atomRefs2="a4 a5" order="1"></bond><bond id="b5" atomRefs2="a2 a6" order="1"></bond><bond id="b6" atomRefs2="a5 a7" order="1"></bond></bondArray>
    // <bondArray><bond id="b1" atomRefs2="a2 a1" order="2"></bond><bond id="b2" atomRefs2="a3 a2" order="1"></bond><bond id="b3" atomRefs2="a2 a4" order="1"></bond><bond id="b4" atomRefs2="a4 a5" order="1"></bond><bond id="b5" atomRefs2="a4 a6" order="1"></bond></bondArray>
    // END MSREADY MRV SOURCE
    let valence = assign_for_match(topology)?;
    let mut candidates = Vec::new();
    for (o_row, oxygen) in topology.atoms.iter().enumerate() {
        if oxygen.element() != Element::O
            || oxygen.formal_charge() != 0
            || oxygen.is_aromatic()
            || oxygen.radical_electrons() != 0
            || hydrogen_count(topology, &valence, o_row)? != 1
        {
            continue;
        }
        for o_neighbor in topology.adjacency.neighbors_of(o_row) {
            let alpha = o_neighbor.atom_index;
            if topology.atoms[alpha].element() != Element::C
                || topology.atoms[alpha].formal_charge() != 0
                || topology.atoms[alpha].is_aromatic()
                || topology.atoms[alpha].radical_electrons() != 0
                || topology.bonds[o_neighbor.bond.index()].order() != BondOrder::Single
                || topology.bonds[o_neighbor.bond.index()].query().is_some()
                || topology.bonds[o_neighbor.bond.index()].is_aromatic()
            {
                continue;
            }
            if topology
                .adjacency
                .neighbors_of(o_row)
                .iter()
                .any(|neighbor| {
                    neighbor.atom_index != alpha
                        && topology.atoms[neighbor.atom_index].atomic_number() != 1
                })
            {
                continue;
            }
            for c_neighbor in topology.adjacency.neighbors_of(alpha) {
                if c_neighbor.atom_index == o_row
                    || topology.atoms[c_neighbor.atom_index].element() != Element::C
                    || topology.atoms[c_neighbor.atom_index].formal_charge() != 0
                    || topology.atoms[c_neighbor.atom_index].is_aromatic()
                    || topology.atoms[c_neighbor.atom_index].radical_electrons() != 0
                    || topology.bonds[c_neighbor.bond.index()].order() != BondOrder::Double
                    || topology.bonds[c_neighbor.bond.index()].query().is_some()
                    || topology.bonds[c_neighbor.bond.index()].is_aromatic()
                {
                    continue;
                }
                let beta = c_neighbor.atom_index;
                if !single_attachments(topology, alpha, &[o_row, beta])
                    || !single_attachments(topology, beta, &[alpha])
                    || attachment_slot_count(topology, &valence, alpha, &[o_row, beta])? < 1
                    || attachment_slot_count(topology, &valence, beta, &[alpha])? < 2
                {
                    continue;
                }
                let mut core_bonds = vec![o_neighbor.bond.index(), c_neighbor.bond.index()];
                core_bonds.extend(
                    topology
                        .adjacency
                        .neighbors_of(o_row)
                        .iter()
                        .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 1)
                        .map(|neighbor| neighbor.bond.index()),
                );
                core_bonds.sort_unstable();
                core_bonds.dedup();
                for transfer in hydrogen_transfers(topology, &valence, o_row, beta)? {
                    if candidates.len() >= limit {
                        return Err(failure(
                            Kind::TransformationLimit,
                            format!("rule candidate count exceeds limit {limit}"),
                        ));
                    }
                    candidates.push(RuleCandidate {
                        rule: RuleId::EnolToCarbonyl,
                        atom_rows: vec![beta, alpha, o_row],
                        bond_rows: core_bonds.clone(),
                        bond_changes: vec![
                            BondChange {
                                row: c_neighbor.bond.index(),
                                order: BondOrder::Single,
                                endpoints: None,
                            },
                            BondChange {
                                row: o_neighbor.bond.index(),
                                order: BondOrder::Double,
                                endpoints: None,
                            },
                        ],
                        charge_changes: Vec::new(),
                        hydrogen_transfer: Some(transfer),
                    });
                }
            }
        }
    }
    Ok(candidates)
}

fn find_ynol_to_ketene(
    topology: &TopologyBlock,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2136).mrv; SHA-256 48379a213859aa8d27ef975d23bb1ba9ed0d2441bede8e7f8dd69858ef022bcc
    // MSReady❗❗: the R slot is concretized to the distal C atom and its branches remain in place. The source archive is not executable, and canonical serialization is performed per candidate without a comparable source benchmark.
    // <atomArray><atom id="a1" elementType="C" x2="-7.700000000000001" y2="1.5950666666666669"></atom><atom id="a2" elementType="O" x2="-6.104933333333334" y2="1.7050133333333335"></atom><atom id="a3" elementType="R" rgroupRef="1" x2="-8.559973333333335" y2="2.8726133333333337" ligandOrder="a1"></atom><atom id="a4" elementType="H" x2="-5.245044247222097" y2="0.427441120271864"></atom></atomArray>
    // <bondArray><bond id="b1" atomRefs2="a1 a2" order="1"></bond><bond id="b2" atomRefs2="a1 a3" order="3"></bond><bond id="b3" atomRefs2="a2 a4" order="1"></bond></bondArray>
    // <bondArray><bond id="b1" atomRefs2="a1 a2" order="2"></bond><bond id="b2" atomRefs2="a2 a3" order="2"></bond></bondArray>
    // END MSREADY MRV SOURCE
    let valence = assign_for_match(topology)?;
    let mut candidates = Vec::new();
    for (o_row, oxygen) in topology.atoms.iter().enumerate() {
        if oxygen.element() != Element::O
            || oxygen.formal_charge() != 0
            || oxygen.is_aromatic()
            || oxygen.radical_electrons() != 0
            || hydrogen_count(topology, &valence, o_row)? != 1
        {
            continue;
        }
        for o_neighbor in topology.adjacency.neighbors_of(o_row) {
            let central = o_neighbor.atom_index;
            if topology.atoms[central].element() != Element::C
                || topology.atoms[central].formal_charge() != 0
                || topology.atoms[central].is_aromatic()
                || topology.atoms[central].radical_electrons() != 0
                || topology.bonds[o_neighbor.bond.index()].order() != BondOrder::Single
                || topology.bonds[o_neighbor.bond.index()].query().is_some()
                || topology.bonds[o_neighbor.bond.index()].is_aromatic()
            {
                continue;
            }
            if topology
                .adjacency
                .neighbors_of(o_row)
                .iter()
                .any(|neighbor| {
                    neighbor.atom_index != central
                        && topology.atoms[neighbor.atom_index].atomic_number() != 1
                })
            {
                continue;
            }
            if topology
                .adjacency
                .neighbors_of(central)
                .iter()
                .any(|neighbor| {
                    neighbor.atom_index != o_row
                        && (topology.atoms[neighbor.atom_index].element() != Element::C
                            || topology.bonds[neighbor.bond.index()].order() != BondOrder::Triple
                            || topology.bonds[neighbor.bond.index()].query().is_some()
                            || topology.bonds[neighbor.bond.index()].is_aromatic())
                })
            {
                continue;
            }
            let distal_neighbors: Vec<_> = topology
                .adjacency
                .neighbors_of(central)
                .iter()
                .filter(|neighbor| neighbor.atom_index != o_row)
                .collect();
            if distal_neighbors.len() != 1 {
                continue;
            }
            let distal = distal_neighbors[0].atom_index;
            let triple_bond = distal_neighbors[0].bond.index();
            if topology.atoms[distal].formal_charge() != 0
                || topology.atoms[distal].is_aromatic()
                || topology.atoms[distal].radical_electrons() != 0
                || !single_attachments(topology, distal, &[central])
            {
                continue;
            }
            let mut core_bonds = vec![o_neighbor.bond.index(), triple_bond];
            core_bonds.extend(
                topology
                    .adjacency
                    .neighbors_of(o_row)
                    .iter()
                    .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 1)
                    .map(|neighbor| neighbor.bond.index()),
            );
            core_bonds.sort_unstable();
            core_bonds.dedup();
            for transfer in hydrogen_transfers(topology, &valence, o_row, distal)? {
                if candidates.len() >= limit {
                    return Err(failure(
                        Kind::TransformationLimit,
                        format!("rule candidate count exceeds limit {limit}"),
                    ));
                }
                candidates.push(RuleCandidate {
                    rule: RuleId::YnolToKetene,
                    atom_rows: vec![o_row, central, distal],
                    bond_rows: core_bonds.clone(),
                    bond_changes: vec![
                        BondChange {
                            row: o_neighbor.bond.index(),
                            order: BondOrder::Double,
                            endpoints: None,
                        },
                        BondChange {
                            row: triple_bond,
                            order: BondOrder::Double,
                            endpoints: None,
                        },
                    ],
                    charge_changes: Vec::new(),
                    hydrogen_transfer: Some(transfer),
                });
            }
        }
    }
    Ok(candidates)
}

fn find_azide_mesomer(
    topology: &TopologyBlock,
    limit: usize,
) -> Result<Vec<RuleCandidate>, MsReadyFailure> {
    // BEGIN MSREADY MRV SOURCE: mrv/Tauto_Neut (#3652)__Metanode (#3412)__MarvinSketch (#2256).mrv; SHA-256 af6ad0c4d042311cac4834f00950fe96e53c8da53ae06c2768c78047f6291458
    // MSReady❗❗: source charge/bond roles are retained and concrete R attachment is preserved. Static MRV review cannot prove Indigo execution; the engine's canonical-per-match work has no executable source complexity baseline.
    // <atomArray><atom id="a1" elementType="N" formalCharge="-1" x2="-4.6599074823575135" y2="1.7601108531054774"></atom><atom id="a2" elementType="N" formalCharge="1" x2="-3.3300252099952097" y2="0.990179011211512"></atom><atom id="a3" elementType="N" x2="-1.9900638298844755" y2="1.7601108531054774"></atom><atom id="a4" elementType="R" rgroupRef="1" x2="-5.056912337559566" y2="3.2480857970033425" ligandOrder="a1"></atom></atomArray>
    // <bondArray><bond id="b1" atomRefs2="a1 a2" order="1"></bond><bond id="b2" atomRefs2="a2 a3" order="3"></bond><bond id="b3" atomRefs2="a1 a4" order="1"></bond></bondArray>
    // <bondArray><bond id="b1" atomRefs2="a1 a2" order="2"></bond><bond id="b2" atomRefs2="a2 a3" order="2"></bond><bond id="b3" atomRefs2="a1 a4" order="1"></bond></bondArray>
    // END MSREADY MRV SOURCE
    let valence = assign_for_match(topology)?;
    let mut candidates = Vec::new();
    for (n1, first) in topology.atoms.iter().enumerate() {
        if first.element() != Element::N
            || first.formal_charge() != -1
            || first.is_aromatic()
            || first.radical_electrons() != 0
            || first.explicit_hydrogens() != 0
            || !first.tracked_isotopic_hydrogens().is_empty()
            || hydrogen_count(topology, &valence, n1)? != 0
        {
            continue;
        }
        let n1_neighbors = topology.adjacency.neighbors_of(n1);
        let second_neighbors: Vec<_> = n1_neighbors
            .iter()
            .filter(|neighbor| topology.atoms[neighbor.atom_index].element() == Element::N)
            .collect();
        let r_neighbors: Vec<_> = n1_neighbors
            .iter()
            .filter(|neighbor| topology.atoms[neighbor.atom_index].element() != Element::N)
            .collect();
        if second_neighbors.len() != 1
            || r_neighbors.len() != 1
            || topology.atoms[r_neighbors[0].atom_index].atomic_number() == 1
        {
            continue;
        }
        let n2 = second_neighbors[0].atom_index;
        let r = r_neighbors[0].atom_index;
        let n1_n2_bond = &topology.bonds[second_neighbors[0].bond.index()];
        let n1_r_bond = &topology.bonds[r_neighbors[0].bond.index()];
        if n1_n2_bond.order() != BondOrder::Single
            || n1_n2_bond.query().is_some()
            || n1_n2_bond.is_aromatic()
            || n1_r_bond.order() != BondOrder::Single
            || n1_r_bond.query().is_some()
            || n1_r_bond.is_aromatic()
        {
            continue;
        }
        let second = &topology.atoms[n2];
        if second.formal_charge() != 1
            || second.is_aromatic()
            || second.radical_electrons() != 0
            || second.explicit_hydrogens() != 0
            || !second.tracked_isotopic_hydrogens().is_empty()
            || hydrogen_count(topology, &valence, n2)? != 0
        {
            continue;
        }
        let n2_neighbors = topology.adjacency.neighbors_of(n2);
        if n2_neighbors.len() != 2 {
            continue;
        }
        let Some(n3_neighbor) = n2_neighbors
            .iter()
            .find(|neighbor| neighbor.atom_index != n1)
        else {
            continue;
        };
        let n3 = n3_neighbor.atom_index;
        let third = &topology.atoms[n3];
        let n2_n3_bond = &topology.bonds[n3_neighbor.bond.index()];
        if third.element() != Element::N
            || third.formal_charge() != 0
            || third.is_aromatic()
            || third.radical_electrons() != 0
            || third.explicit_hydrogens() != 0
            || !third.tracked_isotopic_hydrogens().is_empty()
            || hydrogen_count(topology, &valence, n3)? != 0
            || topology.adjacency.neighbors_of(n3).len() != 1
            || n2_n3_bond.order() != BondOrder::Triple
            || n2_n3_bond.query().is_some()
            || n2_n3_bond.is_aromatic()
        {
            continue;
        }
        if candidates.len() >= limit {
            return Err(failure(
                Kind::TransformationLimit,
                format!("rule candidate count exceeds limit {limit}"),
            ));
        }
        candidates.push(RuleCandidate {
            rule: RuleId::AzideMesomer,
            atom_rows: vec![n1, n2, n3, r],
            bond_rows: vec![
                second_neighbors[0].bond.index(),
                n3_neighbor.bond.index(),
                r_neighbors[0].bond.index(),
            ],
            bond_changes: vec![
                BondChange {
                    row: second_neighbors[0].bond.index(),
                    order: BondOrder::Double,
                    endpoints: None,
                },
                BondChange {
                    row: n3_neighbor.bond.index(),
                    order: BondOrder::Double,
                    endpoints: None,
                },
            ],
            charge_changes: vec![
                ChargeChange { row: n1, charge: 0 },
                ChargeChange {
                    row: n3,
                    charge: -1,
                },
            ],
            hydrogen_transfer: None,
        });
    }
    Ok(candidates)
}

pub(crate) fn assign_for_match(
    topology: &TopologyBlock,
) -> Result<ValenceAssignment, MsReadyFailure> {
    assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, false).map_err(
        |error| MsReadyFailure {
            evidence: Default::default(),
            kind: Kind::OwnerError,
            stage: Stage::Normalization,
            context: format!("core valence assignment failed during rule matching: {error}"),
            original: None,
        },
    )
}

pub(crate) fn hydrogen_count(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    row: usize,
) -> Result<u32, MsReadyFailure> {
    let atom = topology.atoms.get(row).ok_or_else(|| {
        failure(
            Kind::InvariantViolation,
            format!("hydrogen count atom row {row} is out of range"),
        )
    })?;
    let implicit = valence
        .implicit_hydrogens
        .get(row)
        .copied()
        .ok_or_else(|| {
            failure(
                Kind::InvariantViolation,
                format!("core assignment lacks atom row {row}"),
            )
        })?;
    if implicit < 0 {
        return Err(failure(
            Kind::OwnerError,
            format!("negative inferred hydrogen count at atom row {row}"),
        ));
    }
    let hydrogen_neighbors = topology
        .adjacency
        .neighbors_of(row)
        .iter()
        .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 1)
        .count() as u32;
    let tracked = atom.tracked_isotopic_hydrogens().len() as u32;
    let ordinary_implicit = (implicit as u32).saturating_sub(tracked);
    Ok(u32::from(atom.explicit_hydrogens()) + hydrogen_neighbors + tracked + ordinary_implicit)
}

fn hydrogen_transfers(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    donor: usize,
    acceptor: usize,
) -> Result<Vec<HydrogenTransfer>, MsReadyFailure> {
    let mut carriers = Vec::new();
    for neighbor in topology.adjacency.neighbors_of(donor) {
        if topology.atoms[neighbor.atom_index].atomic_number() == 1 {
            let bond = &topology.bonds[neighbor.bond.index()];
            if bond.order() != BondOrder::Single || bond.query().is_some() || bond.is_aromatic() {
                continue;
            }
            carriers.push(HydrogenTransfer {
                donor,
                acceptor,
                carrier: HydrogenCarrier::ExplicitAtom {
                    atom: neighbor.atom_index,
                    bond: neighbor.bond.index(),
                    isotope: topology.atoms[neighbor.atom_index].isotope(),
                },
            });
        }
    }
    if topology.atoms[donor].explicit_hydrogens() != 0 {
        carriers.push(HydrogenTransfer {
            donor,
            acceptor,
            carrier: HydrogenCarrier::Bracket,
        });
    }
    for (index, isotope) in topology.atoms[donor]
        .tracked_isotopic_hydrogens()
        .iter()
        .copied()
        .enumerate()
    {
        carriers.push(HydrogenTransfer {
            donor,
            acceptor,
            carrier: HydrogenCarrier::TrackedIsotope { index, isotope },
        });
    }
    if u32::try_from(valence.implicit_hydrogens.get(donor).copied().unwrap_or(0))
        .unwrap_or_default()
        > topology.atoms[donor].tracked_isotopic_hydrogens().len() as u32
    {
        carriers.push(HydrogenTransfer {
            donor,
            acceptor,
            carrier: HydrogenCarrier::Implicit,
        });
    }
    carriers.sort_by_key(|transfer| match transfer.carrier {
        HydrogenCarrier::ExplicitAtom { isotope, .. } => (0u8, isotope.unwrap_or(0)),
        HydrogenCarrier::Bracket => (1, 0),
        HydrogenCarrier::TrackedIsotope { isotope, .. } => (2, isotope),
        HydrogenCarrier::Implicit => (3, 0),
    });
    carriers.dedup_by_key(|transfer| match transfer.carrier {
        HydrogenCarrier::ExplicitAtom { isotope, .. } => (0u8, isotope.unwrap_or(0)),
        HydrogenCarrier::Bracket => (1, 0),
        HydrogenCarrier::TrackedIsotope { isotope, .. } => (2, isotope),
        HydrogenCarrier::Implicit => (3, 0),
    });
    Ok(carriers)
}

fn attachment_slot_count(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: usize,
    core_neighbors: &[usize],
) -> Result<usize, MsReadyFailure> {
    // Native R-slot adaptation: preserve concrete branches and count owner H
    // (implicit, bracket or explicit) once as available schematic attachments.
    Ok(topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .filter(|neighbor| {
            !core_neighbors.contains(&neighbor.atom_index)
                && topology.atoms[neighbor.atom_index].atomic_number() != 1
        })
        .count()
        + hydrogen_count(topology, valence, atom)? as usize)
}

fn single_attachments(topology: &TopologyBlock, atom: usize, core_neighbors: &[usize]) -> bool {
    topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .all(|neighbor| {
            core_neighbors.contains(&neighbor.atom_index)
                || (topology.bonds[neighbor.bond.index()].order() == BondOrder::Single
                    && !topology.bonds[neighbor.bond.index()].is_aromatic()
                    && topology.bonds[neighbor.bond.index()].query().is_none())
        })
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
mod rule_tests {
    use super::*;
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

    fn parse(smiles: &str) -> TopologyBlock {
        parse_smiles(smiles, &SmilesParseParams::default())
            .unwrap_or_else(|error| panic!("cannot parse {smiles}: {error}"))
            .topology
    }

    #[test]
    fn each_frozen_rule_appearance_has_a_structural_match() {
        let examples = [
            ("CN(=O)=O", RuleId::NitroInitial),
            ("CC(C)=C(NC)C", RuleId::EnamineToImine),
            ("CC(C)=C(O)C", RuleId::EnolToCarbonyl),
            ("CC#CO", RuleId::YnolToKetene),
            ("C[N-][N+]#N", RuleId::AzideMesomer),
            ("CN(=O)=O", RuleId::NitroFinal),
        ];
        for (smiles, rule) in examples {
            let candidates = find_candidates(&parse(smiles), rule)
                .unwrap_or_else(|error| panic!("{} match failed: {error}", rule.id()));
            assert!(
                !candidates.is_empty(),
                "{} did not match {smiles}",
                rule.id()
            );
        }
    }

    #[test]
    fn hydrogen_carriers_cover_implicit_bracket_explicit_and_tracked_isotope_forms() {
        let implicit = find_candidates(&parse("CC(C)=C(O)C(C)C"), RuleId::EnolToCarbonyl).unwrap();
        assert!(implicit.iter().any(|candidate| matches!(
            candidate.hydrogen_transfer.map(|transfer| transfer.carrier),
            Some(HydrogenCarrier::Implicit)
        )));

        let bracket =
            find_candidates(&parse("CC(C)=C([OH])C(C)C"), RuleId::EnolToCarbonyl).unwrap();
        assert!(bracket.iter().any(|candidate| matches!(
            candidate.hydrogen_transfer.map(|transfer| transfer.carrier),
            Some(HydrogenCarrier::Bracket)
        )));

        let explicit =
            find_candidates(&parse("CC(C)=C(O[H])C(C)C"), RuleId::EnolToCarbonyl).unwrap();
        assert!(explicit.iter().any(|candidate| matches!(
            candidate.hydrogen_transfer.map(|transfer| transfer.carrier),
            Some(HydrogenCarrier::ExplicitAtom { isotope: None, .. })
        )));
        for (isotope, smiles) in [(2, "CC(C)=C(O[2H])C(C)C"), (3, "CC(C)=C(O[3H])C(C)C")] {
            let isotopic = find_candidates(&parse(smiles), RuleId::EnolToCarbonyl).unwrap();
            assert!(isotopic.iter().any(|candidate| matches!(
                candidate.hydrogen_transfer.map(|transfer| transfer.carrier),
                Some(HydrogenCarrier::ExplicitAtom { isotope: Some(found), .. }) if found == isotope
            )));
        }

        for isotope in [2, 3] {
            let mut tracked_topology = parse("CC(C)=C(O)C(C)C");
            let oxygen = tracked_topology
                .atoms
                .iter()
                .position(|atom| atom.element() == Element::O)
                .unwrap();
            tracked_topology.atoms[oxygen].set_no_implicit(true);
            tracked_topology.atoms[oxygen].set_tracked_isotopic_hydrogens(vec![isotope]);
            let tracked = find_candidates(&tracked_topology, RuleId::EnolToCarbonyl).unwrap();
            assert!(tracked.iter().any(|candidate| matches!(
                candidate.hydrogen_transfer.map(|transfer| transfer.carrier),
                Some(HydrogenCarrier::TrackedIsotope { isotope: found, .. }) if found == isotope
            )));
        }
    }

    #[test]
    fn packaged_six_appearance_catalog_matches_frozen_order_and_source_hashes() {
        validate_rule_data().unwrap();
        let data: serde_json::Value = serde_json::from_str(RULE_DATA).unwrap();
        let appearances = data["appearances"].as_array().unwrap();
        assert_eq!(appearances.len(), 6);
        assert_eq!(appearances[0]["appearance_id"], "initial_nitro:03:03");
        assert_eq!(appearances[5]["appearance_id"], "final_nitro:03:03");
        assert_eq!(
            appearances[0]["source_entry"],
            appearances[5]["source_entry"]
        );
        assert_ne!(
            appearances[0]["source_sha256"],
            appearances[5]["source_sha256"]
        );
    }
}
