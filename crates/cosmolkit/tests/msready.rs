use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingReceiver, BindingTypeRole, FunctionStatus, Molecule,
    MsReadyFailure, MsReadyFailureKind, MsReadyParams, MsReadyResult, MsReadyStatus, StateModel,
};

#[test]
fn canonical_signatures_defaults_and_native_registration() {
    let _: fn(&Molecule) -> Result<MsReadyResult, MsReadyFailure> = Molecule::with_msready;
    let _: fn(&Molecule, &MsReadyParams) -> Result<MsReadyResult, MsReadyFailure> =
        Molecule::with_msready_with_params;
    assert_eq!(MsReadyParams::default().max_edits, 1024);
    assert_eq!(MsReadyParams::default().max_matches_per_sweep, 4096);
    assert_eq!(MsReadyParams::default().record_id, None);
    assert_eq!(MsReadyParams::default().original_smiles, None);
    for (name, js, count) in [
        ("with_msready", "withMsready", 0),
        ("with_msready_with_params", "withMsreadyWithParams", 1),
    ] {
        let id = format!("Molecule.{name}");
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == id)
            .unwrap();
        assert_eq!(row.feature, "cap-msready");
        assert_eq!(row.status, FunctionStatus::Native);
        assert_eq!(row.python_name, name);
        assert_eq!(row.javascript_name, js);
        let callable = row.callable.unwrap();
        assert_eq!(callable.receiver, Some(BindingReceiver::Shared));
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(callable.parameters.len(), count);
        if count == 1 {
            assert_eq!(callable.parameters[0].default, BindingDefault::Required);
        }
        assert!(callable.output_type.contains("MsReadyResult"));
        assert!(callable.error_type.unwrap().contains("MsReadyFailure"));
    }
    for (name, role) in [
        ("MsReadyParams", BindingTypeRole::Parameter),
        ("MsReadyFailure", BindingTypeRole::Error),
        ("MsReadyResult", BindingTypeRole::Result),
    ] {
        let id = format!("types.{name}");
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == id)
            .unwrap();
        assert_eq!(row.type_role, Some(role));
    }
    let rows: Vec<_> = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "cap-msready")
        .collect();
    assert_eq!(rows.len(), 26); // all 24 nested vocabulary types plus two methods
    assert!(rows.iter().all(|row| row.status == FunctionStatus::Native));
}

#[test]
fn defaults_forward_to_owner_and_return_detached_values_without_changing_input() {
    let molecule = Molecule::from_smiles("CCO.CCO.CN").unwrap();
    let peer = molecule.clone();
    let topology = molecule.topology().clone();
    let properties = molecule.properties().clone();
    let result = molecule.with_msready().unwrap();
    assert_eq!(
        result,
        molecule
            .with_msready_with_params(&MsReadyParams::default())
            .unwrap()
    );
    assert_eq!(
        result,
        cosmolkit_msready::standardize(&topology, &properties, &MsReadyParams::default()).unwrap()
    );
    let MsReadyStatus::Products(products) = result.status else {
        panic!("expected products")
    };
    assert_eq!(products.len(), 2);
    assert_eq!(
        products
            .iter()
            .map(|product| product.multiplicity)
            .sum::<usize>(),
        3
    );
    for value in [&molecule, &peer] {
        assert!(std::ptr::eq(value.topology(), molecule.topology()));
        assert!(std::ptr::eq(value.properties(), molecule.properties()));
        assert_eq!(value.topology(), &topology);
        assert_eq!(value.properties(), &properties);
    }
}

#[test]
fn explicit_metadata_and_failure_evidence_survive_the_facade() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let peer = molecule.clone();
    let params = MsReadyParams {
        record_id: Some("sample-42".into()),
        original_smiles: Some("archived text".into()),
        ..Default::default()
    };
    let result = molecule.with_msready_with_params(&params).unwrap();
    assert_eq!(result.input_link.record_id, params.record_id);
    assert_eq!(result.input_link.original_smiles, params.original_smiles);
    assert_eq!(result.input_link.source_smiles, None); // this is the graph boundary
    let failure = molecule
        .with_msready_with_params(&MsReadyParams {
            max_edits: 0,
            ..params
        })
        .unwrap_err();
    assert_eq!(failure.kind, MsReadyFailureKind::InvalidInput);
    assert_eq!(failure.original.unwrap().topology, *molecule.topology());
    assert_eq!(
        failure.evidence.input_link.record_id.as_deref(),
        Some("sample-42")
    );
    assert!(std::ptr::eq(molecule.topology(), peer.topology()));
    assert!(std::ptr::eq(molecule.properties(), peer.properties()));
    let _: &(dyn std::error::Error + 'static) = &MsReadyFailure {
        kind: failure.kind,
        stage: failure.stage,
        context: failure.context,
        original: None,
        evidence: failure.evidence,
    };
}
