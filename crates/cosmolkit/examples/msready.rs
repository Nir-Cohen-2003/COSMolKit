use cosmolkit::{Molecule, MsReadyParams, MsReadyStatus};

// cargo run -p cosmolkit --no-default-features --features cap-smiles,cap-msready --example msready
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Organic counterions survive; identical acetate origins aggregate.
    let molecule = Molecule::from_smiles("CC(=O)[O-].[Na+].C[NH3+].CC(=O)[O-]")?;
    let result = molecule.with_msready_with_params(&MsReadyParams {
        record_id: Some("counterions-and-multiplicity".into()),
        ..Default::default()
    })?;
    println!(
        "{} {} ({})",
        result.run_info.profile_id, result.run_info.profile_version, result.run_info.behavior
    );
    match result.status {
        MsReadyStatus::Products(products) => {
            for product in products {
                println!(
                    "{} {} {} mass={:.6} origins={}",
                    product.canonical_smiles,
                    product.inchi_key,
                    product.formula,
                    product.neutral_atom_mass,
                    product.multiplicity,
                );
            }
        }
        MsReadyStatus::NoEligibleComponent => println!("No eligible component"),
    }
    println!(
        "exclusions={} events={}",
        result.exclusions.len(),
        result.events.len()
    );
    println!("Original molecule still has {} atoms", molecule.num_atoms());
    Ok(())
}
