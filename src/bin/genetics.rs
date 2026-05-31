use rosalind_rs::simulation::genetics::AllelePopulation;
use std::env;
use std::fs;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <input_file_path>", args[0]);
        process::exit(1);
    }

    let file_path = &args[1];
    let content = match fs::read_to_string(file_path) {
        Ok(c) => c,
        Err(err) => {
            eprintln!("Error reading file '{}': {}", file_path, err);
            process::exit(1);
        }
    };

    // Extract the space-separated numbers from the input file
    let tokens: Vec<&str> = content.split_whitespace().collect();
    if tokens.len() < 3 {
        eprintln!(
            "Error: Input file must contain three integers (dominant, heterozygous, recessive)"
        );
        process::exit(1);
    }

    let homozygous_dominant = parse_population_count(tokens[0], "homozygous dominant");
    let heterozygous = parse_population_count(tokens[1], "heterozygous");
    let homozygous_recessive = parse_population_count(tokens[2], "homozygous recessive");

    let population = AllelePopulation::new(homozygous_dominant, heterozygous, homozygous_recessive);
    let probability = population.dominant_phenotype_probability();

    // Output formatted to 5 decimal places
    println!("{:.5}", probability);
}

/// Helper function to parse input tokens safely with clear error context
fn parse_population_count(token: &str, label: &str) -> u64 {
    match token.parse::<u64>() {
        Ok(count) => count,
        Err(_) => {
            eprintln!(
                "Error: Failed to parse {} count from token '{}'",
                label, token
            );
            process::exit(1);
        }
    }
}
