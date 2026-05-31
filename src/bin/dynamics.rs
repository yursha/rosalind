use rosalind_rs::simulation::dynamics::{AgeStructuredModel, PopulationState};
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

    let tokens: Vec<&str> = content.split_whitespace().collect();
    if tokens.len() < 4 {
        eprintln!(
            "Error: Input file must contain exactly four integers:\n\
             <elapsed_intervals> <fecundity> <initial_juveniles> <initial_adults>"
        );
        process::exit(1);
    }

    let elapsed_intervals = parse_u64(tokens[0], "elapsed time intervals") as u32;
    let fecundity = parse_u64(tokens[1], "fecundity");
    let initial_juveniles = parse_u64(tokens[2], "initial juveniles");
    let initial_adults = parse_u64(tokens[3], "initial adults");

    let model = AgeStructuredModel::new(fecundity, 1.0);
    let initial_state = PopulationState::new(initial_juveniles, initial_adults);

    let final_state = model.project(initial_state, elapsed_intervals);

    println!("{}", final_state.total());
}

fn parse_u64(token: &str, label: &str) -> u64 {
    match token.parse::<u64>() {
        Ok(count) => count,
        Err(_) => {
            eprintln!("Error: Failed to parse {} from token '{}'", label, token);
            process::exit(1);
        }
    }
}
