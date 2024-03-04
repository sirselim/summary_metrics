use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead};

fn calculate_n50(sequence_lengths: &mut Vec<u64>) -> Option<u64> {
    sequence_lengths.sort(); // Sort sequence lengths in ascending order
    let total_sum: u64 = sequence_lengths.iter().sum(); // Calculate total sum of sequence lengths
    let midpoint = total_sum / 2; // Calculate midpoint

    let mut running_sum = 0;
    for length in sequence_lengths.iter().rev() {
        running_sum += *length;
        if running_sum >= midpoint {
            return Some(*length);
        }
    }
    None
}

fn print_help() {
    println!("Usage: <input_file> <minimum_length_threshold>");
    println!("<input_file>                  Nanopore sequencing summary text file");
    println!("<minimum_length_threshold>    length to filter for statistics");
    println!("Options:");
    println!("  -h, --help                  Print this help message");
    println!("  -v, --version               Print version information");
}

fn print_version() {
    println!("summary_metrics version {}", env!("CARGO_PKG_VERSION"));
}

fn main() -> io::Result<()> {
    // Get the filename and minimum length from command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        match args.len() {
            2 if (args[1] == "-h" || args[1] == "--help") => {
                print_help();
                return Ok(());
            }
            2 if (args[1] == "-v" || args[1] == "--version") => {
                print_version();
                return Ok(());
            }
            _ => {
                eprintln!("Invalid number of arguments. Use -h or --help for usage information.");
                std::process::exit(1);
            }
        }
    }
    let filename = &args[1];
    let min_length: u64 = args[2].parse().expect("Minimum length must be a positive integer");

    // Open the file
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);

    // Count the lines
    let mut total_line_count = 0;
    let mut passes_filtering_index = None;
    let mut barcode_index = None;
    let mut sequence_length_index = None;

    let mut total_barcode_counts: HashMap<String, u32> = HashMap::new();
    let mut passing_barcode_counts: HashMap<String, u32> = HashMap::new();
    let mut total_sequence_data_all: u64 = 0; // Total sequence data for all reads
    let mut total_sequence_data_passed: u64 = 0; // Total sequence data for passed reads
    let mut total_sequence_data_most_common_barcode: u64 = 0; // Total sequence data for passed reads with most common barcode
    let mut long_sequences_most_common_barcode: u64 = 0; // Total sequence data longer than the minimum length for the most common barcode

    let mut sequence_lengths: Vec<u64> = Vec::new(); // Store sequence lengths for N50 calculation

    // Iterate over each line in the file
    for (index, line) in reader.lines().enumerate() {
        let line = line?; // Unwrap the line
        let fields: Vec<&str> = line.split('\t').collect(); // Split the line by tab

        // If it's the first line, process the header
        if index == 0 {
            // Find the index of 'passes_filtering' column in the header
            passes_filtering_index = fields.iter().position(|&x| x == "passes_filtering");
            // Find the index of 'barcode_arrangement' column in the header
            barcode_index = fields.iter().position(|&x| x == "barcode_arrangement");
            // Find the index of 'sequence_length_template' column in the header
            sequence_length_index = fields.iter().position(|&x| x == "sequence_length_template");

            if passes_filtering_index.is_none() || barcode_index.is_none() || sequence_length_index.is_none() {
                eprintln!("One or more columns not found in the header.");
                std::process::exit(1);
            }

            continue; // Skip header
        }

        // Increment total line count
        total_line_count += 1;

        // Count occurrences of the barcode arrangement in all reads
        if let Some(barcode_idx) = barcode_index {
            if let Some(barcode_arrangement) = fields.get(barcode_idx) {
                *total_barcode_counts.entry(barcode_arrangement.to_string()).or_insert(0) += 1;
            }
        }

        // Add sequence length to total sequence data for all reads and collect for N50 calculation
        if let Some(seq_len_idx) = sequence_length_index {
            if let Some(seq_len_str) = fields.get(seq_len_idx) {
                if let Ok(seq_len) = seq_len_str.parse::<u64>() {
                    total_sequence_data_all += seq_len;
                    sequence_lengths.push(seq_len);
                }
            }
        }

        // Check if the line passes the filter
        if let Some(passes_filtering_idx) = passes_filtering_index {
            if let Some(passes_filtering) = fields.get(passes_filtering_idx) {
                if *passes_filtering == "TRUE" {
                    // Count occurrences of the barcode arrangement in passed reads
                    if let Some(barcode_idx) = barcode_index {
                        if let Some(barcode_arrangement) = fields.get(barcode_idx) {
                            *passing_barcode_counts.entry(barcode_arrangement.to_string()).or_insert(0) += 1;
                        }
                    }

                    // Add sequence length to total sequence data for passed reads
                    if let Some(seq_len_idx) = sequence_length_index {
                        if let Some(seq_len_str) = fields.get(seq_len_idx) {
                            if let Ok(seq_len) = seq_len_str.parse::<u64>() {
                                total_sequence_data_passed += seq_len;

                                // Check if barcode is detected and matches the most common one
                                if let Some((most_common_barcode, _)) = passing_barcode_counts.iter().max_by_key(|&(_, count)| count) {
                                    if let Some(barcode_idx) = barcode_index {
                                        if let Some(barcode_arrangement) = fields.get(barcode_idx) {
                                            if barcode_arrangement == &most_common_barcode {
                                                total_sequence_data_most_common_barcode += seq_len;
                                                if seq_len > min_length {
                                                    long_sequences_most_common_barcode += seq_len;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Calculate and print N50
    if let Some(n50) = calculate_n50(&mut sequence_lengths) {
        // Convert N50 from bases to kilobases
        let n50_kb = n50 as f64 / 1000.0;
        println!("N50 (total): {:.2} Kb", n50_kb);
    } else {
        println!("No sequence length data found.");
    }

    // Print total number of reads
    println!("Total reads: {}", total_line_count);

    // Find and print the most common barcode in all reads
    if let Some((max_barcode, max_count)) = total_barcode_counts.iter().max_by_key(|&(_, count)| count) {
        println!("Detected barcode (total): {} (count: {})", max_barcode, max_count);
    } else {
        println!("No barcode data found in the file.");
    }

    // Print total number of reads passing the filter
    println!("Total passed reads: {}", passing_barcode_counts.values().sum::<u32>());

    // Find and print the most common barcode in passed reads
    if let Some((max_barcode, max_count)) = passing_barcode_counts.iter().max_by_key(|&(_, count)| count) {
        println!("Detected barcode (passed): {} (count: {})", max_barcode, max_count);
    } else {
        println!("No barcode data found in the passed reads.");
    }

    // Convert total sequence data from bases to gigabases
    let total_sequence_data_all_gb = total_sequence_data_all as f64 / 1_000_000_000.0;
    let total_sequence_data_passed_gb = total_sequence_data_passed as f64 / 1_000_000_000.0;
    let total_sequence_data_most_common_barcode_gb = total_sequence_data_most_common_barcode as f64 / 1_000_000_000.0;
    let long_sequences_most_common_barcode_gb = long_sequences_most_common_barcode as f64 / 1_000_000_000.0;

    // Print total amount of sequence data in gigabases for all reads, passed reads, and passed reads with most common barcode
    println!("Total output: {:.2} gigabases", total_sequence_data_all_gb);
    println!("Total output (passed): {:.2} gigabases", total_sequence_data_passed_gb);
    println!("Total output (passed, barcode): {:.2} gigabases", total_sequence_data_most_common_barcode_gb);
    println!("Total >= {} bases (passed, barcode): {:.2} gigabases", min_length, long_sequences_most_common_barcode_gb);

    Ok(())
}
