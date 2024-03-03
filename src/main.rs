use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;

fn calculate_n50(lengths: &mut Vec<usize>) -> usize {
    // Sort the lengths vector using radix sort
    lengths.sort();

    let total_bases: usize = lengths.iter().sum();
    let half_total_bases = total_bases / 2;
    let mut cumulative_bases = 0;
    for length in lengths {
        cumulative_bases += *length;
        if cumulative_bases >= half_total_bases {
            return *length;
        }
    }
    0 // N50 not found
}

fn print_help() {
    println!("Usage: <input_file> <minimum_length_threshold>");
    println!("<input_file>                  Nanopore sequencing summary text file");
    println!("<minimum_length_threshold>    length to filter for statistics");
    println!("Options:");
    println!("  -h, --help                  Print this help message");
    println!("  -v, --version               Print version information");
    println!("  --json                      Output result to JSON file");
}

fn print_version() {
    println!("summary_metrics version {}", env!("CARGO_PKG_VERSION"));
}

fn main() {
    // Parse command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 && !args.contains(&"--json".to_string()) {
        match args.len() {
            2 if (args[1] == "-h" || args[1] == "--help") => {
                print_help();
                return;
            }
            2 if (args[1] == "-v" || args[1] == "--version") => {
                print_version();
                return;
            }
            _ => {
                eprintln!("Invalid number of arguments. Use -h or --help for usage information.");
                std::process::exit(1);
            }
        }
    }
    let input_file = &args[1];
    let min_length_threshold: usize = args[2].parse().expect("Invalid minimum length threshold");

    // Open the file and read the header
    let file = File::open(input_file).expect("Failed to open input file");
    let cloned_file = file.try_clone().expect("Failed to clone file");
    let reader = BufReader::new(cloned_file);
    let header_line = reader.lines().next().expect("Failed to read header").unwrap();
    let headers: Vec<&str> = header_line.split('\t').collect();

    // Find column indices based on headers
    let passes_filtering_index = headers.iter().position(|&x| x == "passes_filtering").unwrap();
    let sequence_length_index = headers.iter().position(|&x| x == "sequence_length_template").unwrap();
    let barcode_arrangement_index = headers.iter().position(|&x| x == "barcode_arrangement").unwrap();

    // Process the rest of the file
    let mut total_gigabases = 0.0;
    let mut total_gigabases_with_barcode = 0.0;
    let mut pass_reads = 0;
    let mut pass_with_barcode = 0;
    let mut gigabases_with_barcode_and_length = 0.0;
    let mut most_prevalent_barcode = HashMap::new();
    let mut lengths = Vec::new();

    let mut output_to_json = false;
    if args.contains(&"--json".to_string()) {
        output_to_json = true;
    }

    let cloned_reader = BufReader::new(file.try_clone().expect("Failed to clone reader"));

    for line in cloned_reader.lines().skip(1) {
        let line = line.expect("Failed to read line");
        let fields: Vec<&str> = line.split('\t').collect();

        let passes_filtering = fields[passes_filtering_index] == "TRUE";
        let sequence_length: usize = fields[sequence_length_index].parse().expect("Invalid sequence length");
        let barcode_arrangement = fields[barcode_arrangement_index];

        if passes_filtering {
            total_gigabases += sequence_length as f64 / 1_000_000_000.0;
            lengths.push(sequence_length);
            pass_reads += 1;

            let entry = most_prevalent_barcode.entry(barcode_arrangement.to_string()).or_insert(0);
            *entry += 1;
            if *entry > pass_with_barcode {
                pass_with_barcode = *entry;
                total_gigabases_with_barcode += sequence_length as f64 / 1_000_000_000.0;
            }

            if sequence_length >= min_length_threshold && *entry == pass_with_barcode {
                gigabases_with_barcode_and_length += sequence_length as f64 / 1_000_000_000.0;
            }
        }
    }

    let (most_prevalent_barcode, most_prevalent_count) = most_prevalent_barcode.iter().max_by_key(|(_, &count)| count)
        .expect("No barcodes found");

    if output_to_json {
        // Output to JSON file
        let mut json_output = File::create("output.json").expect("Failed to create JSON file");
        let json_data = serde_json::json!({
            "most_prevalent_barcode": most_prevalent_barcode,
            "most_prevalent_count": most_prevalent_count,
            "pass_reads": pass_reads,
            "pass_with_barcode": pass_with_barcode,
            "total_gigabases": total_gigabases,
            "total_gigabases_with_barcode": total_gigabases_with_barcode,
            "gigabases_with_barcode_and_length": gigabases_with_barcode_and_length,
            "n50": calculate_n50(&mut lengths) as f64 / 1000.0
        });
        writeln!(json_output, "{}", json_data.to_string()).expect("Failed to write JSON data to file");
    } else {
        // Output to console
        // println!("Most prevalent barcode: {} (Count: {})", most_prevalent_barcode, most_prevalent_count);
        println!("Number of reads that pass: {}", pass_reads);
        println!("Number of reads that pass with the detected barcode: {}", pass_with_barcode);
        println!("Total gigabases of reads that pass: {:.2} Gb", total_gigabases);
        println!("Total gigabases of reads that pass with the detected barcode: {:.2} Gb", total_gigabases_with_barcode);
        println!("Total gigabases of reads that pass with the detected barcode and are >= {} bp: {:.2} Gb", min_length_threshold, gigabases_with_barcode_and_length);
        println!("N50: {:.2} Kb", calculate_n50(&mut lengths) as f64 / 1000.0);
    }
}
