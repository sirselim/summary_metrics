use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

fn calculate_n50(lengths: &mut Vec<usize>) -> usize {
    lengths.sort_unstable_by(|a, b| b.cmp(a));
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

fn main() {
    // Parse command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input_file> <read_length>", &args[0]);
        std::process::exit(1);
    }
    let input_file = &args[1];
    let _read_len: usize = args[2].parse().expect("Invalid read length");

    // Read input file and calculate statistics
    let file = File::open(input_file).expect("Failed to open input file");
    let reader = BufReader::new(file);

    let mut total_gigabases = 0.0;
    let mut total_gigabases_with_barcode = 0.0;
    let mut pass_reads = 0;
    let mut pass_with_barcode = 0;
    let mut gigabases_with_barcode_and_length = 0.0;
    let mut most_prevalent_barcode = HashMap::new();
    let mut lengths = Vec::with_capacity(100000); // Pre-allocate space for the vector

    for line in reader.lines().skip(1) {
        let line = line.expect("Failed to read line");
        let fields: Vec<&str> = line.split('\t').collect();
        let barcode = fields[26];
        let length: usize = fields[15].parse().expect("Invalid read length");
        let pass = fields[11] == "TRUE";

        if pass {
            total_gigabases += length as f64 / 1_000_000_000.0;
            lengths.push(length);
            pass_reads += 1;

            let entry = most_prevalent_barcode.entry(barcode.to_string()).or_insert(0);
            *entry += 1;
            if *entry > pass_with_barcode {
                pass_with_barcode = *entry;
                total_gigabases_with_barcode = total_gigabases_with_barcode + length as f64 / 1_000_000_000.0;
            }

            if length >= 15000 && *entry == pass_with_barcode {
                gigabases_with_barcode_and_length += length as f64 / 1_000_000_000.0;
            }
        }
    }

    let (most_prevalent_barcode, most_prevalent_count) = most_prevalent_barcode.iter().max_by_key(|(_, &count)| count)
        .expect("No barcodes found");

    println!("Most prevalent barcode: {} (Count: {})", most_prevalent_barcode, most_prevalent_count);
    println!("Number of reads that pass: {}", pass_reads);
    println!("Number of reads that pass with the detected barcode: {}", pass_with_barcode);
    println!("Total gigabases of reads that pass: {:.2} Gb", total_gigabases);
    println!("Total gigabases of reads that pass with the detected barcode: {:.2} Gb", total_gigabases_with_barcode);
    println!("Total gigabases of reads that pass with the detected barcode and are >= 15000bp: {:.2} Gb", gigabases_with_barcode_and_length);
    println!("N50: {:.2} Kb", calculate_n50(&mut lengths) as f64 / 1000.0);
}