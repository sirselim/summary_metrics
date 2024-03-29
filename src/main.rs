use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead};
use colored::*;
use clap::{App, Arg};
use num_format::{Locale, ToFormattedString};
use prettytable::{Table, row, format};
use std::path::{Path, PathBuf};

// function to locate and extract flowcell id from different locations
fn extract_ids(filename: &str, column_names: Vec<&str>) -> Option<Vec<Option<String>>> {
    let file = File::open(filename).ok()?;
    let mut reader = io::BufReader::new(file);

    // read the header line to find the column indices
    let mut header_line = String::new();
    reader.read_line(&mut header_line).ok()?;

    let columns: Vec<&str> = header_line.trim().split('\t').collect();

    // find the indices of the desired columns
    let mut column_indices = Vec::new();
    for &column_name in &column_names {
        if let Some(column_index) = columns.iter().position(|&x| x == column_name) {
            column_indices.push(column_index);
        } else {
            return None; // Column name not found
        }
    }

    // process the first line of data using the found column indices
    let mut first_data_line = String::new();
    reader.read_line(&mut first_data_line).ok()?;

    let columns_in_first_line: Vec<&str> = first_data_line.trim().split('\t').collect();

    let mut extracted_ids = vec![None; column_indices.len()];
    let mut found_flowcell_id = false;

    // checking if there are valid values to return, if missing we can look elsewhere
    for (i, &column_index) in column_indices.iter().enumerate() {
        if let Some(column_data) = columns_in_first_line.get(column_index) {
            if column_data != &"-" {
                extracted_ids[i] = Some(column_data.to_string());
                // If a valid ID is found, mark it as found
                if column_names[i] == "filename_pod5" {
                    found_flowcell_id = true;
                }
            }
        }
    }

    // if flowcell ID not found in any of the columns, extract from the filename
    if !found_flowcell_id {
        let filename_components: Vec<&str> = filename.split('/').last()?.split('_').collect();
        if let Some(flowcell_id) = filename_components.get(2) {
            let extracted_flowcell_id = flowcell_id.to_string();
            extracted_ids[0] = Some(extracted_flowcell_id);
        }
    }

    Some(extracted_ids)
}

// function to calculate N50
fn calculate_n50(sequence_lengths: &[u64]) -> Option<u64> {
    let mut sorted_lengths = sequence_lengths.to_vec();
    sorted_lengths.sort();
    let total_sum: u64 = sequence_lengths.iter().sum();
    let midpoint = total_sum / 2;

    let mut running_sum = 0;
    for &length in sorted_lengths.iter().rev() {
        running_sum += length;
        if running_sum >= midpoint {
            return Some(length);
        }
    }
    None
}

// function to calculate mean
fn calculate_mean(read_lengths: &[u64]) -> f64 {
    if read_lengths.is_empty() {
        return 0.0;
    }
    let sum: u64 = read_lengths.iter().sum();
    sum as f64 / read_lengths.len() as f64
}

// function to calculate median
fn calculate_median(read_lengths: &mut [u64]) -> f64 {
    read_lengths.sort();
    let mid = read_lengths.len() / 2;
    if read_lengths.len() % 2 == 0 {
        (read_lengths[mid - 1] + read_lengths[mid]) as f64 / 2.0
    } else {
        read_lengths[mid] as f64
    }
}

// function to process all identified barcodes
fn barcode(filename: &str, qscore_threshold: f64, length_threshold: u64, output_csv: bool) -> io::Result<()> {
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);

    let mut barcode_counts: HashMap<String, (u32, u32, u64, u64, u64, f64)> = HashMap::new(); // Modify to include an additional u32 for total reads before filtering

    let mut barcode_index = None;
    let mut mean_qscore_index = None;
    let mut sequence_length_index = None;

    // Iterate over each line in the file
    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if index == 0 {
            // Find the indices of required columns
            mean_qscore_index = fields.iter().position(|&x| x == "mean_qscore_template");
            sequence_length_index = fields.iter().position(|&x| x == "sequence_length_template");
            barcode_index = fields.iter().position(|&x| x == "barcode_arrangement");
            continue;
        }

        // Ensure required indices are found
        if let (Some(barcode_idx), Some(mean_qscore_idx), Some(seq_len_idx)) =
            (barcode_index, mean_qscore_index, sequence_length_index)
        {
            if let (Some(barcode_arrangement), Some(passes_filtering), Some(seq_len_str)) =
                (fields.get(barcode_idx), fields.get(mean_qscore_idx), fields.get(seq_len_idx))
            {
                // Convert sequence length and qscore to appropriate types
                let seq_len: u64 = seq_len_str.parse().unwrap_or(0);
                let qscore: f64 = passes_filtering.parse().unwrap_or(0.0);

                // If qscore meets threshold, update barcode counts
                if qscore >= qscore_threshold {
                    let (reads, total_reads, bases, passed_bases, passed_bases_filtered, passed_bases_gb) = barcode_counts
                        .entry(barcode_arrangement.to_string())
                        .or_insert((0, 0, 0, 0, 0, 0.0)); // Modify to include an additional 0 for total reads before filtering
                    *reads += 1;
                    *total_reads += 1; // Increment total reads before filtering
                    *bases += seq_len;
                    *passed_bases += seq_len;
                    if seq_len >= length_threshold {
                        *passed_bases_filtered += seq_len;
                        *passed_bases_gb = *passed_bases_filtered as f64 / 1_000_000_000.0;
                    }
                } else {
                    // If qscore doesn't meet threshold, update total bases only
                    if let Some((_, total_reads, total_bases, _, _, _)) = barcode_counts.get_mut(&barcode_arrangement[..]) {
                        *total_reads += 1; // Increment total reads before filtering
                        *total_bases += seq_len;
                    }
                }
            }
        }
    }

    // Sort barcodes by name
    let mut sorted_barcodes: Vec<_> = barcode_counts.into_iter().collect();
    sorted_barcodes.sort_by(|(barcode1, _), (barcode2, _)| barcode1.cmp(barcode2));

    // Print the results in a table format - using pretty table crate
    let mut table = Table::new();
    // let format = format::FormatBuilder::new()
    //     .column_separator('|')
    //     .borders('|')
    //     .separators(&[format::LinePosition::Top,
    //                   format::LinePosition::Bottom],
    //                 format::LineSeparator::new('-', '+', '+', '+'))
    //     .padding(1, 1)
    //     .build();
    // table.set_format(format);
    table.set_format(*format::consts::FORMAT_NO_LINESEP_WITH_TITLE);
    table.add_row(row![FY => "Barcode", "Total Reads", "Total Reads (passed)", "Total Bases (Gb)", "Passed Bases (Gb)", &format!("Passed Bases > {} bp (Gb)", length_threshold)]);

    for (barcode, (reads, total_reads, total_bases, passed_bases, _, passed_bases_gb)) in sorted_barcodes {
        let total_gb = total_bases as f64 / 1_000_000_000.0;
        let passed_gb = passed_bases as f64 / 1_000_000_000.0;
        table.add_row(row![
            barcode,
            total_reads.to_formatted_string(&Locale::en),
            reads.to_formatted_string(&Locale::en),
            format!("{:.2}", total_gb),
            format!("{:.2}", passed_gb),
            format!("{:.2}", passed_bases_gb)
        ]);
    }

    if !output_csv {
        table.printstd();
    }

    // TO-DO: implement an option for user to provide output directory
    // Output to CSV if the --csv option is present
    if output_csv {
        let input_path = Path::new(filename);
        let parent_dir = input_path.parent().unwrap_or_else(|| Path::new("")); // Get the parent directory of the input file

        let mut csv_filename = input_path.file_name().unwrap().to_string_lossy().to_string();
        if csv_filename.ends_with(".txt") {
            csv_filename = csv_filename.trim_end_matches(".txt").to_string();
        }
        csv_filename = format!("{}_barcode_summary.csv", csv_filename);

        let csv_file_path = PathBuf::from(parent_dir).join(csv_filename); // Create the full path for the CSV file

        println!("CSV file generated: {}", csv_file_path.display());

        let out = File::create(&csv_file_path)?; // Create the CSV file in the same directory as the input file
        table.to_csv(out)?;
    }

    Ok(())
}

// main
fn main() -> io::Result<()> {
    let matches = App::new("summary_metrics")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about(env!("CARGO_PKG_DESCRIPTION"))
        .arg(
            Arg::with_name("input_file")
                .help("Nanopore sequencing summary text file")
                .required(true)
                .index(1)
        )
        .arg(
            Arg::with_name("length")
                .short('l')
                .long("length")
                .takes_value(true)
                .help("Minimum length threshold in basepairs")
                .default_value("15000")
        )
        .arg(
            Arg::with_name("qscore")
                .short('q')
                .long("qscore")
                .takes_value(true)
                .help("Quality score threshold")
                .default_value("9.0")
        )
        .arg(
            Arg::with_name("function")
                .help(&*["- ", &"barcode".cyan().to_string(), " \
                 will identify all present barcodes and provide summary metrics for them.\n"].concat())
                .required(false)
                .index(2)
                .possible_values(&["barcode"])
        )
        .arg(
            Arg::with_name("csv")
                .short('o')
                .long("csv")
                .help("Output the barcode table to CSV"),
        )
        .get_matches();

    let filename = matches.value_of("input_file").unwrap();

    // Parse other optional arguments if provided
    let min_length = matches.value_of("length").unwrap_or("15000").parse().unwrap_or_else(|_| {
        eprintln!("Minimum length must be a positive integer.");
        std::process::exit(1);
    });

    let qscore: f64 = matches.value_of("qscore").unwrap_or("9.0").parse().unwrap_or_else(|_| {
        eprintln!("Qscore must be a valid floating-point number.");
        std::process::exit(1);
    });

    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);

    let mut total_line_count = 0;
    let mut passing_read_counts = 0;
    let mut mean_qscore_index = None;
    let mut barcode_index = None;
    let mut sequence_length_index = None;
    let mut total_barcode_counts: HashMap<String, u32> = HashMap::new();
    let mut passing_barcode_counts: HashMap<String, u32> = HashMap::new();
    let mut total_sequence_data_all: u64 = 0;
    let mut total_sequence_data_passed: u64 = 0;
    let mut total_sequence_data_most_common_barcode: u64 = 0;
    let mut long_sequences_most_common_barcode: u64 = 0;
    let mut sequence_lengths_all: Vec<u64> = Vec::new();
    let mut sequence_lengths_passed: Vec<u64> = Vec::new();
    let column_names = vec!["filename_pod5", "run_id", "experiment_id", "sample_id"];

    // Parse the --csv option
    let output_csv = matches.is_present("csv");

    if let Some(function) = matches.value_of("function") {
        match function {
            "barcode" => {
                barcode(filename, qscore, min_length, output_csv)?;
                return Ok(());
            }
            _ => {
                eprintln!("Invalid function specified.");
                return Ok(());
            }
        }
    }

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
    
        if index == 0 {
            mean_qscore_index = fields.iter().position(|&x| x == "mean_qscore_template");
            sequence_length_index = fields.iter().position(|&x| x == "sequence_length_template");
    
            if mean_qscore_index.is_none() || sequence_length_index.is_none() {
                eprintln!("Required information not found in the sequencing summary file header. File must contain columns: sequence_length_template, mean_qscore_template");
                std::process::exit(1);
            }
    
            // Check for barcode_index only if it hasn't been set yet
            if barcode_index.is_none() {
                barcode_index = fields.iter().position(|&x| x == "barcode_arrangement");
            }
    
            continue;
        }
    
        total_line_count += 1;
    
        if let Some(barcode_idx) = barcode_index {
            if let Some(barcode_arrangement) = fields.get(barcode_idx) {
                *total_barcode_counts.entry(barcode_arrangement.to_string()).or_insert(0) += 1;
            }
        }
    
        if let Some(seq_len_idx) = sequence_length_index {
            if let Some(seq_len_str) = fields.get(seq_len_idx) {
                if let Ok(seq_len) = seq_len_str.parse::<u64>() {
                    total_sequence_data_all += seq_len;
                    sequence_lengths_all.push(seq_len);
    
                    if let Some(mean_qscore_idx) = mean_qscore_index {
                        if let Some(passes_filtering) = fields.get(mean_qscore_idx) {
                            if passes_filtering.parse::<f64>().unwrap() >= qscore {
                                passing_read_counts += 1;
                                total_sequence_data_passed += seq_len;
                                sequence_lengths_passed.push(seq_len);
    
                                if let Some(barcode_idx) = barcode_index {
                                    if let Some(barcode_arrangement) = fields.get(barcode_idx) {
                                        *passing_barcode_counts.entry(barcode_arrangement.to_string()).or_insert(0) += 1;
    
                                        if barcode_arrangement == passing_barcode_counts.iter().max_by_key(|&(_, count)| count).map(|(b, _)| b).unwrap() {
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

    println!();
    println!("----------------------- Summary Metrics -----------------------");
    if let Some(ids) = extract_ids(filename, column_names) {
        if let Some(flowcell_id) = &ids[0] {
            if let Some(separator_index) = flowcell_id.find('_') {
                println!("Flowcell ID: {}", &flowcell_id[..separator_index].cyan().bold());
            } else {
                println!("Flowcell ID: {}", flowcell_id.cyan().bold());
            }
        } else {
            eprintln!("Warning: Flowcell ID not found");
        }
        // Print other IDs without processing them
        if let Some(run_id) = &ids[1] {
            println!("Run ID: {}", run_id);
        } else {
            eprintln!("Warning: Run ID not found");
        }

        if let Some(experiment_id) = &ids[2] {
            println!("Experiment ID: {}", experiment_id);
        } else {
            eprintln!("Warning: Experiment ID not found");
        }

        if let Some(sample_id) = &ids[3] {
            println!("Sample ID: {}", sample_id);
        } else {
            eprintln!("Warning: Sample ID not found");
        }
    } else {
        eprintln!("Warning: No identifying information located in sequencing summary file.");
    }
    println!();
    println!("Total reads: {}", total_line_count.to_formatted_string(&Locale::en));
    println!("Total passed reads: {}", passing_read_counts.to_formatted_string(&Locale::en));
    println!();

    // Print barcode information only if barcode_arrangement column is present
    if barcode_index.is_some() {
        if let Some((max_barcode, max_count)) = total_barcode_counts.iter().max_by_key(|&(_, count)| count) {
            println!("Detected barcode (total): {} (count: {})", max_barcode.green().bold(), max_count.to_formatted_string(&Locale::en));
        } else {
            println!("No barcode data found in the file.");
        }

        if let Some((max_barcode, max_count)) = passing_barcode_counts.iter().max_by_key(|&(_, count)| count) {
            println!("Detected barcode (passed): {} (count: {})", max_barcode.green().bold(), max_count.to_formatted_string(&Locale::en));
        } else {
            println!("No barcode data found in the passed reads.");
        }
        println!();
    }

    let total_sequence_data_all_gb = total_sequence_data_all as f64 / 1_000_000_000.0;
    let total_sequence_data_passed_gb = total_sequence_data_passed as f64 / 1_000_000_000.0;
    let total_sequence_data_most_common_barcode_gb = total_sequence_data_most_common_barcode as f64 / 1_000_000_000.0;
    let long_sequences_most_common_barcode_gb = long_sequences_most_common_barcode as f64 / 1_000_000_000.0;

    println!("{}", format!("Total output: {:.2} Gb", total_sequence_data_all_gb).bold());
    println!("Total output (passed): {:.2} Gb", total_sequence_data_passed_gb);

    if barcode_index.is_some() {
        println!("Total output (passed, barcode): {:.2} Gb", total_sequence_data_most_common_barcode_gb);
        println!("Total >= {} bp (passed, barcode): {:.2} Gb", min_length, long_sequences_most_common_barcode_gb);
    }
    if let Some(n50) = calculate_n50(&sequence_lengths_all) {
        let n50_kb = n50 as f64 / 1000.0;
        println!("{}", format!("N50 (total): {:.2} Kb", n50_kb).yellow().bold());
    } else {
        println!("No sequence length data found.");
    }
    println!();

    let mean_read_length_all = calculate_mean(&sequence_lengths_all);
    let median_read_length_all = calculate_median(&mut sequence_lengths_all);
    let mean_read_length_passed = calculate_mean(&sequence_lengths_passed);
    let median_read_length_passed = calculate_median(&mut sequence_lengths_passed);

    println!("Mean read length (before filtering): {:.2} bp", mean_read_length_all);
    println!("Median read length (before filtering): {:.2} bp", median_read_length_all);
    println!("Mean read length (after filtering): {:.2} bp", mean_read_length_passed);
    println!("Median read length (after filtering): {:.2} bp", median_read_length_passed);
    println!("---------------------------- Done -----------------------------");

    Ok(())
}
