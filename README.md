# summary_metrics

This is a simple tool designed in Rust analyze Oxford Nanopore Technologies sequencing summary text files and calculate various statistics 

## Features

- Extracts run information
  - flowcell ID
  - run ID
  - sample ID
  - experiment ID
- Calculate total output (Gb) of reads
- Determine the N50 value
- Identify the most prevalent barcode
- Count reads that pass filtering criteria, including
  - pass/fail
  - barcode
  - given length
- calculate basic read length statistics (mean, median)

## Prerequisites

You are welcome to download and try the pre-compiled binaries (info below). If you want to build from source, ensure you have the following installed on your system:

- Rust compiler: [Installation instructions](https://www.rust-lang.org/tools/install)

## Installation

### pre-compiled binaries

- [summary_metrics-0.1.3-linux-x64](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.3-linux-x64.tar.gz)
- [summary_metrics-0.1.3-osx-arm](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.3-osx-arm64.tar.gz)

### From source

To build this tool from source, follow these steps:

1. Clone this repository:

```bash
git clone https://github.com/sirselim/summary_metrics.git
```

Navigate to the project directory:

```bash
cd summary_metrics
```

Build the project using Cargo:

```bash
cargo build --release
```

This will compile the Rust code and generate the executable binary in the target/release directory.

## Usage

Once the binary is built, you can run the tool with the following command:

```bash
./target/release/rust_tool <input_file> <read_length>
```

Replace <input_file> with the path to your input file and <read_length> with the desired read length parameter.

For example:

```bash
./target/release/rust_tool input.txt 15000
```

This command will analyze the `input.txt`` file using a read length of 15000 bp.

You should see output similar to below:

```bash
> ./target/release/summary_metrics ../summary_simulator/sequencing_summary_sim_data.txt 15000

----------------------- Summary Metrics -----------------------
Flowcell ID: PAU02321
Run ID: f28b19103a1ace40d6afbd35638193fb1a25b699
Experiment ID: 2024_Mar03_test
Sample ID: dummy_data

Total reads: 1000000
Total passed reads: 888805

Detected barcode (total): barcode02 (count: 850103)
Detected barcode (passed): barcode02 (count: 755502)

Total output: 9.99 Gb
Total output (passed): 8.88 Gb
Total output (passed, barcode): 7.55 Gb
Total >= 15000 bp (passed, barcode): 3.94 Gb
N50 (total): 15.65 Kb

Mean read length (before filtering): 9993.77 bp
Median read length (before filtering): 7390.00 bp
Mean read length (after filtering): 9990.51 bp
Median read length (after filtering): 7387.00 bp
---------------------------- Done -----------------------------
```

### Processing multiple summary files

If you want to process a collection of summary text files something like below can be useful:

```bash
find ./target/dir -type f -name "sequencing_summary_*.txt" -print0 | xargs -0 -I{} sh -c 'echo "Processing {}"; ./target/release/summary_metrics {} 15000'
```

You can also use `gnu parallel` and provide a number of jobs/threads to process at once:

```bash
find ./experiment_dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} 15000'
```

#### Outputting to table

We can take the approach above and pass the output to the python script `table_generator.py`. This has an argument, `--format`, which takes `md`, `csv` 
and `json` as options, depending on the format selected.

##### Markdown

```bash
# output to markdown
find ./experiment_dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} 15000' | python3 ./table_generator.py --format md > my_output.md
```

##### CSV

```bash
# output to csv
find ./experiment_dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} 15000' | python3 ./table_generator.py --format csv > my_output.csv
```

##### json

```bash
# output to json
find ./experiment_dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} 15000' | python3 ./table_generator.py --format json > my_output.json
```

## To Do

- [X] ~~remove hard coded columns, use header values~~
- [X] ~~add proper help options~~
- [X] ~~basic read statistics~~
- [X] ~~extract info for flowcell, sample, experiment etc.~~
- [ ] refactor the mess!
- [ ] add --json output option
  - [ ] explore further customizable output formats
- [ ] additional filtering options
- [ ] explore visualisation options
- [ ] work on error handling
- [ ] further performance optimisations

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributions

Contributions are welcome! Please feel free to submit pull requests or open issues if you encounter any problems or have suggestions for improvements.
