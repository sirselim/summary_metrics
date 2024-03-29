# summary_metrics

This is a simple tool designed in Rust to analyze Oxford Nanopore Technologies sequencing summary text files and calculate various statistics.

## Features

- Extracts run information (if present)
  - flowcell ID
  - run ID
  - sample ID
  - experiment ID
- Calculate total output (Gb) of reads
- Determine the N50 value
- Identify the most prevalent barcode
- Count reads that pass filtering criteria, including
  - pass/fail (can be user provided q-score, default is 9.0)
  - barcode (detected if present)
  - given length (provide statistics for reads meeting this criteria)
- calculate basic read length statistics (mean, median)
- the barcode function can be provided to generate summary statistics for all present barcodes
  - this information can be output to a csv file using the `--csv` option with `barcode`

## Prerequisites

You are welcome to download and try the pre-compiled binaries (info below). If you want to build from source, ensure you have the following installed on your system:

- Rust compiler: [Installation instructions](https://www.rust-lang.org/tools/install)

## Installation

### pre-compiled binaries

- [summary_metrics-0.1.7-linux-x64](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.7-linux-x64.tar.gz)
- [summary_metrics-0.1.7-osx-arm](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.7-osx-arm64.tar.gz)

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
./target/release/rust_tool <input_file>
```

Replace <input_file> with the path to your input file. This will run the tool with a default `--length` of `15000`bp and deafult `--qscore` of `9.0`. If you want to change these options you can.

For example:

```bash
./target/release/rust_tool input.txt --length 20000 --qscore 10.0
```

This command will analyze the `input.txt` file using a read length of `20000` bp and qscore threshold of `10.0`.

You should see output similar to below:

```bash
> ./target/release/summary_metrics ../summary_simulator/sequencing_summary_sim_data.txt --length 20000 --qscore 10.0

----------------------- Summary Metrics -----------------------
Flowcell ID: PAU02665
Run ID: f28b19102a1ace40d9afbd35634193fb1a25b699
Experiment ID: 2023_Apr21_plate2
Sample ID: BC01_F1

Total reads: 13,019,341
Total passed reads: 10,727,726

Detected barcode (total): barcode01 (count: 11,301,554)
Detected barcode (passed): barcode01 (count: 10,238,841)

Total output: 123.42 Gb
Total output (passed): 103.84 Gb
Total output (passed, barcode): 102.00 Gb
Total >= 15000 bp (passed, barcode): 59.31 Gb
N50 (total): 17.29 Kb

Mean read length (before filtering): 9479.63 bp
Median read length (before filtering): 5871.00 bp
Mean read length (after filtering): 9679.74 bp
Median read length (after filtering): 6038.00 bp
---------------------------- Done -----------------------------
```

### Processing multiple summary files

If you want to process a collection of summary text files something like below can be useful:

```bash
find ./target/dir -type f -name "sequencing_summary_*.txt" -print0 | \
  xargs -0 -I{} sh -c 'echo "Processing {}"; ./target/release/summary_metrics {} --length 15000 --qscore 9.0'
```

You can also use `gnu parallel` and provide a number of jobs/threads to process at once:

```bash
find ./target/dir -type f -name "sequencing_summary_*.txt" | \
  parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000 --qscore 9.0'
```

#### Outputting to table

We can take the approach above and pass the output to the python script `table_generator.py`. This has an argument, `--format`, which takes `md`, `csv`
and `json` as options. Depending on the format selected the output will be converted to a table and passed to `stdout`, if you want it in a file you
can redirect the output (i.e. `> my_output.csv`). See below for specific examples.

##### Markdown

```bash
# output to markdown
find ./target/dir -type f -name "sequencing_summary_*.txt" | \
  parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000' | \
  python3 ./table_generator.py --format md > my_output.md
```

##### CSV

```bash
# output to csv
find ./target/dir -type f -name "sequencing_summary_*.txt" | \
  parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000' | \
  python3 ./table_generator.py --format csv > my_output.csv
```

##### json

```bash
# output to json
find ./target/dir -type f -name "sequencing_summary_*.txt" | \
  parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000' | \
  python3 ./table_generator.py --format json > my_output.json
```

### Detecting barcodes and get summary metrics

If you have data that has been barcoded, or you want to see if there are barcodes present, then you can run `summary_metrics` with the `barcode` option. If present, barcodes will be sorted on their barcode ID and you can examine the amount of reads and sequence data generated for each barcode.

You can generate general summary statistics for all present barcodes with `summary_metrics barcode`, i.e.

```bash
 ./target/release/summary_metrics ./test/test_summary.txt barcode
+--------------+-------------+----------------------+------------------+-------------------+------------------------------+
| Barcode      | Total Reads | Total Reads (passed) | Total Bases (Gb) | Passed Bases (Gb) | Passed Bases > 15000 bp (Gb) |
| barcode01    | 11,301,554  | 10,238,841           | 112.06           | 102.00            | 59.32                        |
| barcode02    | 484         | 69                   | 0.00             | 0.00              | 0.00                         |
| barcode03    | 259         | 43                   | 0.00             | 0.00              | 0.00                         |
| barcode04    | 146         | 52                   | 0.00             | 0.00              | 0.00                         |
| barcode05    | 1,463       | 134                  | 0.01             | 0.00              | 0.00                         |
| barcode06    | 129         | 69                   | 0.00             | 0.00              | 0.00                         |
| barcode07    | 54          | 13                   | 0.00             | 0.00              | 0.00                         |
| barcode08    | 38          | 18                   | 0.00             | 0.00              | 0.00                         |
| barcode09    | 149         | 52                   | 0.00             | 0.00              | 0.00                         |
| barcode10    | 162         | 42                   | 0.00             | 0.00              | 0.00                         |
| barcode11    | 64          | 28                   | 0.00             | 0.00              | 0.00                         |
| barcode12    | 276         | 106                  | 0.00             | 0.00              | 0.00                         |
| barcode13    | 213         | 51                   | 0.00             | 0.00              | 0.00                         |
| barcode14    | 59          | 17                   | 0.00             | 0.00              | 0.00                         |
| barcode15    | 127         | 17                   | 0.00             | 0.00              | 0.00                         |
| barcode16    | 147         | 33                   | 0.00             | 0.00              | 0.00                         |
| barcode17    | 62          | 20                   | 0.00             | 0.00              | 0.00                         |
| barcode18    | 3,704       | 365                  | 0.03             | 0.00              | 0.00                         |
| barcode19    | 172         | 55                   | 0.00             | 0.00              | 0.00                         |
| barcode20    | 116         | 30                   | 0.00             | 0.00              | 0.00                         |
| barcode21    | 146         | 63                   | 0.00             | 0.00              | 0.00                         |
| barcode22    | 109         | 33                   | 0.00             | 0.00              | 0.00                         |
| barcode23    | 255         | 41                   | 0.00             | 0.00              | 0.00                         |
| barcode24    | 187         | 62                   | 0.00             | 0.00              | 0.00                         |
| unclassified | 1,647,593   | 483,017              | 11.07            | 1.81              | 0.73                         |
+--------------+-------------+----------------------+------------------+-------------------+------------------------------+
```

It's simple to 'clean up' the output table. We can use something like `grep`, i.e.

```bash
./target/release/summary_metrics ./test/test_summary.txt barcode | grep -v '0.00'
+--------------+-------------+----------------------+------------------+-------------------+------------------------------+
| Barcode      | Total Reads | Total Reads (passed) | Total Bases (Gb) | Passed Bases (Gb) | Passed Bases > 15000 bp (Gb) |
| barcode01    | 11,301,554  | 10,238,841           | 112.06           | 102.00            | 59.32                        |
| unclassified | 1,647,593   | 483,017              | 11.07            | 1.81              | 0.73                         |
+--------------+-------------+----------------------+------------------+-------------------+------------------------------+
```

#### output to csv

You can output the barcode summary table to a csv file, i.e.

```bash
./target/release/summary_metrics ./test/test_summary.txt barcode --csv
```

This will generate a csv file that has the original file name with `_barcode_summary.csv` appended. It will be in the same directory as the input data.

## To Do

- [X] ~~remove hard coded columns, use header values~~
- [X] ~~add proper help options~~
- [X] ~~basic read statistics~~
- [X] ~~extract info for flowcell, sample, experiment etc.~~
- [X] ~~additional filtering options (user defined qscore)~~
- [X] ~~work on error handling~~
- [X] ~~add option to detect and provide info on present barcodes~~
- [ ] add ability to provide barcode(s) to the tool
- [ ] refactor the mess!
- [ ] add --json output option
  - [ ] explore further customizable output formats
- [ ] explore visualisation options
- [ ] further performance optimisations

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributions

Contributions are welcome! Please feel free to submit pull requests or open issues if you encounter any problems or have suggestions for improvements.
