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
- can detect barcodes and provide information on reads and output per barcode (if barcodes are present)

## Prerequisites

You are welcome to download and try the pre-compiled binaries (info below). If you want to build from source, ensure you have the following installed on your system:

- Rust compiler: [Installation instructions](https://www.rust-lang.org/tools/install)

## Installation

### pre-compiled binaries

- [summary_metrics-0.1.5-linux-x64](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.5-linux-x64.tar.gz)
- [summary_metrics-0.1.5-osx-arm](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.5-osx-arm64.tar.gz)

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
Total >= 20000 bp (passed, barcode): 3.94 Gb
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
find ./target/dir -type f -name "sequencing_summary_*.txt" -print0 | xargs -0 -I{} sh -c 'echo "Processing {}"; ./target/release/summary_metrics {} --length 15000 --qscore 9.0'
```

You can also use `gnu parallel` and provide a number of jobs/threads to process at once:

```bash
find ./target/dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000 --qscore 9.0'
```

#### Outputting to table

We can take the approach above and pass the output to the python script `table_generator.py`. This has an argument, `--format`, which takes `md`, `csv`
and `json` as options. Depending on the format selected the output will be converted to a table and passed to `stdout`, if you want it in a file you
can redirect the output (i.e. `> my_output.csv`). See below for specific examples.

##### Markdown

```bash
# output to markdown
find ./target/dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000' | python3 ./table_generator.py --format md > my_output.md
```

##### CSV

```bash
# output to csv
find ./target/dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000' | python3 ./table_generator.py --format csv > my_output.csv
```

##### json

```bash
# output to json
find ./target/dir -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} --length 15000' | python3 ./table_generator.py --format json > my_output.json
```

### Detecting barcodes

If you have data that has been barcoded, or you want to see if there are barcodes present, then you can run `summary_metrics` with the `detect_barcodes` option, i.e.

```bash
./target/release/summary_metrics ./test/test_summary.txt detect_barcodes
------------ Barcode Detection ------------
Barcode      Reads         Bases           
------------ ------------- ----------------
-                   34,369                0
barcode01       11,301,554  112,064,915,918
barcode02              486        4,143,013
barcode03              275        2,253,332
barcode04              153        1,189,210
barcode05            1,472       13,223,181
barcode06              133          626,559
barcode07               74          541,958
barcode08               51          309,857
barcode09              167        1,132,384
barcode10              168        1,385,959
barcode11               70          323,569
barcode12              286        1,593,870
barcode13              219        1,858,268
barcode14               63          422,996
barcode15              132        1,032,240
barcode16              149        1,031,338
barcode17               64          354,610
barcode18            3,705       33,224,704
barcode19              175        1,035,478
barcode20              117          869,446
barcode21              146          867,648
barcode22              118          667,421
barcode23              256        1,956,280
barcode24              197        1,133,617
unclassified     1,647,593   11,066,611,746
------------------- Done ------------------
```

If present, barcodes will be sorted on their barcode ID and you can examine the amount of reads and sequence data generated for each barcode.

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
