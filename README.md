# summary_metrics

This is a simple tool designed in Rust analyze Oxford Nanopore Technologies sequencing summary text files and calculate various statistics 

## Features

- Calculate total output (Gb) of reads
- Determine the N50 value
- Identify the most prevalent barcode
- Count reads that pass filtering criteria, including
  - pass/fail
  - barcode
  - given length

## Prerequisites

You are welcome to download and try the pre-compiled binaries (info below). If you want to build from source, ensure you have the following installed on your system:

- Rust compiler: [Installation instructions](https://www.rust-lang.org/tools/install)

## Installation

### pre-compiled binaries

- [summary_metrics-0.1.0-linux-x64](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.0-linux-x64.tar.gz)
- [summary_metrics-0.1.0-osx-arm](https://github.com/sirselim/summary_metrics/raw/main/binaries/summary_metrics-0.1.0-osx-arm64.tar.gz)

### from source

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

This command will analyze the input.txt file with a read length of 15000 bp.

You should see output similar to below:

```bash
$ ./target/release/sequencing_analysis test_summary.txt 15000
Most prevalent barcode: barcode01 (Count: 10238841)
Number of reads that pass: 10727726
Number of reads that pass with most prevalent barcode: 10238841
Total gigabases of reads that pass: 103.84
Total gigabases of reads that pass with the detected barcode: 102.00
Total gigabases of reads that pass with the detected barcode and are >= 15000bp: 59.32
N50: 17.49 Kb
```

### processing multiple summary files

If you want to process a collection of summary text files something like below can be useful:

```bash
find ./target/dir -type f -name "sequencing_summary_*.txt" -print0 | xargs -0 -I{} sh -c 'echo "Processing {}"; ./target/release/summary_metrics {} 15000'
```

You can also use `gnu parallel` and provide a number of jobs/threads to process at once:

```bash
find ../22_samples -type f -name "sequencing_summary_*.txt" | parallel -j 24 'echo -e "\nProcessing {}"; ./target/release/summary_metrics {} 15000'
```

## To Do

- [X] ~~remove hard coded columns, use header values~~
- [X] ~~add --json output option~~
- [ ] explore further customizable output formats
- [ ] additional filtering options
- [ ] explore visualisation options
- [ ] work on error handling
- [ ] further performance optimisations

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributions

Contributions are welcome! Please feel free to submit pull requests or open issues if you encounter any problems or have suggestions for improvements.
