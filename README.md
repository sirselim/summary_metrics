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

Before using this tool, ensure you have the following installed on your system:

- Rust compiler: [Installation instructions](https://www.rust-lang.org/tools/install)

## Installation

To build this tool from source, follow these steps:

1. Clone this repository:

```bash
git clone https://github.com/your_username/your_repository.git
```

Navigate to the project directory:

```bash
cd your_repository
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

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributions

Contributions are welcome! Please feel free to submit pull requests or open issues if you encounter any problems or have suggestions for improvements.
