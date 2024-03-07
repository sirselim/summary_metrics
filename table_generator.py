#!/usr/bin/env python3

# author: Miles Benton
# created: 2024/03/07 11:22:34
# last modified: 2024/03/07 13:30:09
# description: a small script to take the output of summary_metrics (https://github.com/sirselim/summary_metrics)
# and generate a table in a user specified format (md, csv, json). The default behaviour is outputting to md.

# imports
import re
import sys
import argparse
import json

# parsing function
def parse_summary_metrics(output):
    entries = re.split(r'Processing', output)[1:]
    results = []
    for entry in entries:
        summary = {}
        lines = entry.strip().split('\n')
        for line in lines:
            match_flowcell = re.search(r'Flowcell ID: (.+)', line)
            if match_flowcell:
                summary['Flowcell ID'] = match_flowcell.group(1)
            match_run_id = re.search(r'Run ID: (.+)', line)
            if match_run_id:
                summary['Run ID'] = match_run_id.group(1)
            match_exp_id = re.search(r'Experiment ID: (.+)', line)
            if match_exp_id:
                summary['Experiment ID'] = match_exp_id.group(1)
            match_sample_id = re.search(r'Sample ID: (.+)', line)
            if match_sample_id:
                summary['Sample ID'] = match_sample_id.group(1)
            match_total_reads = re.search(r'Total reads: (\d+)', line)
            if match_total_reads:
                summary['Total reads'] = int(match_total_reads.group(1))
            match_passed_reads = re.search(r'Total passed reads: (\d+)', line)
            if match_passed_reads:
                summary['Total passed reads'] = int(match_passed_reads.group(1))
            match_detected_barcode = re.search(r'Detected barcode \(total\): (\w+) \(count: (\d+)\)', line)
            if match_detected_barcode:
                barcode, count = match_detected_barcode.group(1), match_detected_barcode.group(2)
                summary['Detected barcode'] = f'{barcode}'
            match_total_output = re.search(r'Total output: (\d+\.\d+) Gb', line)
            if match_total_output:
                summary['Total output (Gb)'] = float(match_total_output.group(1))
            match_passed_output = re.search(r'Total output \(passed\): (\d+\.\d+) Gb', line)
            if match_passed_output:
                summary['Total output passed (Gb)'] = float(match_passed_output.group(1))
            match_passed_barcode_output = re.search(r'Total output \(passed, barcode\): (\d+\.\d+) Gb', line)
            if match_passed_barcode_output:
                summary['Total output passed barcode (Gb)'] = float(match_passed_barcode_output.group(1))
            match_passed_barcode_15000bp_output = re.search(r'Total >= 15000 bp \(passed, barcode\): (\d+\.\d+) Gb', line)
            if match_passed_barcode_15000bp_output:
                summary['Total >= 15000 bp passed barcode (Gb)'] = float(match_passed_barcode_15000bp_output.group(1))
            match_n50 = re.search(r'N50 \(total\): (\d+\.\d+) Kb', line)
            if match_n50:
                summary['N50 (total)'] = float(match_n50.group(1))
            match_mean_before_filtering = re.search(r'Mean read length \(before filtering\): (\d+\.\d+) bp', line)
            if match_mean_before_filtering:
                summary['Mean read length before filtering  (bp)'] = float(match_mean_before_filtering.group(1))
            match_median_before_filtering = re.search(r'Median read length \(before filtering\): (\d+\.\d+) bp', line)
            if match_median_before_filtering:
                summary['Median read length before filtering (bp)'] = float(match_median_before_filtering.group(1))
            match_mean_after_filtering = re.search(r'Mean read length \(after filtering\): (\d+\.\d+) bp', line)
            if match_mean_after_filtering:
                summary['Mean read length after filtering (bp)'] = float(match_mean_after_filtering.group(1))
            match_median_after_filtering = re.search(r'Median read length \(after filtering\): (\d+\.\d+) bp', line)
            if match_median_after_filtering:
                summary['Median read length after filtering (bp)'] = float(match_median_after_filtering.group(1))
        results.append(summary)
    return results

# table formatting function
def print_table(data, output_format='md'):
    if output_format == 'md':
        headers = data[0].keys()
        print('|', end='')
        for header in headers:
            print(f' {header} |', end='')
        print()
        print('|', end='')
        for _ in headers:
            print(' --- |', end='')
        print()
        for entry in data:
            print('|', end='')
            for header in headers:
                print(f' {entry.get(header, "")} |', end='')
            print()
    elif output_format == 'csv':
        headers = data[0].keys()
        print(','.join(headers))
        for entry in data:
            print(','.join(str(entry.get(header, "")) for header in headers))
    elif output_format == 'json':
        print(json.dumps(data, indent=4))

# Parse command line arguments
parser = argparse.ArgumentParser(description='Parse summary metrics.')
parser.add_argument('--format', choices=['md', 'csv', 'json'], default='md', help='Output format (md, csv, or json)')
args = parser.parse_args()

# Read from stdin
output = sys.stdin.read()
parsed_data = parse_summary_metrics(output)
print_table(parsed_data, args.format)
