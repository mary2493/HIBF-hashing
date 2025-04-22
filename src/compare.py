# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

import subprocess
import time
import csv
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent  # Das Verzeichnis, in dem das Skript liegt
DATA_DIR = BASE_DIR / "test" / "data"  # Der "data"-Ordner im Projektverzeichnis
HIBF_BINARY = BASE_DIR / "build" / "HIBF-hashing"  # Das HIBF Binary im "build"-Ordner
COMPARE_DIR = DATA_DIR / "compare"  # Der "compare"-Unterordner
CSV_FILE = COMPARE_DIR / "hibf_benchmark_results.csv" 

# === Output CSV file ===
CSV_FILE = COMPARE_DIR / "hibf_benchmark_results.csv"

# === Experiments: (hash_type, file_list_path, reads_file_path, arguments) ===
experiments = [
    ("kmer",      DATA_DIR / "provided_list_files.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "20"]),
    ("kmer",      DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "20"]),
    ("minimiser", DATA_DIR / "provided_list_files.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "18", "--window", "20"]),
    ("minimiser", DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "18", "--window", "20"]),
    ("syncmer",   DATA_DIR / "provided_list_files.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
    ("syncmer",   DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
]

# === Collect all benchmarking results ===
results = []

# === Helper function ===
def run_and_measure(cmd, description, tag, step):
    cmd_str = [str(part) for part in cmd]
    print(f"\nRunning: {description}")
    print(f" Command: {' '.join(cmd_str)}")

    start = time.time()
    result = subprocess.run(cmd_str, capture_output=True, text=True)
    end = time.time()
    elapsed = round(end - start, 2)

    if result.returncode != 0:
        output = f"ERROR: {result.stderr.strip()}"
    else:
        output = result.stdout.strip()

    results.append([tag, step, elapsed, output])

    print(f" Time: {elapsed:.2f} seconds")
    print(" Output:", output)  


# === Main Execution ===
if __name__ == "__main__":
    print("Benchmarking HIBF-hashing: Comparing kmer, minimiser, and syncmer")

    for hash_type, file_list, reads_file, build_args in experiments:
        tag = f"{hash_type}_{file_list.stem}"
        index_file = COMPARE_DIR / f"test_index_{tag}.bin"  
        search_output = COMPARE_DIR / f"search_{tag}.txt"  

        # Build
        build_cmd = [HIBF_BINARY, "build", "-i", file_list, "-o", index_file, "-m", hash_type] + build_args
        run_and_measure(build_cmd, f"Build for {hash_type} using {file_list.name}", tag, "build")

        # Search
        search_cmd = [HIBF_BINARY, "search", "-i", index_file, "-r", reads_file, "-o", search_output]
        run_and_measure(search_cmd, f"Search for {hash_type} using {file_list.name}", tag, "search")

    # Save all results to CSV
    with open(CSV_FILE, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Tag", "Step", "Time (s)", "Output (first 100 chars)"])
        writer.writerows(results)

    print(f"\n All results saved to: {CSV_FILE}")
