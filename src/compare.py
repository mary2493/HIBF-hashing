# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

import subprocess
import time
from pathlib import Path

# Path to HIBF binary
HIBF_BINARY = Path("/home/mary/develop/HIBF-hashing/build/HIBF-hashing")

# Test cases:
experiments = [
    ("kmer",      Path("/home/mary/develop/HIBF-hashing/test/data/provided_list_files.txt"), Path("/home/mary/develop/HIBF-hashing/test/data/provided/query.fq"),     ["--kmer", "20"]),
    ("kmer",      Path("/home/mary/develop/HIBF-hashing/test/data/file_list.txt"),           Path("/home/mary/develop/HIBF-hashing/test/data/reads.fasta"),  ["--kmer", "20"]),
    ("minimiser", Path("/home/mary/develop/HIBF-hashing/test/data/provided_list_files.txt"), Path("/home/mary/develop/HIBF-hashing/test/data/provided/query.fq"),     ["--kmer", "18", "--window", "20"]),
    ("minimiser", Path("/home/mary/develop/HIBF-hashing/test/data/file_list.txt"),           Path("/home/mary/develop/HIBF-hashing/test/data/reads.fasta"),  ["--kmer", "18", "--window", "20"]),
    ("syncmer",   Path("/home/mary/develop/HIBF-hashing/test/data/provided_list_files.txt"), Path("/home/mary/develop/HIBF-hashing/test/data/provided/query.fq"),     ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
    ("syncmer",   Path("/home/mary/develop/HIBF-hashing/test/data/file_list.txt"),           Path("/home/mary/develop/HIBF-hashing/test/data/reads.fasta"),  ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
]

# Function to run a command and measure its execution time and print the output
def run_and_measure(cmd, description):
    cmd_str = [str(part) for part in cmd]  # Alle Teile in Strings umwandeln
    print(f"\nRunning: {description}")
    print(f" Command: {' '.join(cmd_str)}")
    start = time.time()
    result = subprocess.run(cmd_str, capture_output=True, text=True)
    end = time.time()
    print(f" Time: {end - start:.2f} seconds")
    if result.returncode != 0:
        print(" Error:", result.stderr)
    else:
        print(" Output:", result.stdout)


if __name__ == "__main__":
    print("Benchmarking HIBF-hashing: Comparing kmer, minimiser, and syncmer\n")

    for hash_type, file_list, reads_file, build_args in experiments:
        tag = f"{hash_type}_{Path(file_list).stem}"
        index_file = f"test_index_{tag}.bin"
        search_output = f"search_{tag}.txt"

        # Build command
        build_cmd = [str(HIBF_BINARY), "build", "-i", file_list, "-o", index_file, "-m", hash_type] + build_args
        run_and_measure(build_cmd, f"Build for {hash_type} using {file_list}")

        # Search command
        search_cmd = [str(HIBF_BINARY), "search", "-i", index_file, "-r", reads_file, "-o", search_output]
        run_and_measure(search_cmd, f"Search for {hash_type} using {file_list}")

