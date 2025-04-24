# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

import subprocess
import time
from pathlib import Path

# this script is used to compare the performance of different hash types
# Setup directories and paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "test" / "data"
HIBF_BINARY = BASE_DIR / "build" / "HIBF-hashing"
COMPARE_DIR = BASE_DIR / "test" / "benchmark" / "output"
TXT_FILE = COMPARE_DIR / "hibf_benchmark_results.txt"

# Experiments for building and searching without errors
experiments = [
    ("minimiser",      DATA_DIR / "list.txt", DATA_DIR / "query.fq", ["--kmer", "20", "--window", "24"]),
    # ("kmer",      DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "20"]),
    ("minimiser", DATA_DIR / "list.txt", DATA_DIR / "query.fq",      ["--kmer", "20", "--window", "24"]),
    # ("minimiser", DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "18", "--window", "20"]),
    ("syncmer",   DATA_DIR / "list.txt", DATA_DIR / "query.fq",      ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
    # ("syncmer",   DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
]

def run_and_measure(cmd, description):
    print(f"\nRunning: {description}")
    print(f" Command: {' '.join(str(c) for c in cmd)}")
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.time()

    success = result.returncode == 0
    elapsed = round(end - start, 2)
    output = result.stdout.strip() if success else f"ERROR: {result.stderr.strip()}"
    print("Success" if success else "Failed")
    return elapsed, output

def clean_hits(text):
    return "\n".join(line for line in text.splitlines() if not line.strip().startswith("The following hits were found:"))

if __name__ == "__main__":
    grouped_results = {}

    # Build and normal search (error=0)
    for hash_type, file_list, reads_file, build_args in experiments:
        key = (file_list.name, reads_file.name)
        if key not in grouped_results:
            grouped_results[key] = {}

        tag = f"{hash_type}_{file_list.stem}"
        index_file = COMPARE_DIR / f"test_index_{tag}.bin"
        search_output = COMPARE_DIR / f"search_{tag}.txt"

        # build
        build_cmd = [HIBF_BINARY, "build", hash_type, "-i", file_list, "-o", index_file] + build_args
        build_time, _ = run_and_measure(build_cmd, f"Build for {hash_type} using {file_list.name}")

        # search
        search_cmd = [HIBF_BINARY, "search", "-i", index_file, "-r", reads_file, "-o", search_output]
        search_time, search_output_text = run_and_measure(search_cmd, f"Search for {hash_type} using {file_list.name}")

        grouped_results[key][hash_type] = {
            "build_time": build_time,
            "search_time": search_time,
            "hits": "(error=0)\n" + clean_hits(search_output_text)
        }

    # Only search with --error 2 for kmer/minimiser (no new builds)
    extra_error_experiments = [
        ("kmer",      DATA_DIR / "list.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "20"]),
        ("minimiser", DATA_DIR / "list.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "20", "--window", "24"]),
    ]

    for hash_type, file_list, reads_file, _ in extra_error_experiments:
        key = (file_list.name, reads_file.name)
        tag = f"{hash_type}_{file_list.stem}_error2"
        index_file = COMPARE_DIR / f"test_index_{hash_type}_{file_list.stem}.bin"
        search_output = COMPARE_DIR / f"search_{tag}.txt"

        search_cmd = [HIBF_BINARY, "search", "-i", index_file, "-r", reads_file, "-o", search_output, "--error", "2"]
        search_time, search_output_text = run_and_measure(search_cmd, f"Search for {hash_type} using {file_list.name} with error=2")

        grouped_results[key][f"{hash_type}_error2"] = {
            "build_time": "-",
            "search_time": search_time,
            "hits": "(error=2)\n" + clean_hits(search_output_text)
        }

    # Save results to file
    with open(TXT_FILE, "w") as f:
        for (file_list_name, reads_name), data in grouped_results.items():
            f.write(f"Data name of seq: {file_list_name}\n")
            f.write(f"Data name of reads: {reads_name}\n\n")

            headers = [
                "kmer (err=0)",
                "minimiser (err=0)",
                "syncmer (err=0)",
                "kmer (err=2)",
                "minimiser (err=2)"
            ]

            f.write("".ljust(25))
            f.write("\t".join(h.ljust(25) for h in headers) + "\n")
            f.write("Time\n")

            for kind in ["build_time", "search_time"]:
                label = "B" if kind == "build_time" else "S"
                times = [f"{label}:{data[ht][kind]:.2f}s" if ht in data and data[ht][kind] != "-" else "-" for ht in ["kmer", "minimiser", "syncmer", "kmer_error2", "minimiser_error2"]]
                f.write("".ljust(25))
                f.write("\t".join(t.ljust(25) for t in times) + "\n")

            f.write("\n" + "Hits".ljust(25) + "\n")

            hits_by_line = [[] for _ in range(5)]
            max_lines = 0

            for i, ht in enumerate(["kmer", "minimiser", "syncmer", "kmer_error2", "minimiser_error2"]):
                if ht in data:
                    lines = data[ht]["hits"].splitlines()
                    hits_by_line[i] = lines
                    max_lines = max(max_lines, len(lines))

            for i in range(max_lines):
                row = []
                for lines in hits_by_line:
                    row.append(lines[i] if i < len(lines) else "")
                f.write("".ljust(25))
                f.write("\t".join(cell.ljust(25) for cell in row) + "\n")

            f.write("\n" + "="*120 + "\n\n")
    print(f"\nResults were saved to: {TXT_FILE}")
