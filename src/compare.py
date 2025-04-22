# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

import subprocess
import time
import csv
from pathlib import Path

# === Pfade ===
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "test" / "data"
HIBF_BINARY = BASE_DIR / "build" / "HIBF-hashing"
COMPARE_DIR = DATA_DIR / "compare"
CSV_FILE = COMPARE_DIR / "hibf_benchmark_results.csv"

COMPARE_DIR.mkdir(parents=True, exist_ok=True)  # Ordner sicherstellen

# === Experimente ===
experiments = [
    ("kmer",      DATA_DIR / "provided_list_files.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "20"]),
    ("kmer",      DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "20"]),
    ("minimiser", DATA_DIR / "provided_list_files.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "18", "--window", "20"]),
    ("minimiser", DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "18", "--window", "20"]),
    ("syncmer",   DATA_DIR / "provided_list_files.txt", DATA_DIR / "provided/query.fq",     ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
    ("syncmer",   DATA_DIR / "file_list.txt",           DATA_DIR / "reads.fasta",           ["--kmer", "15", "--syncmer_s", "11", "--syncmer_t", "2"]),
]

# === Ergebnisse gruppiert speichern ===
grouped_results = {}

# === Helper Function ===
def run_and_measure(cmd, description, method, step, file_key):
    cmd_str = [str(part) for part in cmd]
    print(f"\nRunning: {description}")
    print(f" Command: {' '.join(cmd_str)}")

    start = time.time()
    result = subprocess.run(cmd_str, capture_output=True, text=True)
    end = time.time()
    elapsed = round(end - start, 2)

    output = result.stdout.strip() if result.returncode == 0 else f"ERROR: {result.stderr.strip()}"

    # Initialisiere Gruppe bei Bedarf
    if file_key not in grouped_results:
        grouped_results[file_key] = {
            "files": file_key,
            "kmer": {}, "minimiser": {}, "syncmer": {}
        }

    # Speichere Ergebnis
    grouped_results[file_key][method][step] = {
        "time": elapsed,
        "output": output
    }

    print(f" Time: {elapsed:.2f} seconds")
    print(" Output:", output)


# === Main Execution ===
if __name__ == "__main__":
    print("Benchmarking HIBF-hashing: Comparing kmer, minimiser, and syncmer")

    for hash_type, file_list, reads_file, build_args in experiments:
        tag = f"{hash_type}_{file_list.stem}"
        file_key = (str(file_list), str(reads_file))
        index_file = COMPARE_DIR / f"test_index_{tag}.bin"
        search_output = COMPARE_DIR / f"search_{tag}.txt"

        # Build
        build_cmd = [HIBF_BINARY, "build", "-i", file_list, "-o", index_file, "-m", hash_type] + build_args
        run_and_measure(build_cmd, f"Build for {hash_type} using {file_list.name}", hash_type, "build", file_key)

        # Search
        search_cmd = [HIBF_BINARY, "search", "-i", index_file, "-r", reads_file, "-o", search_output]
        run_and_measure(search_cmd, f"Search for {hash_type} using {file_list.name}", hash_type, "search", file_key)

    # === Ausgabe als CSV-Datei ===
    with open(CSV_FILE, "w", newline="") as f:
        writer = csv.writer(f, delimiter=";")

        for file_key, data in grouped_results.items():
            seq_file, read_file = file_key
            writer.writerow(["Data name seq:", seq_file])
            writer.writerow(["Data name reads:", read_file])

            for method in ["kmer", "minimiser", "syncmer"]:
                writer.writerow([f"{method}:"])
                build = data[method].get("build", {})
                search = data[method].get("search", {})

                writer.writerow(["Build:", build.get("time", "N/A")])
                writer.writerow(["Search:", search.get("time", "N/A")])
                writer.writerow(["The following hits were found:"])
                writer.writerow([search.get("output", "")])

    print(f"\n Alle Ergebnisse gespeichert unter: {CSV_FILE}")
