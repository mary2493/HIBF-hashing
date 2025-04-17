import subprocess
import time
import os

# === CONFIGURATION ===
exec = "/home/mary/develop/HIBF-hashing/build/HIBF-hashing"     # ggf. Pfad zu deinem Build-Binary
file_list = "/home/mary/develop/HIBF-hashing/test/data/file_list.txt"
reads = "/home/mary/develop/HIBF-hashing/test/data/reads.fasta"
index_prefix = "hibf_index"
output_dir = "benchmark_results"
os.makedirs(output_dir, exist_ok=True)

methods = ["kmer", "minimiser"]  # syncmer kannst du später einfach ergänzen
results = {}

# === BENCHMARKING LOOP ===
for method in methods:
    build_times = []
    search_times = []
    print(f"\n--- Running benchmarks for {method.upper()} ---")

    for run in range(3):
        print(f"  Run {run+1}/3")

        index_path = f"{output_dir}/{index_prefix}_{method}_{run}.hibf"
        out_path = f"{output_dir}/output_{method}_{run}.txt"

        # BUILD
        build_cmd = [
            exec,
            "build",
            "--input", file_list,
            "--output", index_path,
            "--hash", method,
            "--kmer", "20",        # ggf. anpassen
            "--window", "22"       # für minimiser
        ]

        start = time.time()
        subprocess.run(build_cmd, check=True)
        end = time.time()
        build_times.append(end - start)

        # SEARCH
        search_cmd = [
            exec,
            "search",
            "--reads", reads,
            "--index", index_path,
            "--output", out_path,
            "--hash", method,
            "--kmer", "20",
            "--window", "22"
        ]

        start = time.time()
        subprocess.run(search_cmd, check=True)
        end = time.time()
        search_times.append(end - start)

    # Speichern der Zeiten
    results[method] = {
        "build_avg": sum(build_times) / len(build_times),
        "search_avg": sum(search_times) / len(search_times),
        "build_all": build_times,
        "search_all": search_times
    }

# === AUSGABE DER ERGEBNISSE ===
print("\n=== AVERAGE TIMES ===")
for method in methods:
    print(f"{method.upper():<10} | Build: {results[method]['build_avg']:.3f}s | Search: {results[method]['search_avg']:.3f}s")
