import os
import re
import csv
import gzip

# Directory containing the log files
log_dir = "benchmark_logs"

# Output CSV file
output_csv = "benchmark_logs.csv"

filename_pattern = re.compile(r"M(\d+)_S(\d+)_t(\d+)\.(fastreer|vcf2dis)\.log")
time_pattern = re.compile(r"Time:\s*(?:(\d+):)?(\d+):(\d+(?:\.\d+)?)")
memory_pattern = re.compile(r"Memory:\s*(\d+)\s*KB")

rows = []

def open_log_file(filepath):
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt")
    else:
        return open(filepath, "r")

for file in os.listdir(log_dir):
    match = filename_pattern.match(file)
    if not match:
        continue

    snps_million, samples, threads, tool = match.groups()
    snps = int(snps_million) * 1_000_000
    samples = int(samples)
    threads = int(threads)
    tool = tool.lower()

    filepath = os.path.join(log_dir, file)
    with open_log_file(filepath) as f:
        content = f.read()

    time_match = time_pattern.search(content)
    memory_match = memory_pattern.search(content)

    if not time_match or not memory_match:
        continue

    h = int(time_match.group(1)) if time_match.group(1) else 0
    m = int(time_match.group(2))
    s = float(time_match.group(3))
    total_seconds = h * 3600 + m * 60 + s
    
    memory_kb = int(memory_match.group(1))
    memory_mb = memory_kb / 1024

    rows.append({
        "tool": tool,
        "snps": snps,
        "samples": samples,
        "threads": threads,
        "runtime_seconds": total_seconds,
        "memory_mb": memory_mb,
    })

rows.sort(key=lambda x: (x["samples"], x["snps"], x["threads"], x["tool"]))

with open(output_csv, "w", newline="") as csvfile:
    fieldnames = ["tool", "snps", "samples", "threads", "runtime_seconds", "memory_mb"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"CSV summary written to {output_csv}")
