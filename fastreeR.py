#!/usr/bin/env python3

#
# fastreeR https://github.com/gkanogiannis/fastreeR
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

import threading
import argparse
import subprocess
import zipfile
import sys
import os

# Determine JAR directory
JAR_DIR = os.environ.get("FASTREER_JAR_DIR") or os.path.join(os.path.dirname(__file__), "inst/java")
MEM_MB = "256"
JAVA_PARAMS = ["-Djava.awt.headless", "-XX:+UseG1GC", "-XX:+UseStringDeduplication", "-Xmx"+str(MEM_MB)+"M"]
MAIN_CLASS="com.gkano.bioinfo.javautils.JavaUtils"

def check_java_version(min_major=11):
    try:
        result = subprocess.run(["java", "-version"], stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        version_output = result.stderr.splitlines()[0] if result.stderr else result.stdout.splitlines()[0]
        # Example: 'java version "11.0.20"' or 'openjdk version "17.0.9"'
        if '"' in version_output:
            version_str = version_output.split('"')[1]
            major_version = int(version_str.split('.')[0]) if version_str.startswith("1.") is False else int(version_str.split('.')[1])
            if major_version < min_major:
                print(f"[fastreeR] x Java version {version_str} is too old (need >= {min_major})", file=sys.stderr)
                sys.exit(1)
            else:
                print(f"[fastreeR] ✓ Java version OK: {version_str}", file=sys.stderr)
        else:
            print("[fastreeR] x Unable to parse Java version output.", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"[fastreeR] x Failed to check Java version: {e}", file=sys.stderr)
        sys.exit(1)

def build_classpath(jar_dir):
    if not os.path.isdir(jar_dir):
        print(f"Error: library path '{jar_dir}' does not exist or is not a directory.", file=sys.stderr)
        sys.exit(1)
    jars = [os.path.join(jar_dir, f) for f in os.listdir(jar_dir) if f.endswith(".jar")]
    if not jars:
        print(f"Error: no .jar files found in the library path '{jar_dir}'", file=sys.stderr)
        sys.exit(1)
    separator = ";" if os.name == "nt" else ":"  # Windows uses semicolon, Unix uses colon
    return separator.join(jars)

def run_java_tool(tool_name, params, jar_dir, mem_MB=MEM_MB, output_path=None, verbose=False, extraVerbose=False, pipe_stderr=False, progress_every=100, stdin=None):
    global JAVA_PARAMS
    extra = os.environ.get("FASTREE_JAVA_PARAMS", "")
    if extra:
        JAVA_PARAMS += extra.strip().split()
    JAVA_PARAMS += ["-Xmx"+str(mem_MB)+"M"]
    classpath = build_classpath(jar_dir)
    cmd = ["java"] + JAVA_PARAMS + ["-cp", classpath, MAIN_CLASS, tool_name] + params
    if extraVerbose:
        check_java_version(min_major=11)
        print(f"[fastreeR] JAVA_PARAMS: {' '.join(JAVA_PARAMS)}", file=sys.stderr)
        print(f"[fastreeR] Using JAR directory: {jar_dir}", file=sys.stderr)
        print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    try:
        process = subprocess.Popen(
            cmd,
            stdin=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE if pipe_stderr else None,
            text=True,
            bufsize=1
        )
        
        # Thread to print stderr live (Java progress)
        stderr_thread = None

        if pipe_stderr:
            def stream_stderr(stderr_pipe):
                for line in stderr_pipe:
                    sys.stderr.write(line)
             # Start stderr thread
            if verbose:
                stderr_thread = threading.Thread(target=stream_stderr, args=(process.stderr,))
            else:
                # Still need to consume stderr to avoid blocking, but silently
                stderr_thread = threading.Thread(target=lambda p: [None for _ in p], args=(process.stderr,))
            stderr_thread.start()

        # Stream stdout to file or terminal
        # Write stdout to file line-by-line and print progress every N lines
        line_count = 0
        out_stream = open(output_path, "w") if output_path else sys.stdout
        try:
            for line in process.stdout:
                out_stream.write(line)
                line_count += 1
                if verbose and progress_every and line_count % progress_every == 0:
                    print(f"  ... wrote {line_count} lines", file=sys.stderr, flush=True)
        finally:
            if output_path:
                out_stream.close()
            
        process.wait()
        if stderr_thread:
            stderr_thread.join()
            
        if process.returncode != 0:
            print(f"Java process exited with error code {process.returncode}", file=sys.stderr)
            sys.exit(process.returncode)
        
        if output_path:
            print(f"Wrote {line_count} lines to {output_path}", file=sys.stderr)
        else:
            print(f"Wrote {line_count} lines to stdout", file=sys.stderr)
        
    except subprocess.CalledProcessError as e:
        print(f"Java error: {e}", file=sys.stderr)
        sys.exit(e.returncode)

def main():
    parser = argparse.ArgumentParser(
        description="fastreeR CLI: Calculate distance matrices and phylogenetic trees from VCF or FASTA files\n\n"
                    "Citation:\n"
                    "Anestis Gkanogiannis (2016)\n."
                    "A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes.\n"
                    "BMC Bioinformatics 17, 311 (2016)\n"
                    "https://doi.org/10.1186/s12859-016-1186-3\n"
                    "https://github.com/gkanogiannis/fastreeR\n"
    )
    parser.add_argument("--lib", type=str, default=JAR_DIR, help=f"Path to JAR library folder (default: {JAR_DIR})")
    parser.add_argument("--mem", type=int, default=MEM_MB, help=f"Max RAM for JVM in MB (default: {MEM_MB})")
    parser.add_argument("--pipe-stderr", action="store_true", help="Pipe Java stderr to CLI (default: direct passthrough to terminal)")
    parser.add_argument("--version", action="store_true", help="Print version information and exit")
    parser.add_argument("--check", action="store_true", help="Test Java and backend availability")
    parser.add_argument("--extraVerbose", action="store_true", help="Print extra messages on stderr (default: false)")

    subparsers = parser.add_subparsers(dest="command", required=False)

    def add_common_input_output(parser_obj, allow_multiple_inputs=True):
        parser_obj.add_argument("inputs", nargs="*", help="Positional input files")
        parser_obj.add_argument("-i", "--input", dest="named_inputs", action="append", help="Input file(s)")
        parser_obj.add_argument("-o", "--output", help="Output file path (default: stdout)")

    def add_common_vcf_args(p):
        p.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1)")
        p.add_argument("--ignoreHets", action="store_true", help="Ignore heterozygous loci (default: false)")
        p.add_argument("--onlyHets", action="store_true", help="Use only heterozygous loci (default: false)")
        p.add_argument("--ignoreMissing", action="store_true", help="Ignore missing loci (default: false)")
        p.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")

    # Subcommand for VCF-based distance matrix
    parser_vcf2dist = subparsers.add_parser("VCF2DIST", help="Compute distance matrix from VCF(s)")
    add_common_input_output(parser_vcf2dist)
    add_common_vcf_args(parser_vcf2dist)

    # Subcommand for VCF-based tree
    parser_vcf2tree = subparsers.add_parser("VCF2TREE", help="Compute tree from VCF(s)")
    add_common_input_output(parser_vcf2tree)
    add_common_vcf_args(parser_vcf2tree)
    
    # Subcommand for distance matrix to newick tree
    parser_dist2tree = subparsers.add_parser("DIST2TREE", help="Compute tree from distance matrix")
    parser_dist2tree.add_argument("input_file", nargs="?", help="Input dist file")
    parser_dist2tree.add_argument("-i", "--input", dest="named_input", help="Optional input dist file (overrides positional)")
    parser_dist2tree.add_argument("-o", "--output", help="Output file path (default: stdout)")
    parser_dist2tree.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")
   
   # Subcommand for FASTA-based distance matrix
    parser_fasta2dist = subparsers.add_parser("FASTA2DIST", help="Compute distance matrix from FASTA(s)")
    add_common_input_output(parser_fasta2dist)
    parser_fasta2dist.add_argument("-k", "--kmerSize", type=int, default=4, help="Kmer size for D2S calculation (default: 4)")
    parser_fasta2dist.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser_fasta2dist.add_argument("-n", "--normalize", action="store_true", help="Use normalization (default: false)")
    parser_fasta2dist.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")
   
    args = parser.parse_args()

    if args.check:
        global JAVA_PARAMS
        extra = os.environ.get("FASTREE_JAVA_PARAMS", "")
        if extra:
            JAVA_PARAMS += extra.strip().split()
        JAVA_PARAMS += ["-Xmx"+str(args.mem)+"M"]
        jar_dir = os.environ.get("FASTREER_JAR_DIR") or args.lib
        classpath = build_classpath(jar_dir)
        try:
            cmd = ["java"] + JAVA_PARAMS + ["-version"]
            print(f"[fastreeR] Running Java check: {' '.join(cmd)}", file=sys.stderr)
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print(result.stdout.strip())
            print(result.stderr.strip(), file=sys.stderr)
            cmd = ["java"] + JAVA_PARAMS + ["-cp", classpath, MAIN_CLASS , "VCF2TREE"]
            print(f"[fastreeR] Running Java check: {' '.join(cmd)}", file=sys.stderr)
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print(result.stdout.strip())
            print(result.stderr.strip(), file=sys.stderr)
            print("[fastreeR] v Java check succeeded", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"[fastreeR] x Java check failed (exit code {e.returncode})", file=sys.stderr)
            print(e.stderr, file=sys.stderr)
            sys.exit(e.returncode)
        sys.exit(0)

    if args.version:
        print_version_from_jar(args.lib)
        return

    def resolve_inputs(args, allow_multiple=True):
        combined = []
        if args.named_inputs:
            combined.extend(args.named_inputs)
        if hasattr(args, "inputs"):
            combined.extend(args.inputs)
        if not combined:
            print("Error: No input files provided.", file=sys.stderr)
            sys.exit(1)
        if not allow_multiple and len(combined) > 1:
            print("Error: Only one input file allowed for this command.", file=sys.stderr)
            sys.exit(1)
        return combined

    if args.command in ("VCF2DIST", "VCF2TREE"):
        input_files = resolve_inputs(args)
        params = []
        if args.verbose or args.extraVerbose: params.append("--verbose")
        if args.ignoreHets: params.append("--ignoreHets")
        if args.onlyHets: params.append("--onlyHets")
        if args.ignoreMissing: params.append("--ignoreMissing")
        params.extend(["-t", str(args.threads)])
        for f in input_files:
            params.extend(["-i", f])
        run_java_tool(args.command, params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)

    elif args.command == "DIST2TREE":
        input_file = args.named_input or args.input_file
        if not input_file:
            print("Error: No input distance matrix provided.", file=sys.stderr)
            sys.exit(1)
        params = []
        if args.verbose or args.extraVerbose: params.append("--verbose")
        params.append(input_file)
        run_java_tool("DIST2TREE", params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)
    
    elif args.command == "FASTA2DIST":
        input_files = resolve_inputs(args)
        params = []
        if args.verbose or args.extraVerbose: params.append("--verbose")
        if args.normalize: params.append("--normalize")
        params.extend(["-k", str(args.kmerSize), "-t", str(args.threads)])
        for f in input_files:
            params.extend(["-i", f])
        run_java_tool("FASTA2DIST", params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)
            
    else:
        print("Unknown command", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

def print_version_from_jar(jar_dir):
    jars = [os.path.join(jar_dir, f) for f in os.listdir(jar_dir) if f.endswith(".jar")]
    for jar in jars:
        with zipfile.ZipFile(jar, 'r') as zipf:
            for name in zipf.namelist():
                if "pom.properties" in name:
                    with zipf.open(name) as props:
                        for line in props:
                            decoded = line.decode().strip()
                            if decoded.startswith("version="):
                                version = decoded.split("=")[1]
                                print(f"fastreeR version: {version}")
                                print_citation()
                                return
    print("Version info not found in any jar.", file=sys.stderr)
    print_citation()

def print_citation():
    print("\nCitation:")
    print("Anestis Gkanogiannis (2016).")
    print("A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes.")
    print("BMC Bioinformatics 17, 311 (2016)")
    print("https://doi.org/10.1186/s12859-016-1186-3")
    print("https://github.com/gkanogiannis/fastreeR\n\n")

if __name__ == "__main__":
    main()
