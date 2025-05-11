#!/usr/bin/env python3

import threading
import tempfile
import argparse
import subprocess
import zipfile
import sys
import os

JAR_DIR = os.path.join(os.path.dirname(__file__), "../inst/java")
JAVA_PARAMS = ["-Djava.awt.headless=true", "-XX:+UseG1GC", "-XX:+UseStringDeduplication"]
MEM_GB = "1"

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

def run_java_tool(tool_name, params, jar_dir, mem_GB=MEM_GB, output_path=None, verbose=False, pipe_stderr=False, progress_every=100, stdin=None):
    classpath = build_classpath(jar_dir)
    cmd = ["java"]
    cmd.extend(JAVA_PARAMS)
    cmd.extend(["-Xms"+str(mem_GB)+"G", "-Xmx"+str(mem_GB)+"G"])
    cmd.extend(["-cp", classpath, "ciat.agrobio.javautils.JavaUtils", tool_name])
    cmd.extend(params)
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    try:
        
        process = subprocess.Popen(
            cmd,
            stdin=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE if pipe_stderr  else None,
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
                    "📄 Citation:\n"
                    "Anestis Gkanogiannis (2016)\n."
                    "A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes.\n"
                    "BMC Bioinformatics 17, 311 (2016)\n"
                    "https://doi.org/10.1186/s12859-016-1186-3\n"
                    "https://github.com/gkanogiannis/fastreeR\n"
    )
    parser.add_argument(
        "--lib", 
        type=str, 
        default=JAR_DIR,
        help=f"Path to the folder containing JAR libraries (default: {JAR_DIR})"
    )
    parser.add_argument(
        "--mem", 
        type=int, 
        default=MEM_GB,
        help=f"Max RAM for JVM in GB (default: {MEM_GB})"
    )
    parser.add_argument(
        "--pipe-stderr",
        action="store_true",
        help="Pipe stderr and forward from Python (default: direct passthrough to terminal)"
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version information and exit"
    )

    subparsers = parser.add_subparsers(dest="command", required=False)

    # Subcommand for VCF-based distance matrix
    parser_vcf2dist = subparsers.add_parser("VCF2DIST", help="Compute distance matrix from VCF")
    parser_vcf2dist.add_argument("vcf_file", nargs="?", help="Path to input VCF file (positional)")
    parser_vcf2dist.add_argument("-i", "--input", help="Optional input VCF file (overrides positional)")
    #parser_vcf2dist.add_argument("-d", "--distance", choices=["cosine"], default="cosine", help="Distance metric")
    parser_vcf2dist.add_argument("-o", "--output", help="Path to output file for the distance matrix; if omitted, prints to stdout")
    parser_vcf2dist.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser_vcf2dist.add_argument("--ignoreHets", action="store_true", help="Ignore heterozygous loci (default: false)")
    parser_vcf2dist.add_argument("--onlyHets", action="store_true", help="Use only heterozygous loci (default: false)")
    parser_vcf2dist.add_argument("--ignoreMissing", action="store_true", help="Ignore missing loci (default: false)")
    parser_vcf2dist.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")
    parser_vcf2dist.add_argument("--tmpdir",default=None,help="Optional directory to store temporary file when reading from stdin")

    # Subcommand for VCF-based tree
    parser_vcf2tree = subparsers.add_parser("VCF2TREE", help="Compute tree from VCF")
    parser_vcf2tree.add_argument("vcf_file", nargs="?", help="Path to input VCF file (positional)")
    parser_vcf2tree.add_argument("-i", "--input", help="Optional input VCF file (overrides positional)")
    parser_vcf2tree.add_argument("-o", "--output", help="Path to output file for the tree; if omitted, prints to stdout")
    parser_vcf2tree.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser_vcf2tree.add_argument("--ignoreHets", action="store_true", help="Ignore heterozygous loci (default: false)")
    parser_vcf2tree.add_argument("--onlyHets", action="store_true", help="Use only heterozygous loci (default: false)")
    parser_vcf2tree.add_argument("--ignoreMissing", action="store_true", help="Ignore missing loci (default: false)")
    parser_vcf2tree.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")
    parser_vcf2tree.add_argument("--tmpdir",default=None,help="Optional directory to store temporary file when reading from stdin")
   
    # Subcommand for distance matrix to newick tree
    parser_dist2tree = subparsers.add_parser("DIST2TREE", help="Compute tree from distance matrix")
    parser_dist2tree.add_argument("dist_file", nargs="?", help="Path to input distance matrix file (positional)")
    parser_dist2tree.add_argument("-i", "--input", help="Optional input distance matrix file (overrides positional)")
    parser_dist2tree.add_argument("-o", "--output", help="Path to output file for the tree; if omitted, prints to stdout")
    parser_dist2tree.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")
   
    # Subcommand for FASTA-based distance matrix
    parser_fasta2dist = subparsers.add_parser("FASTA2DIST", help="Compute distance matrix from FASTA(s)")
    parser_fasta2dist.add_argument("fasta_file", nargs="?", help="Path to input fasta file (positional)")
    parser_fasta2dist.add_argument("-i", "--input", help="Optional input fasta file (overrides positional)")
    parser_fasta2dist.add_argument("-o", "--output", help="Path to output file for the distance matrix; if omitted, prints to stdout")
    parser_fasta2dist.add_argument("-k", "--kmerSize", type=int, default=4, help="Kmer size for D2S calculation (default: 4)")
    parser_fasta2dist.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser_fasta2dist.add_argument("-n", "--normalize", action="store_true", help="Use normalization")
    parser_fasta2dist.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")
   
    args, unknown = parser.parse_known_args()

    if args.version:
        print_version_from_jar(args.lib)
        return

    if args.command == "VCF2DIST":
        handle_vcf_tool(args, "VCF2DIST")
    
    elif args.command == "VCF2TREE":
        handle_vcf_tool(args, "VCF2TREE")
    
    elif args.command == "DIST2TREE":
        # Determine input source
        stdin_pipe = None
        input_dist = args.input if args.input else args.dist_file
        
        if input_dist == "-":
            print("Reading distance matrix from stdin...", file=sys.stderr)
            stdin_pipe = sys.stdin
        elif not input_dist:
            print("Error: no input dist file provided. Use positional argument, -i/--input, or pipe with -i -", file=sys.stderr)
            sys.exit(1)
        params = []
        if args.verbose:
            params.append("--verbose")
        params.append(input_dist)
        #print(params, file=sys.stderr)
        run_java_tool("DIST2TREE", params, args.lib, args.mem, args.output, args.verbose, stdin_pipe)
    
    elif args.command == "FASTA2DIST":
        input_fasta = args.input if args.input else args.fasta_file
        use_temp_input = False
        if input_fasta == "-":
            if sys.stdin.isatty():
                print("Error: -i - specified but no input is piped to stdin.", file=sys.stderr)
                sys.exit(1)
            print("Reading FASTA from stdin...", file=sys.stderr)
            temp_input = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta",
                                                    dir=args.tmpdir if hasattr(args, "tmpdir") and args.tmpdir else None)
            print(f"Storing stdin content in temporary file: {temp_input.name}", file=sys.stderr)
            for line in sys.stdin:
                temp_input.write(line)
            temp_input.close()
            input_fasta = temp_input.name
            use_temp_input = True
        elif not input_fasta:
            print("Error: no input FASTA file provided. Use positional argument, -i/--input, or pipe with -i -", file=sys.stderr)
            sys.exit(1)
        params = []
        if args.normalize:
            params.append("--normalize")
        if args.verbose:
            params.append("--verbose")
        params.extend(["-k", str(args.kmerSize)])
        params.extend(["-t", str(args.threads)])
        params.extend(["-i", input_fasta])
        run_java_tool("FASTA2DIST", params, args.lib, args.mem, args.output, args.verbose)
        if use_temp_input:
            try:
                os.unlink(input_fasta)
                print(f"Temporary input file {input_fasta} deleted.", file=sys.stderr)
            except Exception as e:
                print(f"Warning: Failed to delete temp file {input_fasta}: {e}", file=sys.stderr)
            
    else:
        parser.print_help()

def handle_vcf_tool(args, tool_name):
    input_vcf = args.input if args.input else args.vcf_file
    use_temp_input = False
    if input_vcf == "-":
        if sys.stdin.isatty():
            print("Error: -i - specified but no input is piped to stdin.", file=sys.stderr)
            sys.exit(1)
        print("Reading VCF from stdin...", file=sys.stderr)
        temp_input = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf",
                                                 dir=args.tmpdir if hasattr(args, "tmpdir") and args.tmpdir else None)
        print(f"Storing stdin content in temporary file: {temp_input.name}", file=sys.stderr)
        for line in sys.stdin:
            temp_input.write(line)
        temp_input.close()
        input_vcf = temp_input.name
        use_temp_input = True
    elif not input_vcf:
        print("Error: no input VCF file provided. Use positional argument, -i/--input, or pipe with -i -", file=sys.stderr)
        sys.exit(1)
    params = []
    if args.verbose:
        params.append("--verbose")
    if args.ignoreHets:
        params.append("--ignoreHets")
    if args.onlyHets:
        params.append("--onlyHets")
    if args.ignoreMissing:
        params.append("--ignoreMissing")
    params.extend(["-t", str(args.threads)])
    params.append(input_vcf)
    run_java_tool(tool_name, params, args.lib, args.mem, args.output, args.verbose)
    if use_temp_input:
        try:
            os.unlink(input_vcf)
            print(f"Temporary input file {input_vcf} deleted.", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Failed to delete temp file {input_vcf}: {e}", file=sys.stderr)

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
    print("\n📄 Citation:")
    print("Anestis Gkanogiannis (2016).")
    print("A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes.")
    print("BMC Bioinformatics 17, 311 (2016)")
    print("https://doi.org/10.1186/s12859-016-1186-3")
    print("https://github.com/gkanogiannis/fastreeR\n\n")

if __name__ == "__main__":
    main()
