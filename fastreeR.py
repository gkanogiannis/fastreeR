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
import re
import json
import gzip
import shutil
import tempfile
from contextlib import contextmanager

FASTREER_VERSION = "2.1.3"

# Default paths for BioFM resources (can be overridden via environment variables or CLI args)
DEFAULT_BIOFM_MODEL = os.environ.get("BIOFM_MODEL", "m42-health/BioFM-265M")
DEFAULT_REFERENCE_GENOME = os.environ.get("BIOFM_REFERENCE_GENOME", None)
DEFAULT_GENE_ANNOTATION = os.environ.get("BIOFM_GENE_ANNOTATION", None)

# Determine JAR directory
JAR_DIR = os.environ.get("FASTREER_JAR_DIR") or os.path.join(os.path.dirname(__file__), "inst/java")
MEM_MB = "256"
JAVA_PARAMS = ["-Djava.awt.headless", "-XX:+UseG1GC", "-XX:+UseStringDeduplication", "-Xmx"+str(MEM_MB)+"M"]
MAIN_CLASS="com.gkano.bioinfo.javautils.JavaUtils"

def check_java_version(min_major=11):
    try:
        result = subprocess.run(["java", "-version"], stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        version_output = result.stderr.splitlines()[0] if result.stderr else result.stdout.splitlines()[0]
        # Match version patterns like "11.0.2", "1.8.0", or "17"
        match = re.search(r'version\s+"?(\d+(?:\.\d+)*)"?', version_output)
        if match:
            version_str = match.group(1)
            major_version = int(version_str.split('.')[0]) if not version_str.startswith("1.") else int(version_str.split('.')[1])
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

def is_gzipped(filepath):
    """Check if a file is gzip compressed based on magic bytes."""
    try:
        with open(filepath, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except (IOError, OSError):
        return False


def get_decompressed_suffix(filepath):
    """Get the appropriate suffix for a decompressed file."""
    base = os.path.basename(filepath)
    # Remove .gz extension if present
    if base.endswith('.gz'):
        base = base[:-3]
    # Determine suffix based on remaining extension
    if base.endswith('.vcf'):
        return '.vcf'
    elif base.endswith('.gff') or base.endswith('.gff3'):
        return '.gff3'
    elif base.endswith('.fa') or base.endswith('.fasta') or base.endswith('.fna'):
        return '.fasta'
    else:
        # Default: use whatever extension remains
        _, ext = os.path.splitext(base)
        return ext if ext else '.txt'


@contextmanager
def decompress_if_gzipped(filepath, verbose=False):
    """
    Context manager that decompresses a gzipped file to a temporary file if needed.
    Yields the path to use (original or temporary decompressed file).
    Cleans up temporary file on exit.
    """
    if not is_gzipped(filepath):
        yield filepath
        return

    suffix = get_decompressed_suffix(filepath)
    temp_fd, temp_path = tempfile.mkstemp(suffix=suffix)

    try:
        if verbose:
            print(f"[VCF2EMB] Decompressing {filepath} to temporary file...", file=sys.stderr)

        with gzip.open(filepath, 'rb') as f_in:
            with os.fdopen(temp_fd, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        if verbose:
            print(f"[VCF2EMB] Decompressed to: {temp_path}", file=sys.stderr)

        yield temp_path
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            os.unlink(temp_path)
            if verbose:
                print(f"[VCF2EMB] Cleaned up temporary file: {temp_path}", file=sys.stderr)


@contextmanager
def prepare_input_files(vcf_path, reference_genome, gene_annotation, verbose=False):
    """
    Context manager that prepares all input files, decompressing gzipped files as needed.
    Yields a tuple of (vcf_path, reference_genome_path, gene_annotation_path).
    """
    with decompress_if_gzipped(vcf_path, verbose) as vcf_ready:
        with decompress_if_gzipped(reference_genome, verbose) as ref_ready:
            with decompress_if_gzipped(gene_annotation, verbose) as ann_ready:
                yield vcf_ready, ref_ready, ann_ready


def generate_embeddings_biofm(vcf_path, output_path, reference_genome, gene_annotation,
                               model_name=DEFAULT_BIOFM_MODEL, output_format="TSV",
                               variant_key_format="CHROM_POS_REF_ALT", max_variants=100,
                               device=None, verbose=False):
    """
    Generate variant embeddings from a VCF file using BioFM-265M model.

    Supports gzipped input files (.vcf.gz, .fa.gz, .fasta.gz, .gff.gz, .gff3.gz).
    Gzipped files are automatically decompressed to temporary files during processing.

    Args:
        vcf_path: Path to input VCF file (can be gzipped)
        output_path: Path to output embeddings file (or None for stdout)
        reference_genome: Path to reference genome FASTA file (can be gzipped)
        gene_annotation: Path to gene annotation GFF3 file (can be gzipped)
        model_name: HuggingFace model name or local path
        output_format: Output format - "TSV" or "HUGGINGFACE"
        variant_key_format: How to format variant keys - "CHROM_POS", "CHROM_POS_REF_ALT", or "VCF_ID"
        max_variants: Maximum number of variants to process (None for all)
        device: Device to use ("cuda", "cpu", or None for auto)
        verbose: Print progress messages
    """
    try:
        import torch
        from biofm_eval import AnnotatedModel, AnnotationTokenizer, Embedder, VCFConverter
    except ImportError as e:
        print(f"Error: Required packages not installed. Please install with:", file=sys.stderr)
        print(f"  pip install biofm-eval torch", file=sys.stderr)
        print(f"Import error: {e}", file=sys.stderr)
        sys.exit(1)

    if not reference_genome or not os.path.exists(reference_genome):
        print(f"Error: Reference genome file required. Provide via --reference or BIOFM_REFERENCE_GENOME env var.", file=sys.stderr)
        sys.exit(1)

    if not gene_annotation or not os.path.exists(gene_annotation):
        print(f"Error: Gene annotation file required. Provide via --annotation or BIOFM_GENE_ANNOTATION env var.", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"[VCF2EMB] Loading BioFM model: {model_name}", file=sys.stderr)

    # Determine device
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    if verbose:
        print(f"[VCF2EMB] Using device: {device}", file=sys.stderr)

    # Load model and tokenizer
    model = AnnotatedModel.from_pretrained(
        model_name,
        torch_dtype=torch.bfloat16 if device == "cuda" else torch.float32,
    )
    model = model.to(device)
    tokenizer = AnnotationTokenizer.from_pretrained(model_name)

    # Initialize embedder
    embedder = Embedder(model, tokenizer)

    # Prepare input files (decompress gzipped files if needed)
    with prepare_input_files(vcf_path, reference_genome, gene_annotation, verbose) as (vcf_ready, ref_ready, ann_ready):
        if verbose:
            print(f"[VCF2EMB] Setting up VCF converter with:", file=sys.stderr)
            print(f"  Reference genome: {ref_ready}", file=sys.stderr)
            print(f"  Gene annotation: {ann_ready}", file=sys.stderr)

        # Set up VCF converter
        vcf_converter = VCFConverter(
            gene_annotation_path=ann_ready,
            reference_genome_path=ref_ready
        )

        if verbose:
            print(f"[VCF2EMB] Processing VCF file: {vcf_ready}", file=sys.stderr)

        # Convert VCF to annotated dataset
        convert_kwargs = {"vcf_path": vcf_ready}
        if max_variants is not None:
            convert_kwargs["max_variants"] = max_variants

        annotated_dataset = vcf_converter.vcf_to_annotated_dataset(**convert_kwargs)

    if verbose:
        print(f"[VCF2EMB] Extracting embeddings for {len(annotated_dataset)} variants...", file=sys.stderr)

    # Extract embeddings
    result = embedder.get_dataset_embeddings(annotated_dataset)
    embeddings = result["embeddings"]

    if verbose:
        print(f"[VCF2EMB] Generated embeddings with shape: {embeddings.shape}", file=sys.stderr)

    # Get variant IDs from the dataset
    variant_ids = []
    for item in annotated_dataset:
        # Extract variant info - format depends on VCFConverter output
        if hasattr(item, 'variant_id'):
            vid = item.variant_id
        elif isinstance(item, dict) and 'variant_id' in item:
            vid = item['variant_id']
        else:
            # Fallback: try to construct from available fields
            chrom = getattr(item, 'chrom', None) or item.get('chrom', 'unknown')
            pos = getattr(item, 'pos', None) or item.get('pos', 0)
            ref = getattr(item, 'ref', None) or item.get('ref', 'N')
            alt = getattr(item, 'alt', None) or item.get('alt', 'N')

            if variant_key_format == "CHROM_POS":
                vid = f"{chrom}:{pos}"
            elif variant_key_format == "VCF_ID":
                vid = getattr(item, 'id', None) or item.get('id', f"{chrom}:{pos}:{ref}:{alt}")
            else:  # CHROM_POS_REF_ALT
                vid = f"{chrom}:{pos}:{ref}:{alt}"

        variant_ids.append(vid)

    # Write output
    out_stream = open(output_path, "w") if output_path else sys.stdout
    try:
        if output_format.upper() == "HUGGINGFACE":
            # HuggingFace JSON format
            output_data = {
                "model_name": model_name,
                "embedding_dim": embeddings.shape[1],
                "variant_key_format": variant_key_format,
                "num_variants": len(variant_ids),
                "variants": [
                    {"id": vid, "embedding": emb.tolist()}
                    for vid, emb in zip(variant_ids, embeddings)
                ]
            }
            json.dump(output_data, out_stream, indent=2)
            out_stream.write("\n")
        else:
            # TSV format
            embedding_dim = embeddings.shape[1]
            header = "#VARIANT_ID\t" + "\t".join([f"DIM_{i}" for i in range(embedding_dim)])
            out_stream.write(header + "\n")

            for vid, emb in zip(variant_ids, embeddings):
                line = vid + "\t" + "\t".join([f"{v:.6f}" for v in emb])
                out_stream.write(line + "\n")

        if verbose:
            if output_path:
                print(f"[VCF2EMB] Wrote {len(variant_ids)} variant embeddings to {output_path}", file=sys.stderr)
            else:
                print(f"[VCF2EMB] Wrote {len(variant_ids)} variant embeddings to stdout", file=sys.stderr)

    finally:
        if output_path:
            out_stream.close()


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
        # p.add_argument("--ignoreHets", action="store_true", help="Ignore heterozygous loci (default: false)")
        # p.add_argument("--onlyHets", action="store_true", help="Use only heterozygous loci (default: false)")
        # p.add_argument("--ignoreMissing", action="store_true", help="Ignore missing loci (default: false)")
        p.add_argument("-v", "--verbose", action="store_true", help="Print progress messages on stderr (default: false)")

    def add_embedding_args(p):
        p.add_argument("-e", "--embeddings", type=str, default=None,
                       help="Path to variant embeddings file for embedding-based distance calculation")
        p.add_argument("--embeddings-format", type=str, default=None, choices=["TSV", "HUGGINGFACE"],
                       help="Embeddings file format: TSV or HUGGINGFACE (auto-detected if not specified)")
        p.add_argument("--variant-key", type=str, default="CHROM_POS_REF_ALT",
                       choices=["CHROM_POS", "CHROM_POS_REF_ALT", "VCF_ID"],
                       help="Variant key format for embedding lookup (default: CHROM_POS_REF_ALT)")

    # Subcommand for VCF-based distance matrix
    parser_vcf2dist = subparsers.add_parser("VCF2DIST", help="Compute distance matrix from VCF(s)")
    add_common_input_output(parser_vcf2dist)
    add_common_vcf_args(parser_vcf2dist)
    add_embedding_args(parser_vcf2dist)

    # Subcommand for VCF-based tree
    parser_vcf2tree = subparsers.add_parser("VCF2TREE", help="Compute tree from VCF(s)")
    add_common_input_output(parser_vcf2tree)
    add_common_vcf_args(parser_vcf2tree)
    add_embedding_args(parser_vcf2tree)
    parser_vcf2tree.add_argument("-b", "--bootstrap", type=int, default=0,
                                 help="Number of bootstrap replicates to perform (default: 0, no bootstrapping)")

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

    # Subcommand for generating variant embeddings from VCF using BioFM
    parser_vcf2emb = subparsers.add_parser("VCF2EMB",
        help="Generate variant embeddings from VCF using BioFM genomic language model",
        description="Generate variant embeddings from a VCF file using the BioFM-265M genomic "
                    "language model. Requires biofm-eval package (pip install biofm-eval) and "
                    "reference genome/annotation files.")
    parser_vcf2emb.add_argument("input_file", nargs="?", help="Input VCF file")
    parser_vcf2emb.add_argument("-i", "--input", dest="named_input", help="Input VCF file (overrides positional)")
    parser_vcf2emb.add_argument("-o", "--output", help="Output embeddings file (default: stdout)")
    parser_vcf2emb.add_argument("-r", "--reference", type=str, default=DEFAULT_REFERENCE_GENOME,
                                help="Path to reference genome FASTA file (or set BIOFM_REFERENCE_GENOME env var)")
    parser_vcf2emb.add_argument("-a", "--annotation", type=str, default=DEFAULT_GENE_ANNOTATION,
                                help="Path to gene annotation GFF3 file (or set BIOFM_GENE_ANNOTATION env var)")
    parser_vcf2emb.add_argument("-m", "--model", type=str, default=DEFAULT_BIOFM_MODEL,
                                help=f"HuggingFace model name or local path (default: {DEFAULT_BIOFM_MODEL})")
    parser_vcf2emb.add_argument("-f", "--format", type=str, default="TSV", choices=["TSV", "HUGGINGFACE"],
                                help="Output format: TSV or HUGGINGFACE JSON (default: TSV)")
    parser_vcf2emb.add_argument("--variant-key", type=str, default="CHROM_POS_REF_ALT",
                                choices=["CHROM_POS", "CHROM_POS_REF_ALT", "VCF_ID"],
                                help="Variant key format in output (default: CHROM_POS_REF_ALT)")
    parser_vcf2emb.add_argument("--max-variants", type=int, default=None,
                                help="Maximum number of variants to process (default: all)")
    parser_vcf2emb.add_argument("--device", type=str, default=None, choices=["cuda", "cpu"],
                                help="Device for model inference (default: auto-detect)")
    parser_vcf2emb.add_argument("-v", "--verbose", action="store_true",
                                help="Print progress messages on stderr (default: false)")

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

            cmd = ["java"] + JAVA_PARAMS + ["-cp", classpath, MAIN_CLASS]
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
        # print version info for wrapper + backend and exit
        print_version(args.lib)
        sys.exit(0)

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
        if args.verbose or args.extraVerbose:
            params.append("--verbose")
        # if args.ignoreHets: params.append("--ignoreHets")
        # if args.onlyHets: params.append("--onlyHets")
        # if args.ignoreMissing: params.append("--ignoreMissing")
        params.extend(["-t", str(int(args.threads))])
        # forward bootstrap only when requesting tree generation
        if args.command == "VCF2TREE" and getattr(args, 'bootstrap', 0) and int(args.bootstrap) > 0:
            params.extend(["--bootstrap", str(int(args.bootstrap))])
        # forward embedding options if provided
        if getattr(args, 'embeddings', None):
            params.extend(["-e", args.embeddings])
        if getattr(args, 'embeddings_format', None):
            params.extend(["--embeddings-format", args.embeddings_format])
        if getattr(args, 'variant_key', None) and args.variant_key != "CHROM_POS_REF_ALT":
            params.extend(["--variant-key", args.variant_key])
        for f in input_files:
            params.extend(["-i", f])
        run_java_tool(args.command, params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)

    elif args.command == "DIST2TREE":
        input_file = args.named_input or args.input_file
        if not input_file:
            print("Error: No input distance matrix provided.", file=sys.stderr)
            sys.exit(1)
        params = []
        if args.verbose or args.extraVerbose:
            params.append("--verbose")
        params.append(input_file)
        run_java_tool("DIST2TREE", params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)

    elif args.command == "FASTA2DIST":
        input_files = resolve_inputs(args)
        params = []
        if args.verbose or args.extraVerbose:
            params.append("--verbose")
        if args.normalize:
            params.append("--normalize")
        params.extend(["-k", str(args.kmerSize), "-t", str(int(args.threads))])
        for f in input_files:
            params.extend(["-i", f])
        run_java_tool("FASTA2DIST", params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)

    elif args.command == "VCF2EMB":
        input_file = args.named_input or args.input_file
        if not input_file:
            print("Error: No input VCF file provided.", file=sys.stderr)
            sys.exit(1)
        generate_embeddings_biofm(
            vcf_path=input_file,
            output_path=args.output,
            reference_genome=args.reference,
            gene_annotation=args.annotation,
            model_name=args.model,
            output_format=args.format,
            variant_key_format=args.variant_key,
            max_variants=args.max_variants,
            device=args.device,
            verbose=args.verbose or args.extraVerbose
        )

    else:
        print("Unknown command", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

def get_backend_version_from_jar(jar_dir):
    try:
        jars = [os.path.join(jar_dir, f) for f in os.listdir(jar_dir) if f.endswith(".jar")]
    except FileNotFoundError:
        return None

    for jar in jars:
        try:
            with zipfile.ZipFile(jar, 'r') as zipf:
                for name in zipf.namelist():
                    if "pom.properties" in name:
                        with zipf.open(name) as props:
                            for line in props:
                                decoded = line.decode().strip()
                                if decoded.startswith("version="):
                                    return decoded.split("=", 1)[1]
        except Exception:
            pass
    return None

def print_version(jar_dir):
    backend_version = get_backend_version_from_jar(jar_dir)
    if backend_version:
        print(f"fastreeR version: {FASTREER_VERSION} (backend: {backend_version})")
    else:
        print(f"fastreeR version: {FASTREER_VERSION} (backend: unknown)")
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
