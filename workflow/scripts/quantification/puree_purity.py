#!/usr/bin/env python3
"""
PUREE Tumor Purity Estimation
Estimates tumor purity (cancer cell fraction) from bulk gene expression data
"""

import sys
import os
import subprocess
import pandas as pd
import argparse

# Logging
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = sys.stderr

def check_puree_installation():
    """Check if PUREE is installed and accessible"""
    try:
        result = subprocess.run(
            ["python3", "-c", "import puree"],
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode == 0:
            cat("PUREE is installed\n")
            return True
        else:
            cat("PUREE is not installed. Attempting to use standalone script...\n")
            return False
    except Exception as e:
        cat(f"Error checking PUREE: {e}\n")
        return False

def run_puree_standalone(expression_file, output_file, gene_id_type="ENSEMBL"):
    """
    Run PUREE using standalone predict_purity.py script
    """
    cat("Running PUREE using standalone script...\n")

    # Check if predict_purity.py exists
    puree_script = "predict_purity.py"
    if not os.path.exists(puree_script):
        # Try to find in common locations
        possible_paths = [
            "/usr/local/bin/predict_purity.py",
            os.path.expanduser("~/PUREE/predict_purity.py"),
            "./PUREE/predict_purity.py"
        ]
        for path in possible_paths:
            if os.path.exists(path):
                puree_script = path
                break
        else:
            raise FileNotFoundError(
                "PUREE predict_purity.py script not found. "
                "Please install PUREE from https://github.com/skandlab/PUREE"
            )

    cat(f"Using PUREE script: {puree_script}\n")

    # Build command
    cmd = [
        "python3",
        puree_script,
        "--data_path", expression_file,
        "--output", output_file,
        "--gene_identifier_type", gene_id_type
    ]

    cat(f"Command: {' '.join(cmd)}\n")

    # Run PUREE
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False
    )

    if result.returncode != 0:
        cat(f"PUREE error:\n{result.stderr}\n")
        raise RuntimeError(f"PUREE failed with return code {result.returncode}")

    cat("PUREE completed successfully\n")
    cat(result.stdout)

    return output_file

def run_puree_api(expression_file, output_file, gene_id_type="ENSEMBL"):
    """
    Run PUREE using Python API
    """
    try:
        from puree import PUREE

        cat("Running PUREE using Python API...\n")

        # Read expression data
        cat(f"Reading expression data from {expression_file}\n")
        expr_df = pd.read_csv(expression_file, sep='\t', index_col=0)

        cat(f"Expression data dimensions: {expr_df.shape[0]} samples x {expr_df.shape[1]} genes\n")

        # Check orientation - PUREE expects samples as rows, genes as columns
        # Our count matrix typically has genes as rows, samples as columns
        cat("Transposing expression data (PUREE expects samples as rows)...\n")
        expr_df_transposed = expr_df.transpose()

        # Create temporary file for PUREE
        temp_file = snakemake.temp[0] if len(snakemake.temp) > 0 else "temp_expression.tsv"
        expr_df_transposed.to_csv(temp_file, sep='\t')

        cat(f"Transposed data dimensions: {expr_df_transposed.shape[0]} samples x {expr_df_transposed.shape[1]} genes\n")

        # Run PUREE
        cat("Estimating tumor purity...\n")
        p = PUREE()
        result = p.get_output(temp_file, gene_id_type)

        # Get purities
        purities_df = result["output"]
        logs = result["logs"]

        cat("PUREE logs:\n")
        cat(logs)

        # Write output
        purities_df.to_csv(output_file, sep='\t', index=False)
        cat(f"Tumor purities written to {output_file}\n")

        # Clean up temp file
        if os.path.exists(temp_file):
            os.remove(temp_file)

        return output_file

    except ImportError:
        cat("PUREE Python API not available, falling back to standalone script...\n")
        return run_puree_standalone(expression_file, output_file, gene_id_type)

def main():
    # Get parameters from Snakemake
    expression_file = snakemake.input["expression"]
    output_file = snakemake.output["purities"]
    gene_id_type = snakemake.params.get("gene_id_type", "ENSEMBL")

    cat("=" * 78 + "\n")
    cat("PUREE Tumor Purity Estimation\n")
    cat("=" * 78 + "\n\n")

    cat(f"Input expression file: {expression_file}\n")
    cat(f"Output purities file: {output_file}\n")
    cat(f"Gene ID type: {gene_id_type}\n\n")

    # Check if input file exists
    if not os.path.exists(expression_file):
        raise FileNotFoundError(f"Expression file not found: {expression_file}")

    # Try to determine gene ID type from expression data
    cat("Examining expression data...\n")
    expr_df = pd.read_csv(expression_file, sep='\t', nrows=5, index_col=0)
    gene_ids = expr_df.index.tolist()

    cat(f"Sample gene IDs: {gene_ids[:5]}\n")

    # Auto-detect gene ID type if not specified
    if gene_id_type == "auto":
        # Check if gene IDs look like ENSEMBL (ENSxxxxx) or HGNC (letters)
        if any(gene_id.startswith("ENS") for gene_id in gene_ids):
            gene_id_type = "ENSEMBL"
            cat("Auto-detected gene ID type: ENSEMBL\n")
        else:
            gene_id_type = "HGNC"
            cat("Auto-detected gene ID type: HGNC\n")

    cat(f"\nUsing gene ID type: {gene_id_type}\n")

    # Run PUREE
    try:
        # First try Python API
        from puree import PUREE
        cat("PUREE Python API found, using API mode...\n\n")
        run_puree_api(expression_file, output_file, gene_id_type)
    except ImportError:
        # Fall back to standalone script
        cat("PUREE Python API not found, using standalone script...\n\n")
        run_puree_standalone(expression_file, output_file, gene_id_type)

    # Read and display results
    cat("\n" + "=" * 78 + "\n")
    cat("Results Summary\n")
    cat("=" * 78 + "\n\n")

    purities_df = pd.read_csv(output_file, sep='\t')
    cat(f"Tumor purity estimated for {len(purities_df)} samples\n\n")

    # Rename columns for clarity
    if len(purities_df.columns) == 1:
        purities_df.columns = ['tumor_purity']
    else:
        # First column is purity
        purities_df = purities_df.iloc[:, [0]]
        purities_df.columns = ['tumor_purity']

    cat("Tumor Purity Statistics:\n")
    cat(f"  Mean: {purities_df['tumor_purity'].mean():.4f}\n")
    cat(f"  Median: {purities_df['tumor_purity'].median():.4f}\n")
    cat(f"  Min: {purities_df['tumor_purity'].min():.4f}\n")
    cat(f"  Max: {purities_df['tumor_purity'].max():.4f}\n")
    cat(f"  Std: {purities_df['tumor_purity'].std():.4f}\n\n")

    cat("Sample purities:\n")
    for idx, row in purities_df.iterrows():
        cat(f"  Sample {idx}: {row['tumor_purity']:.4f}\n")

    cat("\n" + "=" * 78 + "\n")
    cat("PUREE Analysis Complete!\n")
    cat("=" * 78 + "\n")

if __name__ == "__main__":
    main()
