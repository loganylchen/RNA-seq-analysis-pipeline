#!/usr/bin/env python3
"""
Analyze SpliceTools results (RIMedley, SEMedley, SpliceCompare).
Creates summary reports and visualizations.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Parse Snakemake inputs
splicetools_dir = snakemake.input["splicetools_dir"]
output_dir = snakemake.output["output_dir"]
fdr_threshold = snakemake.params["fdr_threshold"]

# Create output directory
os.makedirs(output_dir, exist_ok=True)

def read_rmats_format(file_path):
    """Read rMATS format files from SpliceTools."""
    if not os.path.exists(file_path):
        return None
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def parse_ri_medley(splicetools_dir):
    """Parse RIMedley results."""
    results = {}
    
    # Look for RI output files
    ri_files = [
        "RI.sig.txt",
        "RI_Medley_Summary.txt",
        "RI.MATS.JCEC.txt"
    ]
    
    for fname in ri_files:
        fpath = os.path.join(splicetools_dir, fname)
        if os.path.exists(fpath):
            try:
                df = pd.read_csv(fpath, sep='\t')
                results[fname] = df
                print(f"Loaded {fname}: {len(df)} entries")
            except Exception as e:
                print(f"Could not read {fname}: {e}")
    
    return results

def parse_se_medley(splicetools_dir):
    """Parse SEMedley results."""
    results = {}
    
    # Look for SE output files
    se_files = [
        "SE.sig.txt",
        "SE_Medley_Summary.txt",
        "SE.MATS.JCEC.txt"
    ]
    
    for fname in se_files:
        fpath = os.path.join(splicetools_dir, fname)
        if os.path.exists(fpath):
            try:
                df = pd.read_csv(fpath, sep='\t')
                results[fname] = df
                print(f"Loaded {fname}: {len(df)} entries")
            except Exception as e:
                print(f"Could not read {fname}: {e}")
    
    return results

def parse_splice_compare(splicetools_dir):
    """Parse SpliceCompare results."""
    results = {}
    
    # Look for comparison files
    compare_files = [f for f in os.listdir(splicetools_dir) 
                     if 'compare' in f.lower() or 'summary' in f.lower()]
    
    for fname in compare_files:
        if fname.endswith('.txt'):
            fpath = os.path.join(splicetools_dir, fname)
            try:
                df = pd.read_csv(fpath, sep='\t')
                results[fname] = df
                print(f"Loaded {fname}: {len(df)} entries")
            except Exception as e:
                print(f"Could not read {fname}: {e}")
    
    return results

def plot_ri_summary(ri_results, output_dir, fdr_threshold):
    """Create visualizations for retained intron analysis."""
    
    if "RI.MATS.JCEC.txt" in ri_results:
        df = ri_results["RI.MATS.JCEC.txt"]
        
        # Filter significant
        if 'FDR' in df.columns:
            df['Significant'] = df['FDR'] < fdr_threshold
            sig_df = df[df['Significant']]
            
            # Summary stats
            summary = {
                'Total RI Events': len(df),
                'Significant RI Events': len(sig_df),
                'Increased Retention': len(sig_df[sig_df['IncLevelDifference'] > 0]),
                'Decreased Retention': len(sig_df[sig_df['IncLevelDifference'] < 0])
            }
            
            pd.DataFrame([summary]).to_csv(
                os.path.join(output_dir, 'RI_summary.csv'), index=False
            )
            
            # Volcano plot
            fig, ax = plt.subplots(figsize=(8, 6))
            
            colors = []
            for _, row in df.iterrows():
                if row['FDR'] < fdr_threshold:
                    if row['IncLevelDifference'] > 0:
                        colors.append('#D55E00')  # Increased
                    else:
                        colors.append('#0072B2')  # Decreased
                else:
                    colors.append('lightgray')
            
            ax.scatter(df['IncLevelDifference'], -np.log10(df['FDR']), 
                      c=colors, alpha=0.6, s=20)
            ax.axhline(-np.log10(fdr_threshold), color='gray', 
                      linetype='--', linewidth=1)
            ax.axvline(0, color='gray', linetype='--', linewidth=1)
            ax.set_xlabel('IncLevel Difference (dPSI)', fontsize=12)
            ax.set_ylabel('-log10(FDR)', fontsize=12)
            ax.set_title('Retained Intron Events - Volcano Plot', 
                        fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'RI_volcano.pdf'))
            plt.close()
            
            print(f"RI Summary: {summary}")

def plot_se_summary(se_results, output_dir, fdr_threshold):
    """Create visualizations for skipped exon analysis."""
    
    if "SE.MATS.JCEC.txt" in se_results:
        df = se_results["SE.MATS.JCEC.txt"]
        
        # Filter significant
        if 'FDR' in df.columns:
            df['Significant'] = df['FDR'] < fdr_threshold
            sig_df = df[df['Significant']]
            
            # Summary stats
            summary = {
                'Total SE Events': len(df),
                'Significant SE Events': len(sig_df),
                'Increased Inclusion': len(sig_df[sig_df['IncLevelDifference'] > 0]),
                'Increased Skipping': len(sig_df[sig_df['IncLevelDifference'] < 0])
            }
            
            pd.DataFrame([summary]).to_csv(
                os.path.join(output_dir, 'SE_summary.csv'), index=False
            )
            
            # Volcano plot
            fig, ax = plt.subplots(figsize=(8, 6))
            
            colors = []
            for _, row in df.iterrows():
                if row['FDR'] < fdr_threshold:
                    if row['IncLevelDifference'] > 0:
                        colors.append('#D55E00')  # Increased inclusion
                    else:
                        colors.append('#0072B2')  # Increased skipping
                else:
                    colors.append('lightgray')
            
            ax.scatter(df['IncLevelDifference'], -np.log10(df['FDR']), 
                      c=colors, alpha=0.6, s=20)
            ax.axhline(-np.log10(fdr_threshold), color='gray', 
                      linetype='--', linewidth=1)
            ax.axvline(0, color='gray', linetype='--', linewidth=1)
            ax.set_xlabel('IncLevel Difference (dPSI)', fontsize=12)
            ax.set_ylabel('-log10(FDR)', fontsize=12)
            ax.set_title('Skipped Exon Events - Volcano Plot', 
                        fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'SE_volcano.pdf'))
            plt.close()
            
            print(f"SE Summary: {summary}")

def create_combined_summary(ri_results, se_results, output_dir):
    """Create combined summary across analyses."""
    
    summary_data = []
    
    for event_type, results in [('RI', ri_results), ('SE', se_results)]:
        key = f"{event_type}.MATS.JCEC.txt"
        if key in results:
            df = results[key]
            if 'FDR' in df.columns:
                summary_data.append({
                    'Event Type': event_type,
                    'Total Events': len(df),
                    'Significant (FDR<0.05)': len(df[df['FDR'] < 0.05]),
                    'Mean dPSI': df['IncLevelDifference'].mean(),
                    'Median dPSI': df['IncLevelDifference'].median()
                })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(os.path.join(output_dir, 'combined_summary.csv'), 
                         index=False)
        
        # Bar plot
        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(summary_df))
        width = 0.35
        
        ax.bar(x - width/2, summary_df['Total Events'], width, 
               label='Total', alpha=0.8)
        ax.bar(x + width/2, summary_df['Significant (FDR<0.05)'], width, 
               label='Significant', alpha=0.8)
        
        ax.set_xlabel('Event Type', fontsize=12)
        ax.set_ylabel('Number of Events', fontsize=12)
        ax.set_title('SpliceTools Event Summary', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(summary_df['Event Type'])
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'combined_summary.pdf'))
        plt.close()

# Main analysis
print("Analyzing SpliceTools results...")

# Parse results
ri_results = parse_ri_medley(splicetools_dir)
se_results = parse_se_medley(splicetools_dir)
compare_results = parse_splice_compare(splicetools_dir)

# Create visualizations
if ri_results:
    plot_ri_summary(ri_results, output_dir, fdr_threshold)

if se_results:
    plot_se_summary(se_results, output_dir, fdr_threshold)

# Combined summary
create_combined_summary(ri_results, se_results, output_dir)

print(f"\nSpliceTools analysis complete! Results saved to: {output_dir}")
