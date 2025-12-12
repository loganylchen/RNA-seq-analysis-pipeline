#!/usr/bin/env python3
"""
Convert GTF file to refFlat format for Picard tools.

refFlat format is a tab-separated format with the following columns:
1. geneName - Gene name
2. name - Transcript name
3. chrom - Chromosome name
4. strand - "+" or "-"
5. txStart - Transcription start position (0-based)
6. txEnd - Transcription end position (1-based)
7. cdsStart - Coding region start (0-based)
8. cdsEnd - Coding region end (1-based)
9. exonCount - Number of exons
10. exonStarts - Comma-separated list of exon start positions (0-based)
11. exonEnds - Comma-separated list of exon end positions (1-based)
"""

import sys
from collections import defaultdict


def parse_attributes(attr_string):
    """Parse GTF attributes into a dictionary."""
    attributes = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        try:
            key, value = item.split(None, 1)
            attributes[key] = value.strip('"')
        except ValueError:
            continue
    return attributes


def gtf_to_refflat(gtf_file, output_file):
    """Convert GTF to refFlat format."""
    # Dictionary to store transcripts
    transcripts = defaultdict(lambda: {
        'chrom': '',
        'strand': '',
        'gene_name': '',
        'transcript_id': '',
        'exons': [],
        'cds': [],
    })
    
    # Read GTF file
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            attr_dict = parse_attributes(attributes)
            
            # Skip if no transcript_id
            if 'transcript_id' not in attr_dict:
                continue
            
            transcript_id = attr_dict['transcript_id']
            gene_name = attr_dict.get('gene_name', attr_dict.get('gene_id', 'unknown'))
            
            # Initialize transcript info
            if not transcripts[transcript_id]['chrom']:
                transcripts[transcript_id]['chrom'] = chrom
                transcripts[transcript_id]['strand'] = strand
                transcripts[transcript_id]['gene_name'] = gene_name
                transcripts[transcript_id]['transcript_id'] = transcript_id
            
            # Collect exons and CDS
            if feature == 'exon':
                transcripts[transcript_id]['exons'].append((int(start) - 1, int(end)))  # Convert to 0-based
            elif feature == 'CDS':
                transcripts[transcript_id]['cds'].append((int(start) - 1, int(end)))  # Convert to 0-based
    
    # Write refFlat format
    with open(output_file, 'w') as out:
        for transcript_id, data in sorted(transcripts.items()):
            if not data['exons']:
                continue
            
            # Sort exons and CDS
            exons = sorted(data['exons'], key=lambda x: x[0])
            cds = sorted(data['cds'], key=lambda x: x[0])
            
            # Get transcript boundaries
            tx_start = exons[0][0]
            tx_end = exons[-1][1]
            
            # Get CDS boundaries (if any)
            if cds:
                cds_start = cds[0][0]
                cds_end = cds[-1][1]
            else:
                # Non-coding transcript
                cds_start = tx_end
                cds_end = tx_end
            
            # Format exon starts and ends
            exon_starts = ','.join(str(e[0]) for e in exons) + ','
            exon_ends = ','.join(str(e[1]) for e in exons) + ','
            
            # Write refFlat line
            out.write('\t'.join([
                data['gene_name'],
                transcript_id,
                data['chrom'],
                data['strand'],
                str(tx_start),
                str(tx_end),
                str(cds_start),
                str(cds_end),
                str(len(exons)),
                exon_starts,
                exon_ends
            ]) + '\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.gtf> <output.refFlat>", file=sys.stderr)
        sys.exit(1)
    
    gtf_file = sys.argv[1]
    output_file = sys.argv[2]
    
    gtf_to_refflat(gtf_file, output_file)
    print(f"Converted {gtf_file} to {output_file}", file=sys.stderr)
