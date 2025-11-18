


import sys
from collections import defaultdict

sys.stderr = open(snakemake.log[0], "w")


def parse_gtf(gtf_file):
    """Parse GTF file and group exons by transcript."""
    transcripts = defaultdict(lambda: {
        'chrom': None,
        'strand': None,
        'exons': [],
        'gene_id': None,
        'transcript_id': None,
    })
    
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            if feature != 'exon':
                continue
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, val = attr.split(' ', 1)
                    attr_dict[key] = val.strip('"')
            
            transcript_id = attr_dict.get('transcript_id')
            if not transcript_id:
                continue
            
            transcripts[transcript_id]['chrom'] = chrom
            transcripts[transcript_id]['strand'] = strand
            transcripts[transcript_id]['gene_id'] = attr_dict.get('gene_id')
            transcripts[transcript_id]['transcript_id'] = transcript_id
            transcripts[transcript_id]['exons'].append((int(start) - 1, int(end)))
    
    return transcripts


def gtf_to_bed12(gtf_file, bed_file):
    """Convert GTF to BED12."""
    transcripts = parse_gtf(gtf_file)
    
    with open(bed_file, 'w') as f:
        for transcript_id, data in sorted(transcripts.items()):
            if not data['exons']:
                continue
            
            exons = sorted(data['exons'])
            chrom_start = exons[0][0]
            chrom_end = exons[-1][1]
            
            block_count = len(exons)
            block_sizes = ','.join(str(e[1] - e[0]) for e in exons)
            block_starts = ','.join(str(e[0] - chrom_start) for e in exons)
            
            bed_line = [
                data['chrom'],
                str(chrom_start),
                str(chrom_end),
                transcript_id,
                '0',
                data['strand'],
                str(chrom_start),
                str(chrom_end),
                '0',
                str(block_count),
                block_sizes,
                block_starts
            ]
            
            f.write('\t'.join(bed_line) + '\n')




    
gtf_to_bed12(snakemake.input[0], snakemake.output[0])
