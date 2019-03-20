import argparse
import os
import subprocess

import bam_utils
import annotation

parser = argparse.ArgumentParser()

parser.add_argument('input_bam', type=str,
        help='input bam')
parser.add_argument('input_annotated_vaf', type=str,
        help='Annotated vaf file that is output from annotation station')
parser.add_argument('reference_fasta', type=str,
        help='Reference to use')

parser.add_argument('--input-header', action='store_true',
        help='Whether input tsv file has header or not')
parser.add_argument('--output', type=str,
        default='output.tsv', help='output fp')

args = parser.parse_args()

def main():
    bam_utils.index_reference(args.reference_fasta)
    bam_utils.index_bam(args.input_bam)

    rna_editing_annotations = annotation.get_rna_editing_annotations(args.input_bam,
            args.input_annotated_vaf, args.reference_fasta, has_header=args.input_header)

    out_f = open(args.output, 'w')
    f = open(args.input_annotated_vaf)

    annotation_names = [annotation_name 
            for annotation_name, _ in list(rna_editing_annotations.values())[0]]

    if args.input_header:
        line = f.readline().strip()
        out_f.write(line + '\t' + '\t'.join(annotation_names) + '\n')

    for line in f:
        line = line.strip()
        chrom, pos = line.split('\t', 2)[:2]

        position_annotations = rna_editing_annotations[(chrom, pos)]
        out_f.write(line + '\t' + '\t'.join([str(v)
                for _, v in position_annotations]) + '\n')

    f.close()
    out_f.close()

if __name__ == '__main__':
    main()
