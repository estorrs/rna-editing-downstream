import os
import re
import subprocess
import uuid
from collections import Counter

import bam_utils

def get_position_to_data_dict(input_rna_editing_vaf, has_header=False):
    f = open(input_rna_editing_vaf)
    
    if has_header:
        f.readline()
    
    input_types = ['dna.a', 'dna.t', 'rna.a', 'rna.t']
    fields = ['depth', 'ref_vaf', 'minor_vaf', 'a_vaf', 'c_vaf', 'g_vaf', 't_vaf', 'n_vaf']
    annotations = ['primary_transcript', 'gene', 'strand', 'region', 'non_verbose_region', 'info', 
                   'repeat_name', 'repeat_class', 'repeat_family', 'blat_percent_passing']

    position_to_data_dict = {}

    for line in f:
        pieces = line.strip().split('\t')
        chrom, pos, ref_base = pieces[0], pieces[1], pieces[2]
        
        data_dict = {'ref_base': ref_base}
        data_dict.update({k:{} for k in input_types})
        for i, input_type in enumerate(input_types):
            for j, field in enumerate(fields):
                index = 3 + (len(fields) * i) + j
                val = pieces[index]
                if field == 'depth':
                    val = int(val)
                elif val.isdecimal():
                    val = float(val)
                
                data_dict[input_type][field] = val
    
        for i, annotation in enumerate(annotations):
            index = 3 + (len(fields) * len(input_types)) + i
            val = pieces[index]
            if val.isdecimal():
                val = float(val)
            elif val.isdigit():
                val = int(val)

            data_dict[annotation] = val
        
        position_to_data_dict[(chrom, pos)] = data_dict
    
    return position_to_data_dict

def write_positions_bed(chrom_pos_tups, output_fp):
    f = open(output_fp, 'w')
    for chrom, pos in chrom_pos_tups:
        f.write(f'{chrom}\t{pos}\t{pos}\n')
    f.close()

def write_regions_file(chrom_start_stop_tups, output_fp):
    f = open(output_fp, 'w')
    for chrom, start, stop in chrom_start_stop_tups:
        f.write(f'{chrom}:{start}-{stop}\n')
    f.close()

def get_read_to_reference_seq(read_tups, regions_fp, reference_fasta):
    tool_args = ['samtools', 'faidx', '-r', regions_fp, reference_fasta]
    output = subprocess.check_output(tool_args).decode('utf-8')
    read_to_reference_sequence = bam_utils.get_reads_to_sequences_from_fasta_stream(output)
    
    read_tups_to_reference_seqs = {}
    for chrom, pos, cigar, seq in read_tups:
        end = bam_utils.get_covering_reference_coords(int(pos), cigar, seq)[1]
        reference_seq = read_to_reference_sequence[f'{chrom}:{pos}-{end}']
        read_tups_to_reference_seqs[(chrom, pos, cigar, seq)] = reference_seq
    
    return read_tups_to_reference_seqs

def get_average_base_changes_on_changed_reads(chrom, pos, data_dict, read_collection_reads):
    totals = []
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        
        if not bam_utils.is_match(int(read_start), int(pos), read_cigar, read_seq, reference_seq):
            totals.append(bam_utils.count_mismatches(read_cigar, read_seq, reference_seq))

    if totals:
        return sum(totals) / len(totals)
    return '.'

def get_average_base_changes_on_all_reads(chrom, pos, data_dict, read_collection_reads):
    totals = []
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        
        totals.append(bam_utils.count_mismatches(read_cigar, read_seq, reference_seq))

    if totals:
        return sum(totals) / len(totals)
    return '.'

def get_average_rna_editing_base_changes_on_changed_reads(chrom, pos, data_dict, read_collection_reads):
    totals = []
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        strand = data_dict['strand']
        
        if strand != '+' and strand != '-':
            return '.'
        
        if not bam_utils.is_match(int(read_start), int(pos), read_cigar, read_seq, reference_seq):
            totals.append(bam_utils.count_valid_rna_editing_mismatches(
                    read_cigar, read_seq, reference_seq, strand))

    if totals:
        return sum(totals) / len(totals)
    return '.'

def get_average_rna_editing_base_changes_on_all_reads(chrom, pos, data_dict, read_collection_reads):
    totals = []
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        strand = data_dict['strand']
        
        if strand != '+' and strand != '-':
            return '.'
        
        totals.append(bam_utils.count_valid_rna_editing_mismatches(
                read_cigar, read_seq, reference_seq, strand))

    if totals:
        return sum(totals) / len(totals)
    return '.'

def get_percent_altered_reads_valid_editing_changes(chrom, pos, data_dict, read_collection_reads):
    count = 0
    total = 0
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        strand = data_dict['strand']
        
        if strand != '+' and strand != '-':
            return '.'
        
        if not bam_utils.is_match(int(read_start), int(pos), read_cigar, read_seq, reference_seq):
            if bam_utils.is_valid_rna_editing_site(int(read_start), int(pos),
                                                   read_cigar, read_seq, reference_seq, strand):
                count += 1
            total += 1

    if total:
        return count / total
    return '.'

def get_editing_type(chrom, pos, data_dict, read_collection_reads):
    types = []
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        strand = data_dict['strand']

        if strand != '+' and strand != '-':
            return '.', '.'

        if not bam_utils.is_match(int(read_start), int(pos), read_cigar, read_seq, reference_seq):
            if bam_utils.is_valid_rna_editing_site(int(read_start), int(pos),
                    read_cigar, read_seq, reference_seq, strand):
                types.append(bam_utils.get_rna_editing_type(int(read_start), int(pos),
                        read_cigar, read_seq, reference_seq, strand))
    counter = Counter(types)

    if counter:
        most_common = counter.most_common(1)[0]
        percentage = most_common[1] / len(types)
        return (most_common[0], percentage)
    return '.', '.'


def get_percent_reads_valid_editing_changes(chrom, pos, data_dict, read_collection_reads):
    count = 0
    total = 0
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        strand = data_dict['strand']
        
        if strand != '+' and strand != '-':
            return '.'
        
        count += bam_utils.count_valid_rna_editing_mismatches(
                    read_cigar, read_seq, reference_seq, strand)
        total += bam_utils.count_mismatches(read_cigar, read_seq, reference_seq)

    if total:
        return count / total
    return '.'

def get_is_likely_editing_site(annotation_dict, blat_percent_passing=None):
    if any([annotation_dict['AVG_EDITING_CHANGES_PER_ALTERED_READ'] == '.',
            annotation_dict['AVG_CHANGES_PER_ALTERED_READ'] == '.',
            annotation_dict['ALL_CHANGES_%_VALID_EDITING'] == '.',
            annotation_dict['EDITING_TYPE'] == '.',
            annotation_dict['ALTERED_SITE_%_VALID_EDITING'] == '.']):
        return 'FALSE'
    # make sure percentage of edited reads passes threshold
    read_editing_ratio = annotation_dict['AVG_EDITING_CHANGES_PER_ALTERED_READ'] / annotation_dict['AVG_CHANGES_PER_ALTERED_READ']
    if read_editing_ratio < .8:
        return 'FALSE'
#     # make sure overall editing ratio passes threshold
#     overall_editing_ratio = annotation_dict['ALL_CHANGES_%_VALID_EDITING']
#     if overall_editing_ratio < .8:
#         return 'FALSE'
    # make sure site editing ratio passes threshold
    site_editing_ratio = annotation_dict['ALTERED_SITE_%_VALID_EDITING']
    if site_editing_ratio < .8:
        return 'FALSE'

    editing_type = annotation_dict['EDITING_TYPE']
    percentage = annotation_dict['EDITING_TYPE_%_SUPPORT']
    if percentage < .8:
        return 'FALSE'

    if blat_percent_passing is not None:
        if blat_percent_passing == '.':
            return 'FALSE'
        if float(blat_percent_passing) < .5:
            return 'FALSE'

    return 'TRUE'

def get_annotations(read_collection, position_to_data_dict):
    position_to_annotations = {}
    for (chrom, pos), data_dict in position_to_data_dict.items():
        reads = read_collection.get_reads(chrom, int(pos))
        annotations = {}
        
        annotations['AVG_CHANGES_PER_ALTERED_READ'] = get_average_base_changes_on_changed_reads(chrom,
                pos, data_dict, reads)
        
        annotations['AVG_CHANGES_PER_READ'] = get_average_base_changes_on_all_reads(chrom,
                pos, data_dict, reads)
        
        annotations['AVG_EDITING_CHANGES_PER_ALTERED_READ'] = get_average_rna_editing_base_changes_on_changed_reads(
                chrom, pos, data_dict, reads)
                
        annotations['AVG_EDITING_CHANGES_PER_READ'] = get_average_rna_editing_base_changes_on_all_reads(
          chrom, pos, data_dict, reads)
        
        annotations['ALTERED_SITE_%_VALID_EDITING'] = get_percent_altered_reads_valid_editing_changes(
          chrom, pos, data_dict, reads)

        annotations['ALL_CHANGES_%_VALID_EDITING'] = get_percent_reads_valid_editing_changes(
          chrom, pos, data_dict, reads)

        editing_type, percentage =  get_editing_type(
          chrom, pos, data_dict, reads)
        annotations['EDITING_TYPE'] = editing_type
        annotations['EDITING_TYPE_%_SUPPORT'] = percentage

        annotations['IS_EDITING_SITE'] = get_is_likely_editing_site(annotations,
                blat_percent_passing=data_dict['blat_percent_passing'])

        annotations['IS_EDITING_SITE_NO_BLAT'] = get_is_likely_editing_site(annotations,
                blat_percent_passing=None)

        sorted_annotations = sorted(list(annotations.items()), key=lambda x: x[0])
        print(sorted_annotations)
        
        position_to_annotations[(chrom, pos)] = sorted_annotations

    return position_to_annotations

def get_rna_editing_annotations(input_bam_fp, input_annotated_vaf_fp, reference_fasta_fp,
        has_header=True):
    position_to_data_dict = get_position_to_data_dict(input_annotated_vaf_fp, has_header=has_header)

    chrom_pos_tups = sorted(set(position_to_data_dict.keys()))
    u_id = str(uuid.uuid4())
    temp_positions_fp = f'temp.positions.{u_id}.bed'
    write_positions_bed(chrom_pos_tups, temp_positions_fp)
    read_tups = bam_utils.get_chrom_start_cigar_seq_read_tups(input_bam_fp, temp_positions_fp)

    u_id = str(uuid.uuid4())
    temp_regions_fp = f'temp.regions.{u_id}.txt'
    chrom_start_stop_tups = [(chrom, *bam_utils.get_covering_reference_coords(int(pos), cigar, seq))
            for chrom, pos, cigar, seq in read_tups]
    write_regions_file(chrom_start_stop_tups, temp_regions_fp)
    read_tups_to_reference_seq = get_read_to_reference_seq(read_tups, temp_regions_fp, reference_fasta_fp)

    rc = bam_utils.ReadCollection(chrom_pos_tups)
    for chrom, pos, cigar, seq in read_tups:
        rc.put_read(chrom, pos, cigar, seq,
                reference_sequence=read_tups_to_reference_seq[(chrom, pos, cigar, seq)])

    annotations = get_annotations(rc, position_to_data_dict)

    os.remove(temp_positions_fp)
    os.remove(temp_regions_fp)

    return annotations
