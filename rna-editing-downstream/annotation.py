from collections import Counter
import os
import re
import subprocess
import sys
import uuid

import bam_utils

def get_position_to_data_dict(input_rna_editing_vaf, has_header=False):
    f = open(input_rna_editing_vaf)
    
    if has_header:
        f.readline()
    
    input_types = ['dna', 'rna']
    fields = ['depth', 'ref_vaf', 'minor_vaf', 'a_vaf', 'c_vaf', 'g_vaf', 't_vaf', 'n_vaf']
    annotations = ['primary_transcript', 'gene', 'strand', 'coordinates', 'region', 'non_verbose_region', 'info', 
                   'repeat_name', 'repeat_class', 'repeat_family', 'blat_percent_passing']

    position_to_data_dict = {}

    for line in f:
        pieces = line.strip().split('\t')
        chrom, pos, ref_base, alt_base = pieces[0], pieces[1], pieces[2], pieces[3]
        
        data_dict = {'ref_base': ref_base, 'alt_base': alt_base}
        data_dict.update({k:{} for k in input_types})
        for i, input_type in enumerate(input_types):
            for j, field in enumerate(fields):
                index = 4 + (len(fields) * i) + j
                val = pieces[index]
                if field == 'depth':
                    val = int(val)
                else:
                    try:
                        val = float(val)
                    except ValueError:
                        pass

                data_dict[input_type][field] = val
    
        for i, annotation in enumerate(annotations):
            index = 4 + (len(fields) * len(input_types)) + i
            val = pieces[index]
            import sys
            print(i, index, val, file=sys.stderr)
            if val.isdigit():
                val = int(val)
            else:
                try:
                    val = float(val)
                except ValueError:
                    pass

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
    for chrom, pos, cigar, seq, _, _ in read_tups:
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

def get_editing_level(alt_base, minor_vaf, a_vaf, c_vaf, g_vaf, t_vaf):
    vafs = [v for v in [a_vaf, c_vaf, g_vaf, t_vaf]
            if v <= minor_vaf]
    max_vaf = max(vafs)

    # if sum of lessser vafs is more than 20 percent of the max vaf then return early
    little_vaf_sum = 0
    for vaf in vafs:
        if max_vaf != vaf:
            little_vaf_sum += vaf
    if max_vaf == 0.0 or little_vaf_sum / max_vaf > .2:
        return '.', 0.0

    if 'a' == alt_base.lower():
        return 'A', a_vaf
    if 'c' == alt_base.lower():
        return 'C', c_vaf
    if 'g' == alt_base.lower():
        return 'G', g_vaf
    if 't' == alt_base.lower():
        return 'T', t_vaf

    return '.', 0.0

def get_editing_type_and_percentage(data_dict):
    strand = data_dict['strand']
    ref = data_dict['ref_base']

    d = data_dict['rna']

    alt = data_dict['alt_base']
    alt, vaf = get_editing_level(alt, d['minor_vaf'], d['a_vaf'], d['c_vaf'], d['g_vaf'], d['t_vaf'])

    tup = (alt.lower(), ref.lower())

    if strand != '.' and tup in bam_utils.VALID_RNA_EDITING_CHANGES[strand]:
        return bam_utils.RNA_EDITING_TYPES[strand][tup], vaf
    else:
        return '.', '.'

def get_neighborhood_valid_editing_percentage(chrom, pos, data_dict, read_collection_reads,
        editing_type):
    count = 0
    total = 0
    for read_chrom, read_start, read_cigar, read_seq, read_dict in read_collection_reads:
        reference_seq = read_dict['reference_sequence']
        sequence_quality = read_dict['sequence_quality']
        mapping_quality = read_dict['mapping_quality']
        strand = data_dict['strand']

        if strand != '+' and strand != '-':
            return '.'
        if editing_type == '.':
            return '.'

        if bam_utils.is_valid_rna_editing_site(int(read_start), int(pos),
                read_cigar, read_seq, reference_seq, sequence_quality, mapping_quality, strand,
                editing_type, min_base_quality=25, min_mapping_quality=20):
            count += bam_utils.count_valid_rna_editing_mismatches(
                    read_cigar, read_seq, reference_seq, sequence_quality, strand, editing_type,
                    min_base_quality=25)
            total += bam_utils.count_mismatches(read_cigar, read_seq, reference_seq, sequence_quality,
                    min_base_quality=25)

    if total:
        return count / total
    return '.'

def call_editing_site(editing_type, neighborhood_valid_editing_percentage,
        blat_percent_passing=None):
    if blat_percent_passing is None:
        blat_percent_passing = 1.0

    if '.' in [editing_type, neighborhood_valid_editing_percentage,
            blat_percent_passing]:
        return 'FALSE'

    valid_neighborhood_editing_levels = neighborhood_valid_editing_percentage >= .8
    valid_blat_result = blat_percent_passing >= .5

    if valid_neighborhood_editing_levels and valid_blat_result:
        return 'TRUE'
    return 'FALSE'


def get_annotations(read_collection, position_to_data_dict):
    position_to_annotations = {}
    for (chrom, pos), data_dict in position_to_data_dict.items():
        reads = read_collection.get_reads(chrom, int(pos))
        annotations = {}

        editing_type, percentage = get_editing_type_and_percentage(data_dict)
        annotations['EDITING_TYPE'] = editing_type
        annotations['EDITING_%'] = percentage

        neighborhood_valid_editing_percentage = get_neighborhood_valid_editing_percentage(
                chrom, pos, data_dict, reads, editing_type)

        annotations['NEIGHBORHOOD_VALID_EDITING_%'] = neighborhood_valid_editing_percentage

        annotations['IS_EDITING_SITE'] = call_editing_site(editing_type,
                neighborhood_valid_editing_percentage,
                blat_percent_passing=data_dict['blat_percent_passing'])
        annotations['IS_EDITING_SITE_NO_BLAT'] = call_editing_site(editing_type,
                neighborhood_valid_editing_percentage,
                blat_percent_passing=None)

        sorted_annotations = sorted(list(annotations.items()), key=lambda x: x[0])

        position_to_annotations[(chrom, pos)] = sorted_annotations

    return position_to_annotations

def get_rna_editing_annotations(input_bam_fp, input_annotated_vaf_fp, reference_fasta_fp,
        has_header=True):
    position_to_data_dict = get_position_to_data_dict(input_annotated_vaf_fp, has_header=has_header)

    chrom_pos_tups = sorted(set(position_to_data_dict.keys()))

    # chunk to hopefully not run out of memory
    chunked_chrom_pos_tups = []
    prev = 0
    for i in range(5000, len(chrom_pos_tups) + 5000, 5000):
        chunked_chrom_pos_tups.append(chrom_pos_tups[prev:i])
        prev = i

    annotations = {}
    for chrom_pos_tups_chunk in chunked_chrom_pos_tups:
        if chrom_pos_tups_chunk:

            u_id = str(uuid.uuid4())
            temp_positions_fp = f'temp.positions.{u_id}.bed'
            write_positions_bed(chrom_pos_tups_chunk, temp_positions_fp)
            read_tups = bam_utils.get_chrom_start_cigar_seq_qual_read_tups(input_bam_fp, temp_positions_fp)

            u_id = str(uuid.uuid4())
            temp_regions_fp = f'temp.regions.{u_id}.txt'
            chrom_start_stop_tups = [(chrom, *bam_utils.get_covering_reference_coords(int(pos), cigar, seq))
                    for chrom, pos, cigar, seq, _, _ in read_tups]
            write_regions_file(chrom_start_stop_tups, temp_regions_fp)
            read_tups_to_reference_seq = get_read_to_reference_seq(read_tups, temp_regions_fp, reference_fasta_fp)

            rc = bam_utils.ReadCollection(chrom_pos_tups_chunk)
            for chrom, pos, cigar, seq, qual_seq, map_qual in read_tups:
                rc.put_read(chrom, pos, cigar, seq,
                        reference_sequence=read_tups_to_reference_seq[(chrom, pos, cigar, seq)],
                        sequence_quality=qual_seq, mapping_quality=map_qual)

            annotations.update(get_annotations(rc, position_to_data_dict))

            os.remove(temp_positions_fp)
            os.remove(temp_regions_fp)

    return annotations
