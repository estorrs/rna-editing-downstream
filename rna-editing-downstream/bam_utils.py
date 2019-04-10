import os
import re
import subprocess

BOTH_COUNTS = set(['M', 'X', '='])
REFERENCE_COUNTS = set(['N', 'D'])
READ_COUNTS = set(['I', 'S', 'H', 'P'])
IDENTIFIER_SPLIT_REGEX = re.compile(r'M|X|=|N|D|I|S|H|P')
COUNT_SPLIT_REGEX = re.compile(r'[0-9]+')

VALID_RNA_EDITING_CHANGES = {
        '-': {('c', 't'), ('a', 'g')},
        '+': {('g', 'a'), ('t', 'c')}
        }

RNA_EDITING_TYPES = {
        '-': {('c', 't'): 'A->I', ('a', 'g'): 'C->U'},
        '+': {('g', 'a'): 'A->I', ('t', 'c'): 'C->U'}
        }

R_RNA_EDITING_TYPES = {
        '-': {'A->I': ('c', 't'), 'C->U': ('a', 'g')},
        '+': {'A->I': ('g', 'a'), 'C->U': ('t', 'c')}
        }

def index_reference(reference_fasta_fp):
    """index the given reference if it doesnt exist"""
    if not os.path.isfile(reference_fasta_fp + '.fai'):
        tool_args = ['samtools', 'faidx', reference_fasta_fp]
        print(subprocess.check_output(tool_args).decode('utf-8'))

def index_bam(bam_fp):
    """index the given bam."""
    if not os.path.isfile(bam_fp + '.bai'):
        tool_args = ['samtools', 'index', bam_fp]
        print(subprocess.check_output(tool_args).decode('utf-8'))

def filter_bam_by_positions(bam_fp, positions_fp, output_fp, threads=1):
    """run bam filter step"""
    tool_args = ['samtools', 'view', '-h',
            '-@', str(threads),
            '-L', positions_fp,
            '-o', output_fp,
            bam_fp]

    print(subprocess.check_output(tool_args).decode('utf-8'))

def get_covering_reference_coords(start, cigar, seq):
    ref_end = int(start)
    counts = re.split(IDENTIFIER_SPLIT_REGEX, cigar)[:-1]
    identifiers = re.split(COUNT_SPLIT_REGEX, cigar)[1:]

    for count, identifier in zip(counts, identifiers):
        ref_end += int(count)

    return int(start), ref_end - 1

def count_mismatches(cigar, read_seq, reference_seq, base_quality_seq, min_base_quality=20):
    read_counter = 0
    ref_counter = 0

    counts = [int(c) for c in re.split(IDENTIFIER_SPLIT_REGEX, cigar)[:-1]]
    identifiers = re.split(COUNT_SPLIT_REGEX, cigar)[1:]

    mismatches = 0
    for count, identifier in zip(counts, identifiers):
        if identifier in BOTH_COUNTS:
            read_nucleotides = read_seq[read_counter:read_counter + count].lower()
            base_qualities = base_quality_seq[read_counter:read_counter + count]
            ref_nucleotides = reference_seq[ref_counter:ref_counter + count].lower()
            for read_base, base_quality, ref_base in zip(read_nucleotides, base_qualities,
                    ref_nucleotides):
                converted_base_quality = ord(base_quality) - 33
                if read_base != ref_base and converted_base_quality >= min_base_quality:
                    mismatches += 1

            read_counter += count
            ref_counter += count

        elif identifier in REFERENCE_COUNTS:
            ref_counter += count
        elif identifier in READ_COUNTS:
            read_counter += count

    return mismatches

# this is bad repition, incorporate into other function
def count_valid_rna_editing_mismatches(cigar, read_seq, reference_seq, base_quality_seq, strand,
        editing_type=None, min_base_quality=20):
    read_counter = 0
    ref_counter = 0

    counts = [int(c) for c in re.split(IDENTIFIER_SPLIT_REGEX, cigar)[:-1]]
    identifiers = re.split(COUNT_SPLIT_REGEX, cigar)[1:]

    mismatches = 0
    for count, identifier in zip(counts, identifiers):
        if identifier in BOTH_COUNTS:
            read_nucleotides = read_seq[read_counter:read_counter + count].lower()
            base_qualities = base_quality_seq[read_counter:read_counter + count]
            ref_nucleotides = reference_seq[ref_counter:ref_counter + count].lower()

            for read_base, base_quality, ref_base in zip(read_nucleotides, base_qualities, ref_nucleotides):
                if editing_type is None:
                    is_editing_change = (read_base, ref_base) in VALID_RNA_EDITING_CHANGES[strand]
                else:
                    expected_tup = R_RNA_EDITING_TYPES[strand][editing_type]
                    is_editing_change = (read_base, ref_base) == expected_tup

                converted_base_quality = ord(base_quality) - 33
                passes_base_quality = converted_base_quality >= min_base_quality

                if read_base != ref_base and is_editing_change and passes_base_quality:
                    mismatches += 1

            read_counter += count
            ref_counter += count

        elif identifier in REFERENCE_COUNTS:
            ref_counter += count
        elif identifier in READ_COUNTS:
            read_counter += count

    return mismatches

def is_valid_rna_editing_site(start, target_pos, cigar, read_seq, reference_seq, strand):
    read_counter = 0
    ref_counter = 0

    counts = [int(c) for c in re.split(IDENTIFIER_SPLIT_REGEX, cigar)[:-1]]
    identifiers = re.split(COUNT_SPLIT_REGEX, cigar)[1:]

    mismatches = 0
    for count, identifier in zip(counts, identifiers):
        if identifier in BOTH_COUNTS:
            for i in range(count):
                if start + i + ref_counter == target_pos:
                    if (read_seq[read_counter + i].lower(),
                            reference_seq[ref_counter + i].lower()) in VALID_RNA_EDITING_CHANGES[strand]:
                        return True
                    else:
                        return False
            read_counter += count
            ref_counter += count

        elif identifier in REFERENCE_COUNTS:
            ref_counter += count

        elif identifier in READ_COUNTS:
            read_counter += count

    return None

def get_rna_editing_type(start, target_pos, cigar, read_seq, reference_seq, strand):
    read_counter = 0
    ref_counter = 0

    counts = [int(c) for c in re.split(IDENTIFIER_SPLIT_REGEX, cigar)[:-1]]
    identifiers = re.split(COUNT_SPLIT_REGEX, cigar)[1:]

    mismatches = 0
    for count, identifier in zip(counts, identifiers):
        if identifier in BOTH_COUNTS:
            for i in range(count):
                if start + i + ref_counter == target_pos:
                    read_base = read_seq[read_counter + i].lower()
                    ref_base = reference_seq[ref_counter + i].lower()
                    if (read_base, ref_base) in VALID_RNA_EDITING_CHANGES[strand]:
                        return RNA_EDITING_TYPES[strand][(read_base, ref_base)]
                    else:
                        return '.'

            read_counter += count
            ref_counter += count

        elif identifier in REFERENCE_COUNTS:
            ref_counter += count

        elif identifier in READ_COUNTS:
            read_counter += count

    return '.'

def is_match(start, target_pos, cigar, read_seq, reference_seq):
    read_counter = 0
    ref_counter = 0

    counts = [int(c) for c in re.split(IDENTIFIER_SPLIT_REGEX, cigar)[:-1]]
    identifiers = re.split(COUNT_SPLIT_REGEX, cigar)[1:]

    mismatches = 0
    for count, identifier in zip(counts, identifiers):
        if identifier in BOTH_COUNTS:
            for i in range(count):
                if start + i + ref_counter == target_pos:
                    if read_seq[read_counter + i].lower() == reference_seq[ref_counter + i].lower():
                        return True
                    else:
                        return False
            read_counter += count
            ref_counter += count

        elif identifier in REFERENCE_COUNTS:
            ref_counter += count

        elif identifier in READ_COUNTS:
            read_counter += count

    return None
    
class ReadCollection(object):
    def __init__(self, position_tups):
        """
        position tups - [(chrom, pos), ...]
        """
        position_tups = [(c, int(p)) for c, p in position_tups]
        self.positions = sorted(position_tups, key=lambda x: (x[0], x[1]))

        self.chrom_to_positions = {c:[] for c, _ in position_tups}
        for chrom, pos in position_tups:
            self.chrom_to_positions[chrom].append((chrom, pos))
        for chrom, ps in self.chrom_to_positions.items():
            self.chrom_to_positions[chrom] = sorted(ps, key=lambda x: (x[0], x[1]))

        self.position_to_reads = {p:[] for p in position_tups}

    def put_read(self, chrom, start, cigar, sequence, **kwargs):
        start, end = get_covering_reference_coords(int(start), cigar, sequence)
        start, end = int(start), int(end)

        for chrom, pos in self.chrom_to_positions[chrom]:
            if pos > start and pos < end:
                self.position_to_reads[(chrom, pos)].append((chrom, start, cigar, sequence, kwargs))

            if pos > end:
                break

    def get_reads(self, chrom, pos):
        return self.position_to_reads.get((chrom, pos), [])


def get_reads_to_sequences_from_fasta(input_fasta_fp):
    f = open(input_fasta_fp)

    active_seq_id = None
    active_seq = ''
    reads_to_sequences = {}
    for line in f:
        line = line.strip()
        if line[0] == '>':
            if active_seq_id is not None:
                reads_to_sequences[active_seq_id] = active_seq
            active_seq_id = line[1:]
        else:
            active_seq += line
    f.close()

    return reads_to_sequences

def get_reads_to_sequences_from_fasta_stream(input_fasta_stream):
    lines = input_fasta_stream.split('\n')

    active_seq_id = None
    active_seq = ''
    reads_to_sequences = {}
    for line in lines:
        if line:
            line = line.strip()
            if line[0] == '>':
                if active_seq_id is not None:
                    reads_to_sequences[active_seq_id] = active_seq
                active_seq_id = line[1:]
                active_seq = ''
            else:
                active_seq += line
    # grab last one
    reads_to_sequences[active_seq_id] = active_seq

    return reads_to_sequences

def get_chrom_start_cigar_seq_qual_read_tups(input_bam_fp, positions_fp, max_depth=500):
    tool_args = ['samtools', 'view',
            '-L', positions_fp,
             input_bam_fp]
    ps_1 = subprocess.Popen(tool_args, stdout=subprocess.PIPE)
    output = subprocess.check_output(('cut', '-f', '3,4,6,10,11'), stdin=ps_1.stdout).decode('utf-8')
    ps_1.wait()

    read_tups = []
    position_to_depth = {}
    for line in output.split('\n'):
        if line:
            pieces = line.split('\t')
            chrom, pos, cigar, seq, qual_seq = pieces[0], pieces[1], pieces[2], pieces[3], pieces[4]

            if (chrom, pos) in position_to_depth:
                position_to_depth[(chrom, pos)] += 1
            else:
                position_to_depth[(chrom, pos)] = 1

            if position_to_depth[(chrom, pos)] < max_depth:
                read_tups.append((chrom, pos, cigar, seq, qual_seq))

    return read_tups

def write_position_fasta(input_bam_fp, positions_fp, output_fasta_fp, max_depth=500):
    """Writes a fasta with the given positions and bam.

    Will also return a dict mapping reads to their sequence"""
    tool_args = ['samtools', 'view',
            '-L', positions_fp,
             input_bam_fp]
    ps_1 = subprocess.Popen(tool_args, stdout=subprocess.PIPE)
    output = subprocess.check_output(('cut', '-f', '3,4,10'), stdin=ps_1.stdout).decode('utf-8')
    ps_1.wait()
    
    read_tups = []
    for line in output.split('\n'):
        # line is chrom\tstart\tsequence\n
        pieces = line.split('\t')
        if len(pieces) == 3:
            read_tups.append((pieces[0], int(pieces[1]), pieces[2]))

    # grab positions from file
    f = open(positions_fp)
    positions = []
    for line in f:
        # line is chrom\tpos....\n
        pieces = line.strip().split('\t')
        positions.append((pieces[0], int(pieces[1])))
    f.close()

    # create read collection
    rc = ReadCollection(positions)
    for chrom, start, seq in read_tups:
        rc.put_read(chrom, start, seq)
    
    f = open(output_fasta_fp, 'w')
    reads_to_data = {}
    for chrom, pos in positions:
        reads = rc.get_reads(chrom, pos)
        for i, (_, _, seq, base) in enumerate(reads):
            f.write(f'>{chrom}:{pos}|{i}\n')
            f.write(seq + '\n')

            reads_to_data[f'{chrom}:{pos}|{i}'] = {
                    'sequence': seq,
                    'sequence_base': base
                    }

            if i >= max_depth:
                break
    f.close()

    return reads_to_data
