import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'

INPUT_BAM = os.path.join(TEST_DATA_DIR,
        'C3N-00435.chr1.filtered.bam')
INPUT_ANNOTATED_VAF = os.path.join(TEST_DATA_DIR,
        'C3N-00435.chr1.tsv')
INPUT_ANNOTATED_VAF = os.path.join(TEST_DATA_DIR,
        'new.test.tsv')
REFERENCE_FASTA = os.path.join(TEST_DATA_DIR,
        'hg38.fa')

def test_rna_downstream_annotation():
    tool_args = ['python', 'rna-editing-downstream/rna_editing_downstream.py',
            '--input-header',
            '--output', 'ouput.tsv',
            INPUT_BAM, INPUT_ANNOTATED_VAF, REFERENCE_FASTA]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert True
