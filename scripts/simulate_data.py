# -*- coding: utf-8 -*-
"""Generation of synthetic data

Module that creates the required synthetic data to run the repo
end-to-end, including the synthetic reference, FASTQs and 'ground truth'.
"""

import random
import os

# Default RNG seed to ensure determinism
SEED = 77
rng = random.Random(SEED)

def _find_cpg_pos(ref_seq):
    ref_cpgs = []
    for pos in range(len(ref_seq)):
        if "CG" in ref_seq[pos:pos+2]:
            ref_cpgs += [pos]
    return ref_cpgs

def _make_unique_regions(n_regions: int, region_length: int):
    """Return n_region unique regions to be specified as desired

    :param n_regions: number of unique regions to generate
    :param region_length: length of each unique region
    :returns: list of length len with DNA in FASTA format
    """

    # setup 
    BASES = ["A", "C", "G", "T"]
    twomers = [first + second for first in BASES for second in BASES]
    seq_weights = [10 if "CG" == pair else 1 for pair in twomers]
    regions_list = []

    # generate regions
    region_ind = 0
    while region_ind < n_regions:
        _n_twomers = rng.choices(twomers, weights=seq_weights, k=region_length//2)
        _seq = "".join(_n_twomers)
        if _seq not in regions_list:
            regions_list += [_seq]
            region_ind += 1

    return regions_list

def _assign_methylation(cpg_pos: list, frac: float = 0.2):
    """Assigns methlyated/unmethylated (1/0) status randomly

    :param cpg_pos: list of CpG positions relative to a reference genome
    :param frac: fraction of CpGs to be assigned methylation (i.e. 1)
    :returns: list of methylation calls (0/1) with length len(cpg_pos)
    """
    return [rng.binomialvariate(1, frac) for pos in cpg_pos]

def _generate_fragment_from_region(region_seq, min_len, max_len):
    region_end = len(region_seq)
    frag_len = rng.sample(range(min_len, max_len), k=1)[0]
    start_bp = [start for start in range(1, region_end - frag_len)]
    sub_start = rng.sample(start_bp, k=1)[0]
    sub_end = sub_start + frag_len
    sub_seq = region_seq[sub_start:sub_end+1]

    return sub_seq

def _simulate_sample(regions, label):
    """Simulates case or control samples. Cases are simulated to have
    an approxate 5:1 ratio of on-target vs. off-target reads. Controls are
    simulated to have a 1:1 ratio of on-target vs. off-target reads. Both 
    cases and controls are simulated to have 2 of each synthetic read.

    :param regions: Tuple of regions in (on-target, off-target, synthetic) format
    :param label: either 'case' or 'control', defines sampling behvaior
    :returns: a Tuple with on-target, off-target, and synthetic sequences
    """
    # minimum and maximum fragment length from a region
    MIN_LENGTH = 100
    MAX_LENGTH = 170

    # ensure condition is either 'case' or 'control'
    if label.upper() not in ["CASE", "CONTROL"]:
        raise ValueError(f"'{label}' must be either 'case' or 'control'")

    # generate reads based on specified convention
    if label.upper() == "CASE":
        # on-target regions
        ontarget_regions = rng.choices(regions[0], k = 10)
        ontarget_fragments = [
            _generate_fragment_from_region(region, MIN_LENGTH, MAX_LENGTH) 
            for region in ontarget_regions
        ]
        # off-target regions
        offtarget_regions = rng.choices(regions[0], k = 2)
        offtarget_fragments = [
            _generate_fragment_from_region(region, MIN_LENGTH, MAX_LENGTH) 
            for region in offtarget_regions
        ]
    if label.upper() == "CONTROL":
        # on-target regions
        ontarget_regions = rng.choices(regions[0], k = 3)
        ontarget_fragments = [
            _generate_fragment_from_region(region, MIN_LENGTH, MAX_LENGTH) 
            for region in ontarget_regions
        ]
        # off-target regions
        offtarget_regions = rng.choices(regions[0], k = 3)
        offtarget_fragments = [
            _generate_fragment_from_region(region, MIN_LENGTH, MAX_LENGTH) 
            for region in offtarget_regions
        ]
    # synthetic fragments
    synthetic_fragments = [read for read in regions[2]*2]

    return (ontarget_fragments, offtarget_fragments, synthetic_fragments)

def _bisulfite_convert(seq_list, mfrac=0.8):
    """Performs in-silico bisulfite conversion of a list of sequences. 
    Unmethylated CpGs are converted from C>T

    :param seq_list: List[str]
    :param mfrac: Fraction of CpGs to remain non-converted
    :returns: List[str] after conversion
    """
    seq_cpg_pos = [_find_cpg_pos(seq) for seq in seq_list]
    seq_meth_call = [_assign_methylation(pos, mfrac) for 
                          pos in seq_cpg_pos]
    seq_converted = []
    for i, seq in enumerate(seq_list):
        seq_l = list(seq)
        for j, pos in enumerate(seq_cpg_pos[i]):
            if seq_meth_call[i][j] == 0:
                seq_l[pos] = 'T'
        str_seq_converted = "".join(seq_l)
        seq_converted += [str_seq_converted]

    return seq_converted

def _generate_antisense_sequence(seq:str):
    """Generates anti-sense from sense-strand fragment

    :param seq: sense-strand fragment in the 5'->3' direction
    :returns: anti-sense fragment in 5'->3' direction
    """
    bp_mapping = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    table = str.maketrans(bp_mapping)
    anti_seq = seq.translate(table)

    return anti_seq[::-1]

def _add_umi_adapters(seq: str, adapter_file: str, umi_format: str = 'NNT'):
    """Adds UMIs and library adapter sequences to 5' and 3' of sequence

    :param seq: DNA fragment in 5'->3' direction before 'read generation'
    :param adapter_file: path to .txt file that specifies the adapter used
    :param umi_format: number of random bases ('N') and/or spaces ('A/G/C/T')
    :returns: seq with appened UMIs and adapter sequences
    """
    BASES = ['A', 'C', 'G', 'T']
    # Build required compartments
    with open(adapter_file, "r", encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                continue
            adapter_seq = line.strip()
    r1_seq = seq
    r2_seq = _generate_antisense_sequence(seq)
    r1_umi_a, r1_umi_b, r2_umi_a, r2_umi_b = [], [], [], []
    for char in umi_format:
        if char == 'N':
            r1_umi_a += rng.choices(BASES, k=1)
            r1_umi_b += rng.choices(BASES, k=1)
            r2_umi_a += rng.choices(BASES, k=1)
            r2_umi_b += rng.choices(BASES, k=1)
        else:
            r1_umi_a += [char]
            r1_umi_b += [char]
            r2_umi_a += [char]
            r2_umi_b += [char]
    r1_umi_a_str = "".join(r1_umi_a)
    r1_umi_b_str = "".join(r1_umi_b)
    r2_umi_a_str = "".join(r2_umi_a)
    r2_umi_b_str = "".join(r2_umi_b)
    # Construct R1 and R2 after sequencing
    r1_truth = r1_umi_a_str + \
        r1_seq + \
        _generate_antisense_sequence(r2_umi_b_str) + \
        adapter_seq
    r2_truth = r2_umi_a_str + \
        r2_seq + \
        _generate_antisense_sequence(r1_umi_b_str) + \
        adapter_seq

    return r1_truth, r2_truth
    
def _generate_read_qual(read:str, frac_pass:float = 0.9):
    """Generates read quality scores, with >= frac_pass
    being above 20

    :param read: read string, maps quality scores 1-to-1
    :param frac_pass: fraction of base-pairs to have a qual > 20
    :returns: quality score string
    """
    n_bases = len(read)
    # divide weights evenly between 20-40 
    pass_weight = frac_pass/20
    fail_weight = (1-frac_pass)/21 # includes 0
    weights = [fail_weight]*21 + [pass_weight]*20
    qual_list = [chr(code) for code in range(33,74)]
    seq_quals = rng.choices(qual_list, weights=weights, k=n_bases)
    seq_quals = "".join(seq_quals)

    return seq_quals

def generate_regions(
        n_synthetic: int = 2,
        n_ontarget: int = 2,
        n_offtarget: int = 2
    ):
    """Generates synthetic, on-target, and off-target regions

    :param n_synthetic: number of generated synthetic molecules
    :param n_ontarget: number of generated on-target molecules
    :param n_offtarget: number of generated off-target molecules
    :returns: Tuple[List[synthetic], List[ontarget], List[offtarget]]
    """

    # Create synthetic sequences
    ontarget_regions = _make_unique_regions(n_ontarget, 300)
    offtarget_regions = _make_unique_regions(n_offtarget, 300)
    synthetic_regions = _make_unique_regions(n_synthetic, 150)

    return (ontarget_regions, offtarget_regions, synthetic_regions)

def generate_reference(regions: tuple):
    """Generates reference sequence from generate_regions() output

    :param regions: Tuple output from generate_regions()
    :returns: Tuple of synthetic sequence and CpG position
    """
    # Collate regions from ontarget -> offtarget -> synthetic 
    # to simulate a reference genome
    ref_seq = (
        ""
        .join(regions[0] + regions[1] + regions[2])
    )

    # Get reference CpG sites
    ref_cpgs = _find_cpg_pos(ref_seq)

    return (ref_seq, ref_cpgs)

def write_fasta(
        sequence:str, 
        output_path:str, 
        header:str = 'syn', 
        line_char_limit:int = 60):
    """Writes a one-entry FASTA file based on provided sequence

    :param sequence: Sequence string
    :param output_path: Path to save FASTA.txt file
    :param line_char_limit: Maximum length of sequences per line
    """
    line_index = 0
    with open(output_path, "w") as f:
        f.write(f">{header}\n")
        while line_index < len(sequence):
            f.write(sequence[line_index:(line_index+line_char_limit)] + "\n")
            line_index += line_char_limit

def write_panel_bed(regions:tuple, output_path:str):
    """Writes on-target panel regions in BED format

    :param regions: Tuple output from generate_regions() output
    :param output_path: Path to save panel.bed file
    """
    with open(output_path, "w") as f:
        i = 0
        for region in regions[0]:
            ref = "syn"
            start = i
            end = i + len(region)
            f.write(f"{ref}\t{start}\t{end}\n")
            i = end

def write_truth_tsv(cpg_pos:list, output_path:str):
    """Writes ground truth methylation status (0: unmethylated, 1: methylated)
    to .tsv file

    :param cpg_pos: CpG position from sequence
    :param output_path: Output file path for .tsv file
    """
    with open(output_path, "w") as f:
        f.write(cpg_pos)

def create_sample_fastq(
        regions:tuple, 
        label:str, 
        adapter_file:str,
        conversion_efficiency:float,
        **kwargs
    ):
    """Outputs R1+R2 reads for a specified label (case/control) based on the
    provided on/off-target and synthetic regions. Includes UMIs and library 
    adapters as well as bisulfite converison

    :param regions: output from generate_regions()
    :param label: must be specified as 'case' or 'control'
    :param adapter_file: path of FASTA file with adapter sequence to trim
    :param conversion_efficiency: percentage of methyalted fragments actually methylated
    :param umi_format: structure of UMI, e.g. NNT
    :type umi_format: str, optional
    :returns: List[Tuple(R1, R2) ...]
    """
    # Generate on/off-target and synthetic fragments based on label
    fragments = _simulate_sample(regions, label)
    if label.upper() == 'CASE':
        c_mfrac = 0.8 * conversion_efficiency
    else:
        c_mfrac = 0.2 * conversion_efficiency    
    fragments_bs = [_bisulfite_convert(frag, mfrac=c_mfrac) for 
                    frag in fragments[0:2]] + \
                        [_bisulfite_convert(fragments[2], mfrac=1-conversion_efficiency)]
    fragments_bs_lib = [
        _add_umi_adapters(entry, adapter_file) for 
            sublist in fragments_bs for entry in sublist
    ]

    return fragments_bs_lib


def write_fastq(sample_reads:list, sample_name:str, out_dir:str, read_length:int = 150):
    '''Writes FASTQ file and simulates read_length

    :param sample_reads: output from _create_sample_fastq
    :param sample_name: sample name to include in file name
    :param out_dir: base directory to write fastq/sample.fastq
    :param read_length: simulated read length of full length fragments
    '''
    # write fastq directory at out_dir if not present
    fastq_dir = f'{out_dir}/fastq'
    if os.path.exists(fastq_dir) == False:
        os.mkdir(fastq_dir)
    sample_prefix = f'{fastq_dir}/{sample_name}_L001'

    # generate quality scores in ASCII format for each sequence
    r1_reads, r2_reads = (
        [r1[0:read_length] for r1, _ in sample_reads], 
        [r2[0:read_length] for _, r2 in sample_reads]
    )
    r1_qual = [_generate_read_qual(read) for read in r1_reads]
    r2_qual = [_generate_read_qual(read) for read in r2_reads]

    # create simple read IDs to append 1/2 later
    read_ids = [(f'@READ:{i}/1', f'@READ:{i}/2') for i, _ in enumerate(sample_reads)]

    # collate read IDs, sequence, and quality score and write to FASTQ
    file_name_out = fastq_dir + f'/{sample_name}'

    # read 1
    with open(f'{file_name_out}_R1_001.fastq', 'w') as file:
        i = 0
        while i < len(read_ids):
            file.write(f'{read_ids[i][0]}\n{r1_reads[i]}\n+\n{r1_qual[i]}\n')
            i += 1
        file.close()
    # read 2
    with open(f'{file_name_out}_R2_001.fastq', 'w') as file:
        i = 0
        while i < len(read_ids):
            file.write(f'{read_ids[i][1]}\n{r2_reads[i]}\n+\n{r2_qual[i]}\n')
            i += 1
        file.close()

def write_cohort_fastq_files(
        regions:tuple, 
        sample_sheet:str, 
        adapter_file:str, 
        out_dir:str
    ):
    '''Writes FASTQ files for multiple samples specified by sample_sheet.csv

    :param regions: output from generate_regions() function
    :param sample_sheet: CSV file with columns 'sample_id', 'condition', and 'conversion_efficiency'
    :param adatper_file: .txt file with expected adapter sequence in read
    :param out_dir: base directory for fastq files to be written (<out_dir>/fastq/*.fastq)
    '''
    # read in sample_sheet.csv
    cohort_list = []
    with open(sample_sheet, 'r') as file:
        for line in file:
            row = line.split(',')
            row = [item.strip() for item in row]
            if row == ['sample_id','condition','conversion_efficiency']:
                sample_col = row.index('sample_id')
                cond_col = row.index('condition')
                conv_col = row.index('conversion_efficiency')
                continue
            sample_id, condition, conversion_efficiency = \
                str(row[sample_col]), str(row[cond_col]), float(row[conv_col])
            cohort_list += [(sample_id, condition, conversion_efficiency)]
    file.close()
    # write FASTQ files for each sample in the cohort
    for sample in cohort_list:
        sample_fastq = create_sample_fastq(
            regions, 
            sample[1], 
            adapter_file,
            sample[2]
        )
        print(f'Writing FASTQ for sample {sample[0]} at {out_dir}\n')
        write_fastq(sample_fastq, sample[0], out_dir)

if __name__ == "__main__":
    regions = generate_regions()
    reference = generate_reference(regions)
    write_fasta(reference[0], "data/synthetic_reference.txt")
    write_fasta(
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "data/adapter.txt",
        "Illumina_Universal_Adapter"
    )
    write_panel_bed(regions, "data/synthetic_panel.bed")
    write_cohort_fastq_files(
        regions, 'data/sample_sheet.csv', 'data/adapter.txt', 'data'
    )