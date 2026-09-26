# -*- coding: utf-8 -*-
"""Generation of synthetic data

Module that creates the required synthetic data to run the repo
end-to-end, including the synthetic reference, FASTQs and 'ground truth'.
"""

import random

# Default RNG seed to ensure determinism
SEED = 77
rng = random.Random(SEED)

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
    ref_cpgs = []
    for pos in range(len(ref_seq)):
        if "CG" in ref_seq[pos:pos+2]:
            ref_cpgs += [pos]

    return (ref_seq, ref_cpgs)

def write_fasta(sequence:str, output_path:str, line_char_limit:int = 60):
    """Writes a one-entry FASTA file based on provided sequence

    :param sequence: Sequence string
    :param output_path: Path to save FASTA.txt file
    :param line_char_limit: Maximum length of sequences per line
    """
    line_index = 0
    with open(output_path, "w") as f:
        f.write(">syn\n")
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



if __name__ == "__main__":
    regions = generate_regions()
    reference = generate_reference(regions)
    write_fasta(reference[0], "data/synthetic_reference.txt")
    write_panel_bed(regions, "data/synthetic_panel.bed")