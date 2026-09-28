/*
* Convert reference FASTA to consider bisulfite conversion via Bismark
*/
process bismarkGenomePrep {

    container 'quay.io/biocontainers/bismark:3.1.0--hfa8f182_0'

    input:
    path fasta_dir

    output:
    path "${fasta_dir}/Bisulfite_Genome"

    script:
    """
    bismark_genome_preparation ${fasta_dir}
    """
}

/*
* Pipeline parameters
*/
params {
    fasta_dir: Path = 'data'
}

/*
* Workflow
*/
workflow {
    bismarkGenomePrep(channel.fromPath(params.fasta_dir))
}