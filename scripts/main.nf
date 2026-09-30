/*
* Convert reference FASTA to consider bisulfite conversion via Bismark
*/
process bismarkGenomePrep {

    container 'quay.io/biocontainers/bismark:3.1.0--hfa8f182_0'

    input:
    path fasta_dir

    output:
    path fasta_dir

    script:
    """
    bismark_genome_preparation ${fasta_dir}
    """
}

/*
* Extract UMIs from FASTQ files and append them to the read names. The
* /1 and /2 mate suffixes are dropped so umi_tools can match read pairs.
* sed edits the header and n;n;n passes the other three lines of each record
* through, so a quality line starting with '@' is never touched
*/
process extractUMI {

    container 'quay.io/biocontainers/umi_tools:1.1.6--py312h0fa9677_0'
    tag "${sample}"

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1.umi.fastq"), path("${sample}_R2.umi.fastq")

    script:
    """
    sed 's#/[12]\$##;n;n;n' ${r1} > ${sample}_R1.renamed.fastq
    sed 's#/[12]\$##;n;n;n' ${r2} > ${sample}_R2.renamed.fastq

    umi_tools extract \\
        --bc-pattern=${params.umi_pattern} \\
        --bc-pattern2=${params.umi_pattern} \\
        --stdin=${sample}_R1.renamed.fastq \\
        --read2-in=${sample}_R2.renamed.fastq \\
        --stdout=${sample}_R1.umi.fastq \\
        --read2-out=${sample}_R2.umi.fastq \\
        --log=${sample}.umi_extract.log
    """
}

/*
* Trim library adapters from the 3' end of both reads. The 3' UMI sits
* directly upstream of the adapter, so it is trimmed along with it
*/
process trimAdapters {

    container 'quay.io/biocontainers/cutadapt:5.2--py312hfabe715_2'
    tag "${sample}"

    input:
    tuple val(sample), path(r1), path(r2)
    path adapter

    output:
    tuple val(sample), path("${sample}_R1.trim.fastq"), path("${sample}_R2.trim.fastq")

    script:
    """
    ADAPTER=\$(grep -v '>' ${adapter})

    cutadapt \\
        -a "NNN\${ADAPTER}" \\
        -A "NNN\${ADAPTER}" \\
        -q 20 \\
        -m 30 \\
        -o ${sample}_R1.trim.fastq \\
        -p ${sample}_R2.trim.fastq \\
        ${r1} ${r2} > ${sample}.cutadapt.log
    """
}

/*
* Align trimmed read pairs to the bisulfite converted reference
*/
process bismarkAlign {

    container 'quay.io/biocontainers/bismark:3.1.0--hfa8f182_0'
    tag "${sample}"

    input:
    tuple val(sample), path(r1), path(r2)
    path genome_dir

    output:
    tuple val(sample), path("${sample}_pe.bam")
    path "${sample}_PE_report.txt"

    script:
    """
    bismark \\
        --genome ${genome_dir} \\
        --basename ${sample} \\
        -1 ${r1} \\
        -2 ${r2}
    """
}

/*
* Coordinate sort and index alignments, required for UMI deduplication
*/
process sortBam {

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    tag "${sample}"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

    script:
    """
    samtools sort -o ${sample}.sorted.bam ${bam}
    samtools index ${sample}.sorted.bam
    """
}

/*
* Collapse PCR duplicates using alignment position and UMI
*/
process dedupUMI {

    container 'quay.io/biocontainers/umi_tools:1.1.6--py312h0fa9677_0'
    tag "${sample}"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.dedup.bam")

    script:
    """
    umi_tools dedup \\
        --paired \\
        --stdin=${bam} \\
        --stdout=${sample}.dedup.bam \\
        --output-stats=${sample}.dedup \\
        --log=${sample}.dedup.log
    """
}

/*
* Name sort alignments so mates are adjacent for methylation calling
*/
process nameSortBam {

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    tag "${sample}"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.nsort.bam")

    script:
    """
    samtools sort -n -o ${sample}.nsort.bam ${bam}
    """
}

/*
* Call methylation for every cytosine from the deduplicated alignments
*/
process methylationExtract {

    container 'quay.io/biocontainers/bismark:3.1.0--hfa8f182_0'
    tag "${sample}"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.nsort.bismark.cov.gz")

    script:
    """
    bismark_methylation_extractor \\
        -p \\
        --comprehensive \\
        --bedGraph \\
        ${bam}
    """
}

/*
* Estimate bisulfite conversion efficiency as 100 - % CpG methylation in
* the synthetic spike-in, which is unmethylated by design
*/
process conversionQC {

    container 'quay.io/biocontainers/bismark:3.1.0--hfa8f182_0'
    tag "${sample}"

    input:
    tuple val(sample), path(cov)
    path spikein

    output:
    tuple val(sample), path(cov), stdout

    script:
    """
    zcat ${cov} | awk '
        NR == FNR { chr = \$1; from = \$2; to = \$3; next }
        \$1 == chr && \$2 > from && \$2 <= to { meth += \$5; total += \$5 + \$6 }
        END { printf "%.2f", (total > 0) ? 100 - 100 * meth / total : 0 }
    ' ${spikein} -
    """
}

/*
* Summarize on-target rate and CpG methylation for a single sample
*/
process sampleSummary {

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    tag "${sample}"

    input:
    tuple val(sample), path(bam), path(cov), val(conversion)
    path panel

    output:
    path "${sample}.summary.tsv"

    script:
    """
    TOTAL=\$(samtools view -c -f 64 ${bam})
    ON_TARGET=\$(samtools view -c -f 64 -L ${panel} ${bam})

    printf "sample\\tconversion_pct\\tdedup_pairs\\ton_target_pairs\\ton_target_frac\\ton_target_cpgs\\ton_target_meth_pct\\n" > ${sample}.summary.tsv
    zcat ${cov} | awk -v s=${sample} -v c=${conversion} -v t=\$TOTAL -v o=\$ON_TARGET '
        NR == FNR { n++; chr[n] = \$1; from[n] = \$2; to[n] = \$3; next }
        {
            for (i = 1; i <= n; i++) {
                if (\$1 == chr[i] && \$2 > from[i] && \$2 <= to[i]) {
                    meth += \$5; total += \$5 + \$6; cpgs++
                }
            }
        }
        END {
            frac = (t > 0) ? o / t : 0
            pct = (total > 0) ? 100 * meth / total : 0
            printf "%s\\t%s\\t%d\\t%d\\t%.3f\\t%d\\t%.2f\\n", s, c, t, o, frac, cpgs, pct
        }' ${panel} - >> ${sample}.summary.tsv
    """
}

/*
* Pipeline parameters
*/
params {
    samplesheet: Path = 'data/nextflow_sample_sheet.txt'
    fasta_dir: Path = 'data'
    adapter: Path = 'data/adapter.txt'
    panel: Path = 'data/synthetic_panel.bed'
    spikein: Path = 'data/synthetic_spikein.bed'
    umi_pattern: String = 'NNX'
    min_conversion: Float = 95.0
}

/*
* Workflow
*/
workflow {

    main:
    // one (sample, R1, R2) tuple per row of the sample sheet
    reads = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

    genome = bismarkGenomePrep(channel.fromPath(params.fasta_dir)).first()

    // preprocessing
    umi_reads = extractUMI(reads)
    trimmed_reads = trimAdapters(umi_reads, file(params.adapter))

    // alignment and deduplication
    aligned = bismarkAlign(trimmed_reads, genome)
    dedup = dedupUMI(sortBam(aligned[0]))

    // methylation calling, samples below the conversion threshold are dropped
    calls = methylationExtract(nameSortBam(dedup))
    qc = conversionQC(calls, file(params.spikein))
    qc_pass = qc.filter { _sample, _cov, conversion -> conversion.toFloat() >= params.min_conversion }
    qc.filter { _sample, _cov, conversion -> conversion.toFloat() < params.min_conversion }
        .subscribe { sample, _cov, conversion ->
            log.warn "${sample} failed QC, conversion efficiency ${conversion}% < ${params.min_conversion}%"
        }

    // cohort summary
    summary = sampleSummary(
        dedup.join(qc_pass),
        file(params.panel)
    )
    cohort_summary = summary.collectFile(name: 'cohort_summary.tsv', keepHeader: true, skip: 1, sort: true)

    publish:
    dedup_bams = dedup
    cohort = cohort_summary
}

output {
    dedup_bams {
        path 'bams'
    }
    cohort {
        path '.'
    }
}
