rule nanopolishIndexing:
    """
    Creates index file that links read ids  with signal level data in
     FAST5 files.
    """
    input:
        fq = fq / "{sample}.fastq",
        f5 = f5
    output:
        index = fq / "{sample}.fastq.index",
        fai = fq / "{sample}.fastq.index.fai",
        gzi = fq / "{sample}.fastq.index.gzi",
        readdb = fq / "{sample}.fastq.index.readdb"
    threads: 1
    container: "https://depot.galaxyproject.org/singularity/nanopolish:0.13.1--ha077697_0"
    log: "results/logs/{sample}_index.log"
    shell:
        """
        (nanopolish index -d {input.f5} {input.fq}) >{log} 2>&1
        """

rule minimap2Align:
    """
    Align reads to the genome.
    """
    input:
        reference = config["ref"],
        fq = fq / "{sample}.fastq"
    output:
        sam = temp("results/minimap2/{sample}.sam")
    threads: 5
    container: "https://depot.galaxyproject.org/singularity/minimap2:2.24--h7132678_1"
    log: "results/logs/{sample}_mapping.log"
    shell:
        """
        (
        minimap2 -t {threads} -a -x map-ont {input.reference} {input.fq} -o {output.sam}
        ) >{log} 2>&1
        """

rule samtools:
    """
    Convert SAM to BAM.
    Sorts BAM.
    Creates sorted BAM index.
    """
    input:
        sam = rules.minimap2Align.output.sam
    output:
        sorted = "results/minimap2/{sample}.bam",
        bai = "results/minimap2/{sample}.bam.bai"
    threads: 6
    container: "docker://quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    log: "results/logs/{sample}_samtools.log"
    shell:
        """
        (
        samtools view -Sb {input.sam} | \
        samtools sort -o {output.sorted}
        samtools index {output.sorted}
        ) >{log} 2>&1
        """

rule nanopolishCallMethyl:
    """
    Call methylation
    """
    input:
        fq = fq / "{sample}.fastq",
        bam = rules.samtools.output.sorted,
        reference = config["ref"],
        bai = rules.nanopolishIndexing.output.index,
        fai = rules.nanopolishIndexing.output.fai,
        gzi = rules.nanopolishIndexing.output.gzi,
        readdb = rules.nanopolishIndexing.output.readdb
    output:
        meth = "results/nanopolish/{sample}_methylation.tsv"
    threads: 5
    container: "https://depot.galaxyproject.org/singularity/nanopolish:0.13.1--ha077697_0"
    log: "results/logs/{sample}_methylation.log"
    shell:
        """
        (
        nanopolish call-methylation -t {threads} \
        -r {input.fq} -b {input.bam} -g {input.reference} > {output.meth}
        ) >{log} 2>&1
        """

rule nanopolishCalcMethylFreq:
    """
    Calculate methylation frequency from called methylation.
    """
    input:
        tsv = rules.nanopolishCallMethyl.output.meth
    output:
        freq = "results/nanopolish/{sample}_frequency.tsv"
    threads: 1
    container: "https://depot.galaxyproject.org/singularity/nanopolish:0.13.1--ha077697_0"
    log: "results/logs/{sample}_frequency.log"
    shell:
        """
        (
        python3 workflow/scripts/nanopolish_scripts/calculate_methylation_frequency.py {input.tsv} > {output.freq}
        ) >{log} 2>&1
        """

