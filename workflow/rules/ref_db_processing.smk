rule map_reads_for_ref_db:
    input:
        reads=get_fastqs,
        idx=get_bwa_index,
    output:
        "results/{date}/ref-db/mapped/ref~{reference}/{sample}.bam",
    log:
        "logs/{date}/ref-db/bwa-mem/ref~{reference}/{sample}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra="",
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.69.0/bio/bwa/mem"


rule extract_reads_of_interest_for_ref_db:
    input:
        "results/{date}/ref-db/mapped/ref~main+human/{sample}.bam",
    output:
        "results/{date}/ref-db/mapped/ref~main+human/nonhuman/{sample}.bam",
    log:
        "logs/{date}/ref-db/extract_reads_of_interest/{sample}.log",
    threads: 1
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract-reads-of-interest.py"


rule order_nonhuman_reads_for_ref_db:
    input:
        "results/{date}/ref-db/mapped/ref~main+human/nonhuman/{sample}.bam",
    output:
        fq1="results/{date}/ref-db/nonhuman-reads/{sample}.1.fastq.gz",
        fq2="results/{date}/ref-db/nonhuman-reads/{sample}.2.fastq.gz",
        bam_sorted="results/{date}/ref-db/nonhuman-reads/{sample}.sorted.bam",
    log:
        "logs/{date}/ref-db/order_nonhuman_reads/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 8
    shell:
        """
        samtools sort  -@ {threads} -n {input} -o {output.bam_sorted} > {log} 2>&1
        samtools fastq -@ {threads} {output.bam_sorted} -1 {output.fq1} -2 {output.fq2} >> {log} 2>&1
        """


rule sort_bam_for_ref_db:
    input:
        expand(
            "results/{{date}}/ref-db/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["adapters"]["amplicon-reference"],
        ),
    output:
        "results/{date}/ref-db/clipped-reads/{sample}.bam",
    log:
        "logs/{date}/ref-db/sort-bam/{sample}.log",
    params:
        extra="-m 4G",
        tmp_dir="/tmp/",
    threads: 8
    wrapper:
        "0.74.0/bio/samtools/sort"


rule clip_primer_for_ref_db:
    input:
        sortbam="results/{date}/ref-db/clipped-reads/{sample}.bam",
        sortindex="results/{date}/ref-db/clipped-reads/{sample}.bam.bai",
        bed=config["adapters"]["amplicon-primers"],
        ref_fasta="resources/genomes/{reference}.fasta".format(
            reference=config["adapters"]["amplicon-reference"]
        ),
    output:
        clippedbam=temp(
            "results/{date}/ref-db/clipped-reads/{sample}.primerclipped.bam"
        ),
        hardclippedbam=temp(
            "results/{date}/ref-db/clipped-reads/{sample}.primerclipped.hard.bam"
        ),
        sorthardclippedbam=temp(
            "results/{date}/ref-db/clipped-reads/{sample}.primerclipped.hard.sorted.bam"
        ),
        fq1="results/{date}/ref-db/clipped-reads/{sample}.1.fastq.gz",
        fq2="results/{date}/ref-db/clipped-reads/{sample}.2.fastq.gz",
    log:
        "logs/{date}/ref-db/primer-clipping/{sample}.log",
    params:
        dir=lambda w, input: os.path.dirname(input.sortbam),
        bam=lambda w, input: input.sortbam.split("/")[-1],
        dir_depth=lambda w, input: "".join(
            ["../"] * (len(input.sortbam.split("/")) - 1)
        ),
    conda:
        "../envs/bamclipper.yaml"
    threads: 10
    shell:
        """
        cd {params.dir}
        bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} > {params.dir_depth}{log} 2>&1
        cd {params.dir_depth}
        fgbio --sam-validation-stringency=LENIENT ClipBam -i {output.clippedbam} -o {output.hardclippedbam} -H true -r {input.ref_fasta} >> {log} 2>&1
        samtools sort  -@ {threads} -n {output.hardclippedbam} -o {output.sorthardclippedbam}  >> {log} 2>&1
        samtools fastq -@ {threads} {output.sorthardclippedbam} -1 {output.fq1} -2 {output.fq2}  >> {log} 2>&1
        """


rule rename_mv_reads_for_ref_db:
    input:
        fastq1=lambda wildcards: get_reads_after_qc_for_ref_db(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc_for_ref_db(wildcards, read="2"),
    output:
        fastq1="results/{date}/ref-db/preprocessed-reads/{sample}.1.fastq.gz",
        fastq2="results/{date}/ref-db/preprocessed-reads/{sample}.2.fastq.gz",
    log:
        "logs/{date}/ref-db/rename-mv-for-ref-db/{sample}.log",
    threads: 1
    conda:
        "../envs/unix.yaml"
    shell:
        """
        cp {input.fastq1} {output.fastq1} > {log} 2>&1
        cp {input.fastq2} {output.fastq2} >> {log} 2>&1
        """
