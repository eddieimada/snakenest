rule VC_create_dictIndex:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
    output:
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        fai="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai"
    threads: 2
    resources:
        mem_mb=4000,
        runtime_min="00:30:00"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk CreateSequenceDictionary -R {input.reference}
        samtools faidx {input.reference}
        """

rule VC_create_intervalList:
    input:
        gtf="resources/Homo_sapiens.GRCh38.104.gtf",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict"
    output:
        tmp=temp("results/variantCalling/intervalList/exome.bed"),
        intervals="results/variantCalling/intervalList/exons.interval_list"
    threads: 2
    resources:
        mem_mb=4000,
        runtime_min="00:30:00"
    conda:
        "../envs/variantCalling.yaml"
    script:
        "../scripts/createInterval.R"

rule VC_create_uBam:
    input:
        f1="results/trim/{sample}/{sample}_{unit}_R1_trimmed.fastq.gz",
        f2="results/trim/{sample}/{sample}_{unit}_R2_trimmed.fastq.gz"
    output:
        ubam="results/variantCalling/ubams/{sample}_{unit}.ubam"
    params:
        compression=config["compression_level"],
        RG="{sample}"
    threads: 2
    resources:
        mem_mb=config["mem_create_uBam"],
        runtime_min=config["rt_create_uBam"]
    benchmark:
        "benchmark/create_uBam/{sample}_{unit}.tsv"
    log:
        "logs/create_uBam/{sample}_{unit}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        SM=$(basename {input.f1} |  sed 's/_L00.*//g')
        gatk FastqToSam \
        -F1 {input.f1} \
        -F2 {input.f2} \
        -O {output.ubam} \
        -RG {params.RG} \
        -SM $SM \
        -PL illumina \
        --COMPRESSION_LEVEL {params.compression} 2>> {log}
        """


rule VC_mergeuBams:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        refIndex="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai",
        ubam="results/variantCalling/ubams/{sample}_{unit}.ubam",
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam"
    output:
        merged=temp("results/variantCalling/mergedBam/{sample}_{unit}.bam")
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_mergeBams"],
        runtime_min=config["rt_mergeBams"]
    benchmark:
        "benchmark/mergeBams/{sample}_{unit}.tsv"
    log:
        "logs/mergeBams/{sample}_{unit}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        MergeBamAlignment \
        --REFERENCE_SEQUENCE {input.reference} \
        --UNMAPPED_BAM {input.ubam} \
        --ALIGNED_BAM {input.bam} \
        --OUTPUT {output.merged} \
        --INCLUDE_SECONDARY_ALIGNMENTS false \
        --VALIDATION_STRINGENCY SILENT \
        --COMPRESSION_LEVEL {params.compression} 2>> {log}
        """

rule VC_concatBam:
    input:
        bams=VC_gather_bams
    output:
        mbam="results/variantCalling/concatBam/{sample}.Aligned.sortedByCoord.out.bam"
    threads: 2
    resources:
        mem_mb=config["mem_concat"],
        runtime_min=config["rt_concat"]
    benchmark:
        "benchmark/mergedBam/{sample}.tsv"
    log:
        "logs/mergedBam/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        INPUT=({input.bams})
        if ((${{#INPUT[@]}} == 1)); then
            cp {input.bams} {output.mbam}
        else
            samtools merge -f -1 {output.mbam} {input.bams} 2>>{log}
        fi
        """

rule VC_markDuplicates:
    input:
        bam="results/variantCalling/concatBam/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        dedup=temp("results/variantCalling/markDuplicates/{sample}.bam"),
        index=temp("results/variantCalling/markDuplicates/{sample}.bai"),
        metrics="results/variantCalling/metrics/{sample}_dedup.metrics"
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_markDuplicates"],
        runtime_min=config["rt_markDuplicates"]
    benchmark:
        "benchmark/markDuplicates/{sample}.tsv"
    log:
        "logs/markDuplicates/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        MarkDuplicates \
        --INPUT {input.bam} \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --METRICS_FILE {output.metrics} \
        --OUTPUT {output.dedup} \
        --COMPRESSION_LEVEL {params.compression} 2>> {log}
        """

rule VC_splitNCigars:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        refIndex="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai",
        bam="results/variantCalling/markDuplicates/{sample}.bam"
    output:
        sacr=temp("results/variantCalling/splitNCigars/{sample}.bam"),
        index=temp("results/variantCalling/splitNCigars/{sample}.bai")
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_splitNCigars"],
        runtime_min=config["rt_splitNCigars"]
    benchmark:
        "benchmark/splitNCigars/{sample}.tsv"
    log:
        "logs/splitNCigars/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        SplitNCigarReads \
        -R {input.reference} \
        -I {input.bam} \
        -O {output.sacr} 2>> {log}
        """

rule VC_baseRecalibrator:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        refIndex="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai",
        bam="results/variantCalling/splitNCigars/{sample}.bam",
        bai="results/variantCalling/splitNCigars/{sample}.bai",
        dbSNP="resources/Homo_sapiens_assembly38.dbsnp138.vcf",
        knowIndels="resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    output:
        table="results/variantCalling/recalibration/{sample}.tbl"
    threads: 2
    resources:
        mem_mb=config["mem_baseRecalibrator"],
        runtime_min=config["rt_baseRecalibrator"]
    benchmark:
        "benchmark/baseRecalibrator/{sample}.tsv"
    log:
        "logs/baseRecalibrator/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -XX:ParallelGCThreads=1 \
        -Xms4000m" \
        BaseRecalibrator \
        -R {input.reference} \
        -I {input.bam} \
        --use-original-qualities \
        -O {output.table} \
        -known-sites {input.dbSNP} \
        -known-sites {input.knowIndels} 2>> {log}
        """

rule VC_applyBQSR:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        refIndex="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai",
        bam="results/variantCalling/splitNCigars/{sample}.bam",
        bai="results/variantCalling/splitNCigars/{sample}.bai",
        table="results/variantCalling/recalibration/{sample}.tbl"
    output:
        recalbam="results/variantCalling/recalibration/{sample}.bam",
        bai="results/variantCalling/recalibration/{sample}.bai"
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_applyBQSR"],
        runtime_min=config["rt_applyBQSR"]
    benchmark:
        "benchmark/applyBQSR/{sample}.tsv"
    log:
        "logs/applyBQSR/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
        -XX:+PrintGCDetails -XX:ParallelGCThreads=1 \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
        ApplyBQSR \
        --add-output-sam-program-record \
        -R {input.reference} \
        -I {input.bam} \
        --use-original-qualities \
        -O {output.recalbam} \
        --bqsr-recal-file {input.table} 2>> {log}
        """

rule VC_haplotypeCaller:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        refIndex="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai",
        bam="results/variantCalling/recalibration/{sample}.bam",
        bai="results/variantCalling/recalibration/{sample}.bai",
        dbSNP="resources/Homo_sapiens_assembly38.dbsnp138.vcf",
        intervals="results/variantCalling/intervalList/exons.interval_list"
    output:
        vcf=temp("results/variantCalling/vcf/{sample}.vcf.gz")
    threads: 2
    resources:
        mem_mb=config["mem_haplotypeCaller"],
        runtime_min=config["rt_haplotypeCaller"]
    benchmark:
        "benchmark/haplotypeCaller/{sample}.tsv"
    log:
        "logs/haplotypeCaller/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=1" \
        HaplotypeCaller \
        -R {input.reference} \
        -I {input.bam} \
        -L {input.intervals} \
        -O {output.vcf} \
        -dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        --dbsnp {input.dbSNP} 2>> {log}
        """

rule VC_filterVCF:
    input:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        refDict="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.dict",
        refIndex="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai",
        vcf="results/variantCalling/vcf/{sample}.vcf.gz"
    output:
        vcf_f="results/variantCalling/vcf_filtered/{sample}.vcf.gz",
        vcfIdx="results/variantCalling/vcf_filtered/{sample}.vcf.gz.tbi"
    threads: 2
    resources:
        mem_mb=config["mem_haplotypeCaller"],
        runtime_min=config["rt_haplotypeCaller"]
    benchmark:
        "benchmark/filterVCF/{sample}.tsv"
    log:
        "logs/filterVCF/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        VariantFiltration \
        --R {input.reference} \
        --V {input.vcf} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O {output.vcf_f}
        """    
