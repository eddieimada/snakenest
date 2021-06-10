rule rseqc_gtf2bed:
    input:
        gtf="resources/Homo_sapiens.GRCh38.104.gtf",
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db")
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_gtf2bed.log",
    conda:
        "../envs/RNAseq.yaml"
    script:
        "../scripts/gtf2bed.py"

rule fastqc:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/fastqc/{sample}/{sample}_{unit}_fastqc.zip"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/fastqc/{sample}_{unit}.log",
    params:
        dir="results/qc/fastqc/{sample}/",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        fastqc -o {params.dir} {input.bam}
        mv results/qc/fastqc/{wildcards.sample}/{wildcards.sample}_{wildcards.unit}.Aligned.sortedByCoord.out_fastqc.zip {output} 
        """

rule rseqc_junction_annotation:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="results/qc/rseqc/{sample}_{unit}.junctionanno",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",
        prefix="results/qc/rseqc/{sample}_{unit}.junctionsat",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.stats.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_stat/{sample}_{unit}.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.infer_experiment.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_infer/{sample}_{unit}.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_innerdis/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.inner_distance_freq",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdistribution.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_readdis/{sample}_{unit}.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_readdup/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.readdup",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_readgc/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.readgc",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule bamqc:
    input:
        expand("results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.infer_experiment.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.stats.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.readdistribution.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",zip, sample=samples_trim, unit=unit_trim),
        expand("results/STAR_2p/{sample}_{unit}Log.final.out",zip, sample=samples_trim, unit=unit_trim)
    output:
        "results/qc/bam_qc.html"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    params:
    	outDir="results/qc",
        name="bam_qc.html"
    log:
        "logs/bam_qc.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "multiqc"
        " --force"
        " -o {params.outDir}"
        " -n {params.name}"
        " {input}"

rule sample_qc:
    input:
        expand("results/counts/featureCounts/{sample}_geneCounts_gencode.tsv.summary", sample=SAMPLES),
        expand("results/counts/featureCounts/{sample}_geneCounts_fc.tsv.summary", sample=SAMPLES),
        expand("results/salmon/{sample}/aux_info/meta_info.json",sample=SAMPLES),
        expand("results/salmon/{sample}/libParams/flenDist.txt",sample=SAMPLES) 
    output:
        "results/qc/samples_qc.html"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    params:
    	outDir="results/qc",
        name="samples_qc.html"
    log:
        "logs/samples_qc.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "multiqc"
        " --force"
        " -o {params.outDir}"
        " -n {params.name}"
        " {input}"

rule reads_qc:
    input:
        expand("results/qc/fastqc/{sample}/{sample}_{unit}_fastqc.zip",zip, sample=samples_trim, unit=unit_trim)
    output:
        "results/qc/reads_qc.html"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    params:
    	outDir="results/qc",
        name="reads_qc.html"
    log:
        "logs/read_qc.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "multiqc"
        " --force"
        " -o {params.outDir}"
        " -n {params.name}"
        " {input}"