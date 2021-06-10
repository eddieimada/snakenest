rule download_resources:
    output:
        reference="resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        transcripts="resources/Homo_sapiens.GRCh38.transcripts.fa.gz",
        tx2gene="resources/tx2gene.tsv.gz",
        gtf="resources/Homo_sapiens.GRCh38.104.gtf",
        paSites="resources/atlas.clusters.2.0.GRCh38.96.bed",
        fc="resources/FC_robust.saf",
        blacklist="resources/blacklist_hg38_GRCh38_v2.1.0.tsv.gz",
        knowFusions="resources/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz",
        proteinDomains="resources/protein_domains_hg38_GRCh38_v2.1.0.gff3",
        dbSNP="resources/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbSNPidx="resources/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        knowIndels="resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        knowIndelsidx="resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
    priority: 3
    threads: 1
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/downloads/resources.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        set -e
        echo "Downloading resources..."
        [[ -f {output.dbSNP} ]] || gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf {output.dbSNP}
        [[ -f {output.dbSNPidx} ]] || gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx {output.dbSNPidx}
        [[ -f {output.knowIndels} ]] || gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz {output.knowIndels}
        [[ -f {output.knowIndelsidx} ]] || gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi {output.knowIndelsidx}
        [[ -f {output.reference} ]] || wget -O resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz {output.reference}
        [[ -f {output.paSites} ]] || wget -O resources/atlas.clusters.2.0.GRCh38.96.bed.gz https://www.dropbox.com/s/a9ffcgsm9in0qzq/atlas.clusters.2.0.GRCh38.96.bed.gz
        [[ -f {output.tx2gene} ]] || wget -O resources/tx2gene.tsv.gz https://www.dropbox.com/s/bltd39gcs7603st/tx2gene_hg38.tsv.gz
        [[ -f {output.fc} ]] || wget -O resources/FC_robust.saf.gz https://www.dropbox.com/s/oy3n86poeoejndc/FC_robust.saf.gz
        [[ -f {output.gtf} ]] || wget -O resources/Homo_sapiens.GRCh38.104.gtf.gz http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
        [[ -f {output.transcripts} ]] || wget -O resources/Homo_sapiens.GRCh38.cdna.all.fa.gz http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        [[ -f {output.transcripts} ]] || wget -O resources/Homo_sapiens.GRCh38.ncrna.fa.gz http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
        [[ -f {output.blacklist} ]] || wget -O resources/arriba_v2.1.0.tar.gz https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz
        tar -xvf resources/arriba_v2.1.0.tar.gz --directory resources/
        mv resources/arriba_v2.1.0/database/blacklist_hg38_GRCh38_v2.1.0.tsv.gz resources/arriba_v2.1.0/database/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz resources/arriba_v2.1.0/database/protein_domains_hg38_GRCh38_v2.1.0.gff3 resources/
        rm -rf resources/arriba_v2.1.0/ resources/arriba_v2.1.0.tar.gz
        [[ -f {output.transcripts} ]] || cat resources/Homo_sapiens.GRCh38.cdna.all.fa.gz resources/Homo_sapiens.GRCh38.ncrna.fa.gz > {output.transcripts}
        echo "Unzipping some files..."
        gunzip resources/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz resources/FC_robust.saf.gz resources/atlas.clusters.2.0.GRCh38.96.bed.gz resources/Homo_sapiens.GRCh38.104.gtf.gz
        echo "Done!"
        """

rule star_index_download:
    output:
        index = "resources/STARindex_hg19/SA"
    priority: 3
    threads: 1
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/STARindex/STARindex.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        set -e
        wget -O resources/STARindex_hg19 https://www.dropbox.com/s/7zf1ad0hokpizfi/STARindex_hg19.tar
        tar -xvf resources/STARindex_hg19 --directory resources/
        """

rule salmon_index_download:
    output:
        "resources/salmon_hg19/ctable.bin"
    threads: 1
    priority: 3
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/salmonIndex/salmonIndex.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        set -e
        wget -O resources/salmon_hg19.tar https://www.dropbox.com/s/mybii6z71v9pt57/salmon_hg19.tar
        tar -xvf resources/salmon_hg19.tar --directory resources/
        """

rule download_ctatLib:
    output:
        ctat = "resources/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
    threads: 1
    priority: 2
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/STARfusion/download_CTAT.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        wget -O resources/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
        tar -xzvf resources/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz --directory resources/
        """
