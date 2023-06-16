SAMPLES = ["BELA_1_2004","SMAL_6_2013"] # Add or remove samples here

rule all:
    input:
        expand("results/preprocess_fastqc/{sample}_output/{sample}.R1_fastqc.html", sample=SAMPLES),
        expand("results/preprocess_fastqc/{sample}_output/{sample}.R2_fastqc.html", sample=SAMPLES),
        expand("results/fastp/{sample}_1.fastp.fastq.gz", sample=SAMPLES),
        expand("results/fastp/{sample}_2.fastp.fastq.gz", sample=SAMPLES),
        expand("results/fastp/{sample}.fastp.json", sample=SAMPLES),
        expand("results/fastp/{sample}.fastp.html", sample=SAMPLES),
        expand("results/bam/{sample}.bam", sample=SAMPLES)
        "bam/{sample}.bam.bai"


rule fastqc:
    input:
        r1 = "data/samples/contemporary_set/{sample}.R1.fastq.gz",
        r2 = "data/samples/contemporary_set/{sample}.R2.fastq.gz"
    conda:
        "environment.yaml"
    output:
        html_report1 = "results/preprocess_fastqc/{sample}_output/{sample}.R1_fastqc.html"
        html_report2 = "results/preprocess_fastqc/{sample}_output/{sample}.R2_fastqc.html"
        #zip_archive = "results/preprocess_fastqc/{sample}_output.zip"
    shell:
        """
        mkdir -p results/preprocess_fastqc/{wildcards.sample}_output
        fastqc {input.r1} {input.r2} --threads 4 --outdir results/preprocess_fastqc/{wildcards.sample}_output
        """

rule fastp:
    input:
        r1 = "data/samples/contemporary_set/{sample}.R1.fastq.gz",
        r2 = "data/samples/contemporary_set/{sample}.R2.fastq.gz"
    conda:
        "environment.yaml"
    output:
        r1_out = "results/fastp/{sample}_1.fastp.fastq.gz",
        r2_out = "results/fastp/{sample}_2.fastp.fastq.gz",
        html_out = "results/fastp/{sample}.fastp.html"
    params:
        thread = 12,
        reads_to_process = 1000
    shell:
        """
        mkdir -p results/fastp/
        fastp \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --out1 {output.r1_out} \
        --out2 {output.r2_out} \
        --html {output.html_out} \
        --thread {params.thread} \
        --overrepresentation_analysis \
        --reads_to_process {params.reads_to_process} \
        --detect_adapter_for_pe
        """

rule bwa:
    input:
        ref="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna",
        r1="results/fastp/{sample}_1.fastp.fastq.gz",
        r2="results/fastp/{sample}_2.fastp.fastq.gz"
    output:
        "results/bam/{sample}.bam"
    log:
        "logs/bwa/{sample}.log"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:Illumina"
    threads: 2
    shell:
        """
        mkdir -p results/bam/
        mkdir -p logs/bwa/
        mkdir -p tmp/bwa/{wildcards.sample}
        bwa mem -R '{params.rg}' -t {threads} -M {input.ref} {input.r1} {input.r2} \
        |samtools sort -m 6G -@{threads} -T tmp/bwa/{wildcards.sample} - >{output}.tmp) 2>{log} \
        && mv {output}.tmp {output}"
         """


rule bwa:
    input:
        ref="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna",
        r1="results/fastp/{sample}_1.fastp.fastq.gz",
        r2="results/fastp/{sample}_2.fastp.fastq.gz"
    output:
        "results/bam/{sample}.bam"
    log:
        "logs/bwa/{sample}.log"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:Illumina"
    threads: 20
    shell:
        """
        mkdir -p results/bam/
        mkdir -p logs/bwa/
        mkdir -p tmp/bwa/{wildcards.sample}
        bwa mem -R '{params.rg}' -t {threads} -M {input.ref} {input.r1} {input.r2} \
        |samtools sort -m 6G -@{threads} -T tmp/bwa/{wildcards.sample} - >{output}.tmp) 2>{log} \
        && mv {output}.tmp {output}"
         """

rule index_bam:
    input:
        "results/bam/{sample}.bam"
    output:
        "results/bam/{sample}.bam.bai"
    log:
        "logs/samtools/index.{sample}.log"
    threads: 1
    shell:
        "samtools index {input} {output}"

rule markdupl:
    input:
        "results/bam/{sample}.bam"
    output:
        bam="results/bam/{sample}.md.bam",
        metrics="results/bam/{sample}.metrics"
    log:
        "logs/picard/{sample}.mkdupl.log"
    params:
        picard="/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar",
        mem="60g"
    threads: 10
    shell:
        "(java -Xmx{params.mem} -jar {params.picard} MarkDuplicates INPUT={input} \
         METRICS_FILE={output.metrics} TMP_DIR=$SNIC_TMP ASSUME_SORTED=true \
         VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT={output.bam}.tmp) 2>{log} \
         && mv {output.bam}.tmp {output.bam} && mv {output.bam}.tmp.bai {output.bam}.bai"

rule qualimap:
    input:
        "results/bam/{sample}.bam"
    output:
        "qualimap_results/{sample}_qualimap"
    params:
        #extra="-nt 8" # optional extra parameters
    log:
        "logs/qualimap/{sample}.log"
    shell:
        """
        mkdir -p results/bam/qualimap/{wildcards.sample}_qualimap
        qualimap bamqc -bam {input} -outdir qualimap_results/{wildcards.sample}_qualimap -c -sd {params.extra} 2> {log}
        """

















# The following rules will not be executed as jobs
localrules: all

rule all:
    input:

        "vep/74Females.SNPs.HF.mac2.txt",
        "vep/74Females.SNPs.HF.mac1.txt"

SAMPLES = ["AfGWo", "EthWo", "Dhole", "BlBJa", "SiSJa"]
FILTER = ["mac2", "mac1"]

# The following rules will not be executed as jobs
localrules: all, makeGVCFList

rule all:
    input:
        expand("outgroup/{prefix}.100S95F14R.chr1-38.{prefix2}.extract.vcf", prefix=SAMPLES, prefix2=FILTER),
        expand("outgroup/{prefix}.74Females.chrX.{prefix2}.extract.vcf", prefix=SAMPLES, prefix2=FILTER)

rule bwa:
    input:
        ref="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        r1="fastq/{prefix}_1.fastq.gz",
        r2="fastq/{prefix}_2.fastq.gz"
    output:
        "bam/{prefix}.bam"
    log:
        "logs/bwa/{prefix}.log"
    params:
        rg=r"@RG\tID:{prefix}\tSM:{prefix}\tLB:{prefix}\tPL:Illumina"
    threads: 20
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} -M {input.ref} {input.r1} {input.r2} \
         |samtools sort -m 6G -@{threads} -T $SNIC_TMP/{wildcards.prefix} - >{output}.tmp) 2>{log} \
         && mv {output}.tmp {output}"

rule index_bam:
    input:
        "bam/{prefix}.bam"
    output:
        "bam/{prefix}.bam.bai"
    log:
        "logs/samtools/index.{prefix}.log"
    threads: 1
    shell:
        "samtools index {input} {output}"

rule mkdupl:
    input:
        "bam/{prefix}.bam"
    output:
        bam="bam/{prefix}.md.bam",
        metrics="bam/{prefix}.metrics"
    log:
        "logs/picard/{prefix}.mkdupl.log"
    params:
        picard="/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar",
        mem="60g"
    threads: 10
    shell:
        "(java -Xmx{params.mem} -jar {params.picard} MarkDuplicates INPUT={input} \
         METRICS_FILE={output.metrics} TMP_DIR=$SNIC_TMP ASSUME_SORTED=true \
         VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT={output.bam}.tmp) 2>{log} \
         && mv {output.bam}.tmp {output.bam} && mv {output.bam}.tmp.bai {output.bam}.bai"

rule haplotype_caller:
    input:
        bam="bam/{prefix}.md.bam",
        ref="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa"
    output:
        "gvcf/{prefix}.g.vcf.gz"
    log:
        "logs/GATK/HaplotypeCaller.{prefix}.log"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="120g"
    threads: 20
    shell:
        "(java -Xmx{params.mem} -Djava.io.tmpdir=$SNIC_TMP -jar {params.gatk} \
        -T HaplotypeCaller -R {input.ref} -I {input.bam} --emitRefConfidence GVCF \
        --variant_index_type LINEAR --variant_index_parameter 128000 -nct {threads} \
        -jdk_deflater -jdk_inflater -o gvcf/{wildcards.prefix}.g.vcf.tmp.gz) 2>{log} \
        && mv gvcf/{wildcards.prefix}.g.vcf.tmp.gz {output} \
        && mv gvcf/{wildcards.prefix}.g.vcf.tmp.gz.tbi {output}.tbi"

rule makeGVCFList:
    input:
        "gvcf/{prefix}.g.vcf.gz"
    output:
        "{prefix}.gvcf.list"
    log:
        "logs/misc/create.GVCFList.{prefix}.log"
    shell:
        "echo {input} >{output}"

rule AllsitesGenotypeGVCFs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        list="{prefix}.gvcf.list"
    output:
        "vcf/{prefix}.vcf"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="120g"
    log:
        "logs/GATK/GenotypeGVCFs.{prefix}.log"
    threads: 20
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=$SNIC_TMP -jar {params.gatk} \
        -T GenotypeGVCFs -R {input.fa} --variant {input.list} -nt {threads} \
        -allSites -o {output}.tmp && mv {output}.tmp {output}"

rule zip:
    input:
        "vcf/{prefix}.vcf"
    output:
        "vcf/{prefix}.vcf.gz"
    log:
        "logs/zip/{prefix}.log"
    threads: 1
    shell:
        "bgzip {input}"

rule tabix:
    input:
        "vcf/{prefix}.vcf.gz"
    output:
        "vcf/{prefix}.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/{prefix}.log"
    threads: 1
    shell:
        "tabix {params} {input}"

rule extractBedSites:
    input:
        vcf="vcf/{species, \w+}.vcf.gz",
        index="vcf/{species, \w+}.vcf.gz.tbi",
        bed="bed/{filter}.bed"
    output:
        "outgroup/{species, \w+}.{filter}.extract.vcf"
    log:
        "logs/vcf/{species, \w+}.{filter}.extract.log"
    threads: 1
    shell:
        "zcat {input.vcf} |grep -v '<NON_REF>' | \
        intersectBed -header -a - -b {input.bed} >{output}"




rule GenotypeGVCFs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        list="{prefix}.gvcf.list"
    output:
        "vcf/{prefix}.vcf"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="1000g"
    log:
        "logs/GATK/GenotypeGVCFs.{prefix}.log"
    threads: 20
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=$SNIC_TMP -jar {params.gatk} \
        -T GenotypeGVCFs -R {input.fa} --variant {input.list} -nt {threads} -o {output}.tmp && mv {output}.tmp {output}"

rule zip:
    input:
        "vcf/{prefix}.vcf"
    output:
        "vcf/{prefix}.vcf.gz"
    log:
        "logs/zip/{prefix}.log"
    threads: 1
    shell:
        "bgzip {input}"

rule tabix:
    input:
        "vcf/{prefix}.vcf.gz"
    output:
        "vcf/{prefix}.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/{prefix}.log"
    threads: 1
    shell:
        "tabix {params} {input}"

rule idepth:
    input:
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz"
    output:
        "coverage/{prefix}.vcf.idepth"
    log:
        "logs/coverage/{prefix}.log"
    threads: 1
    shell:
        "vcftools --gzvcf {input.vcf} --stdout --depth >{output}"

rule extractSNPs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz"
    output:
        "vcf/{prefix}.SNPs.vcf"
    log:
        "logs/GATK/{prefix}.extractSNPs.log"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="60g"
    threads: 10
    shell:
        "java -Xmx{params.mem} -jar {params.gatk} -T SelectVariants \
        -R {input.fa} -V {input.vcf} -selectType SNP -o {output}"


rule filterSNPs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        index="vcf/{prefix}.SNPs.vcf.gz.tbi",
        vcf="vcf/{prefix}.SNPs.vcf.gz"
    output:
        "vcf/{prefix}.SNPs.HF.vcf"
    log:
        "logs/GATK/{prefix}.HardFiltSNPs.log"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="60g",
        filt=r"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    threads: 10
    shell:
        "java -Xmx{params.mem} -jar {params.gatk} -T VariantFiltration \
        -R {input.fa} -V {input.vcf} --filterExpression '{params.filt}' \
        --filterName 'hard_filt' -o {output}"


rule vcfFilter:
    input:
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz",
        dep="coverage/{prefix}.vcf.idepth"
    output:
        "vcf/{prefix}.mac2.vcf"
    log:
        "logs/filter/{prefix}.mac2.log"
    params:
        min=10,
        minal=2,
        maxal=2,
        mac=2,
        maxmiss=0.9,
        GQ=30,
    threads: 1
    shell:
        "m=`awk '(NR>1){{n++; sum+=$3}}END{{mean=sum/n; twice=2*mean; print twice}}' {input.dep}` && "
        "vcftools --gzvcf {input.vcf} --out vcf/{wildcards.prefix} --mac {params.mac}  \
        --max-alleles {params.maxal}  --min-alleles {params.minal}  --min-meanDP {params.min} \
        --max-meanDP $m  --minGQ {params.GQ} --max-missing {params.maxmiss} --remove-filtered-all \
        --recode --recode-INFO-all && mv vcf/{wildcards.prefix}.recode.vcf {output}"

rule keepSingletonFilter:
    input:
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz",
        dep="coverage/{prefix}.vcf.idepth"
    output:
        "vcf/{prefix}.mac1.vcf"
    log:
        "logs/filter/{prefix}.mac1.log"
    params:
        min=10,
        minal=2,
        maxal=2,
        maxmiss=0.9,
        GQ=30,
    threads: 1
    shell:
        "m=`awk '(NR>1){{n++; sum+=$3}}END{{mean=sum/n; twice=2*mean; print twice}}' {input.dep}` && "
        "vcftools --gzvcf {input.vcf} --out vcf/{wildcards.prefix}  \
        --max-alleles {params.maxal}  --min-alleles {params.minal}  --min-meanDP {params.min} \
        --max-meanDP $m  --minGQ {params.GQ} --max-missing {params.maxmiss} --remove-filtered-all \
        --recode --recode-INFO-all && mv vcf/{wildcards.prefix}.recode.vcf {output}"

rule vep:
    input:
        vcf="vcf/{prefix}.vcf.gz",
        index="vcf/{prefix}.vcf.gz.tbi",
    output:
        "vep/{prefix}.txt"
    log:
        "logs/vep/{prefix}.log"
    params:
        sp="canis_familiaris"
    threads: 4
    shell:
        "vep --cache --dir $VEP_CACHE -i {input.vcf} -o {output} \
        --species {params.sp} --fork {threads} --force_overwrite --sift b"
