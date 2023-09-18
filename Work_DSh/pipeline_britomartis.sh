#Dialy log Dasha

#### Working directory ####
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation




###

ssh -AX <username>@rackham.uppmax.uu.se

ml bioinfo-tools dds-cli

dds data ls --project "ngisthlm00193"
     ︵
 ︵ (  )   ︵
(  ) ) (  (  )   SciLifeLab Data Delivery System
 ︶  (  ) ) (    https://delivery.scilifelab.se/
      ︶ (  )    CLI Version 2.2.62
          ︶
INFO     No saved token found, or token has expired, proceeding with authentication
DDS username: dariash
DDS password:
╭─ Error ──────────────────────────────────────────────────────────────────────────────────────────╮
│ Failed to authenticate user: Missing or incorrect credentials                                    │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
(base) [daria@rackham3 ~]$ dds data ls --project "ngisthlm00193"
     ︵
 ︵ (  )   ︵
(  ) ) (  (  )   SciLifeLab Data Delivery System
 ︶  (  ) ) (    https://delivery.scilifelab.se/
      ︶ (  )    CLI Version 2.2.62
          ︶
INFO     No saved token found, or token has expired, proceeding with authentication
DDS username: dariash
DDS password:
INFO     Please enter the one-time authentication code sent to your email address (leave empty to exit):
Authentication one-time code: 11775319
INFO     Listing files for project 'ngisthlm00193'

 Files / directories in project: ngisthlm00193
 └── P27562/



screen -S britomartis

ml bioinfo-tools dds-cli
cd /proj/uppstore2017185/b2014034/private/raw_data/Assmann
dds data get --get-all --project "ngisthlm00193"

#screen -S britomartis_arxiv
#ml bioinfo-tools dds-cli
#cd /proj/uppoff2020002/private/raw_data_backups/Assmann
#dds data get --get-all --project "ngisthlm00193"

screen -r britomartis

/proj/uppoff2020002/private/raw_data_backups/Assmann

###Getting backup copy

ls /proj/uppoff2020002/private/raw_data_backups/Assmann


#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J assmann_download
#SBATCH --mail-user=daria.shipilina@gmail.com

ml bioinfo-tools dds-cli
cd /proj/uppoff2020002/private/raw_data_backups/Assmann
dds data get --get-all --project "ngisthlm00193"


########### DOWNLOAD ATHALIA
## from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/220/545/GCA_905220545.2_ilMelAtha1.2/
mkdir reference_genome_M.athalia
cd reference_genome_M.athalia/

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/220/545/GCA_905220545.2_ilMelAtha1.2/GCA_905220545.2_ilMelAtha1.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/220/545/GCA_905220545.2_ilMelAtha1.2/GCA_905220545.2_ilMelAtha1.2_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/220/545/GCA_905220545.2_ilMelAtha1.2/README.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/220/545/GCA_905220545.1_ilMelAtha1.1/GCA_905220545.1_ilMelAtha1.1_feature_count.txt.gz

### thinking
 --save_output_as_bam


 #### Setting up the pipeline

 ## 1. Get sample NAMES (command line is fine)

 #SCRIPT FOR READING SEQUENCES NAMES FROM STANDARD DELIVERY FOLDER
 #This loop iterates through the delivery folder with standard UPPMAX structure and outputs soft links to fasta files
 for d in /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/*/ ; do
         cd "$d"
         for filename in */*/* ; do
                 echo "$d""$filename"
         done
 done > fasta_pathes_britomartis.txt


### Making test sample list manually

# S55, S53, S15, S20

declare -A groups=(["MUZAV5_43"]="P27562_1015" ["MUZAV5_98"]="P27562_1020" ["MUZAV90_05"]="P27562_1053" ["MUZAV26_94"]="P27562_1055")

read1=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Concervation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R1')
read2=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Concervation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R2')
echo "$files,$files,L001,$read1,$read2"
done

mkdir 00_Benchmarking4x
nano make_sample_list.sh
bash make_sample_list.sh > sample_sheet.csv

## Sample list ready

module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Concervation/00_Mapping_Calling_sarek/00_Benchmarking4x

nextflow run nf-core/sarek --input sample_sheet.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna.gz --skip_tools baserecalibrator --outdir ./results


#  Problem with reference, can't be gunzipped, restarting magical_williams

### Also renamed FOLDER
cd reference_genome_M.athalia/
gunzip GCA_905220545.2_ilMelAtha1.2_genomic.fna.gz

nextflow run nf-core/sarek --input sample_sheet.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -resume magical_williams --skip_tools baserecalibrator --outdir ./results

cp multiqc_report.html multiqc_report_160323.html





# Benchmarking step somehow rewrtitten

### Mapping of all the samples

#!/bin/bash -l
module load bioinfo-tools Nextflow nf-core nf-core-pipelines

#Confiruging Nextflow in the working directory
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll
NXF_OPTS='-Xms1g -Xmx4g'

#Providing full path to new data
#gtf="/crex/proj/uppstore2017185/b2014034/nobackup/Dasha/VanessaRNAseq/genome_data/vcard.gtf"
fasta="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna"
nextflow run nf-core/sarek --input sample_sheet.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna --skip_tools baserecalibrator --outdir ./results


### Making sample list
declare -A groups=(["КAZA_2_2006"]="P27562_1068" ["КAZA_1_2006"]="P27562_1038" ["КALM_9_1981"]="P27562_1056" ["КALM_8_1969"]="P27562_1029" ["КALM_7_1961"]="P27562_1028" ["КALM_6_1961"]="P27562_1025" ["КALM_5_1961"]="P27562_1021" ["КALM_4_1957"]="P27562_1024" ["КALM_3_1955"]="P27562_1031" ["КALM_2_1955"]="P27562_1030" ["КALM_11_1997"]="P27562_1036" ["КALM_10_1983"]="P27562_1057" ["АLTA_2_2015"]="P27562_1069" ["АLTA_1_2015"]="P27562_1051" ["VORO_1_1998"]="P27562_1071" ["VAST_3_2005"]="P27562_1055" ["VAST_2_1999"]="P27562_1066" ["VAST_1_1983"]="P27562_1065" ["URAL_2_1995"]="P27562_1061" ["URAL_1_1991"]="P27562_1062" ["UPPS_3_1969"]="P27562_1022" ["UPPS_2_1958"]="P27562_1023" ["UPPS_1_1951"]="P27562_1026" ["STOC_6_1965"]="P27562_1035" ["STOC_5_1965"]="P27562_1013" ["STOC_4_1965"]="P27562_1008" ["STOC_3_1965"]="P27562_1007" ["STOC_2_1965"]="P27562_1002" ["STOC_1_1965"]="P27562_1001" ["SMAL_6_2013"]="P27562_1044" ["SMAL_5_1998"]="P27562_1037" ["SMAL_4_1996"]="P27562_1054" ["SMAL_3_1967"]="P27562_1033" ["SMAL_2_1967"]="P27562_1020" ["SMAL_1_1967"]="P27562_1005" ["SLOV_1_1990"]="P27562_1059" ["RUSS_6_2008"]="P27562_1050" ["RUSS_5_2008"]="P27562_1042" ["RUSS_4_2008"]="P27562_1041" ["RUSS_3_1999"]="P27562_1060" ["RUSS_2_1999"]="P27562_1046" ["RUSS_1_1998"]="P27562_1067" ["POLA_2_2003"]="P27562_1073" ["POLA_1_2003"]="P27562_1072" ["KRAS_4_2002"]="P27562_1075" ["KRAS_3_2002"]="P27562_1074" ["KRAS_2_2002"]="P27562_1070" ["KRAS_1_2002"]="P27562_1052" ["KALM_1_2018"]="P27562_1049" ["JAPA_2_1994"]="P27562_1063" ["JAPA_1_1994"]="P27562_1053" ["ITAL_1_1946"]="P27562_1034" ["GAST_9_1965"]="P27562_1014" ["GAST_8_1965"]="P27562_1011" ["GAST_7_1965"]="P27562_1010" ["GAST_6_1965"]="P27562_1004" ["GAST_5_1943"]="P27562_1015" ["GAST_4_1941"]="P27562_1019" ["GAST_3_1941"]="P27562_1017" ["GAST_2_1941"]="P27562_1009" ["GAST_14_1969"]="P27562_1032" ["GAST_13_1969"]="P27562_1027" ["GAST_12_1969"]="P27562_1012" ["GAST_11_1965"]="P27562_1018" ["GAST_10_1965"]="P27562_1016" ["GAST_1_1941"]="P27562_1003" ["DALA_1_1965"]="P27562_1006" ["CHEH_3_2008"]="P27562_1047" ["CHEH_2_2004"]="P27562_1048" ["CHEH_1_1989"]="P27562_1058" ["BELA_1_2004"]="P27562_1045" ["BAJK_3_2016"]="P27562_1064" ["BAJK_2_2016"]="P27562_1043" ["BAJK_1_2016"]="P27562_1040" ["BAJK_1_2016"]="P27562_1039")


echo 'patient,sample,lane,fastq_1,fastq_2'
for files in "${!groups[@]}"; do
        #echo "$files - ${groups[$files]}";
        read1=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R1')
        read2=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R2')
        echo "$files,$files,L001,$read1,$read2"
done

### Creating groups
#creating IDs correspondence output
awk -F, 'NR!=1 {print}' M.britomartis_sample_sarek.csv | awk 'BEGIN{ORS=" ";} {print "[\"S"NR"_"$2"\"]=\""$1"\""}'

awk -F, 'BEGIN{ORS=" ";} {print "[\""$8"_"$4"\"]=\""$1"\""}' M.britomartis_sample_sarek_.csv

### Test sarek run - succesful

interactive -A naiss2023-5-52 -p node -t 12:00:00


#### Plan for further analysis
<<todo:
1. Try new trimming on 2 files from the test run within snakemake
2. QC, decide how to trim the rest
3. Assemble mtDNA for all the samples (also snakemake?)

todo


module load bioinfo-tools snakemake

#Test data:
wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz
tar --wildcards -xf v5.4.5.tar.gz --strip 1 "*/data"

mkdir workflow
touch workflow/Snakefile


mkdir snakemake-tutorial
curl -L https://api.github.com/repos/snakemake/snakemake-tutorial-data/tarball -o snakemake-tutorial-data.tar.gz
tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"
mamba env create --name snakemake-tutorial --file environment.yaml
conda activate snakemake-tutorial
snakemake -np mapped_reads/A.bam
snakemake --cores 1 mapped_reads/A.bam
nano Snakefile
snakemake --cores 1 mapped_reads/{A,B}.bam

.
├── data
│   ├── genome.fa
│   ├── genome.fa.amb
│   ├── genome.fa.ann
│   ├── genome.fa.bwt
│   ├── genome.fa.fai
│   ├── genome.fa.pac
│   ├── genome.fa.sa
│   └── samples
│       ├── A.fastq
│       ├── B.fastq
│       └── C.fastq
├── environment.yaml
├── mapped_reads
│   ├── A.bam
│   └── B.bam
├── Snakefile
└── snakemake-tutorial-data.tar.gz

#Writing rules!

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


         FASTP	fastp	0.23.2
         FASTQC	fastqc	0.11.9


fastp --in1 MUZAV26_94-L001_1.fastq.gz --in2 MUZAV26_94-L001_2.fastq.gz --out1 MUZAV26_94-L001_1.fastp.fastq.gz --out2 MUZAV26_94-L001_2.fastp.fastq.gz --json MUZAV26_94-L001.fastp.json --html MUZAV26_94-L001.fastp.html --thread 12 --detect_adapter_for_pe --disable_adapter_trimming --split_by_lines 200000000

printf "%s\n" P27562_1039_S39_L001_R1_001.fastq.gz P27562_1039_S39_L001_R2_001.fastq.gz | while read f; do [[ $f =~ ^BAJK_1_2016-L001.* ]] || ln -s $f BAJK_1_2016-L001_$f ; done
    fastqc --quiet --threads 4 BAJK_1_2016-L001*

    cat <<-END_VERSIONS > versions.yml
    "NFCORE_SAREK:SAREK:FASTQC":
        fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS

  [ ! -f  BAJK_1_2016-L001_1.fastq.gz ] && ln -sf P27562_1039_S39_L001_R1_001.fastq.gz BAJK_1_2016-L001_1.fastq.gz
            [ ! -f  BAJK_1_2016-L001_2.fastq.gz ] && ln -sf P27562_1039_S39_L001_R2_001.fastq.gz BAJK_1_2016-L001_2.fastq.gz
            fastp \
                --in1 BAJK_1_2016-L001_1.fastq.gz \
                --in2 BAJK_1_2016-L001_2.fastq.gz \
                --out1 BAJK_1_2016-L001_1.fastp.fastq.gz \
                --out2 BAJK_1_2016-L001_2.fastp.fastq.gz \
                --json BAJK_1_2016-L001.fastp.json \
                --html BAJK_1_2016-L001.fastp.html \
                 \
                 \
                 \
                --thread 12 \
                --detect_adapter_for_pe \
                --disable_adapter_trimming      --split_by_lines 200000000 \
                2> BAJK_1_2016-L001.fastp.log

            cat <<-END_VERSIONS > versions.yml
            "NFCORE_SAREK:SAREK:FASTP":
                fastp: $(fastp --version 2>&1 | sed -e "s/fastp //g")
            END_VERSIONS



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






module load fastp

###Finding file pathes
find /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x -name "*.fastq.gz"  | sort > LR999924.g.vcf.list

fastp \
    --in1 SMAL_2_1967-L001_1.fastp.fastq.gz \
    --in2 SMAL_2_1967-L001_2.fastp.fastq.gz \
    --out1 SMAL_2_1967-L001_1.fastpX.fastq.gz \
    --out2 SMAL_2_1967-L001_2.fastpX.fastq.gz \
    --json SMAL_2_1967-L001.fastpX.json \
    --html SMAL_2_1967-L001.fastpX.html \
    --thread 12 \
    --detect_adapter_for_pe \
    --split_by_lines 200000000

Some chuncks:


    /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/work/c7/ac12632666d326dcae9679d1a5ea38/0001.MUZAV26_94-L001_1.fastp.fastq.gz
    /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/work/c7/ac12632666d326dcae9679d1a5ea38/0001.MUZAV26_94-L001_2.fastp.fastq.gz

mv /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/work/c7/ac12632666d326dcae9679d1a5ea38/0001.MUZAV26_94-L001_1.fastp.fastq.gz SMAL_2_1967-L001_1.fastp.fastq.gz
mv /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/work/c7/ac12632666d326dcae9679d1a5ea38/0001.MUZAV26_94-L001_2.fastp.fastq.gz SMAL_2_1967-L001_2.fastp.fastq.gz


fastp \
    --in1 SMAL_2_1967-L001_1.fastq.gz \
    --in2 SMAL_2_1967-L001_2.fastq.gz \
    --out1 SMAL_2_1967-L001_1.fastp.fastq.gz \
    --out2 SMAL_2_1967-L001_2.fastp.fastq.gz \
    --json SMAL_2_1967-L001.fastp.json \
    --html SMAL_2_1967-L001.fastp.html \
    --thread 12 \
    --reads_to_process 1000 \
    --detect_adapter_for_pe \



ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1055/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1055_S55_L001_R1_001.fastq.gz SMAL_2_1967-L001_1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1055/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1055_S55_L001_R2_001.fastq.gz SMAL_2_1967-L001_2.fastq.gz

/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1055/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1055_S55_L001_R1_001.fastq.gz,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1055/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1055_S55_L001_R2_001.fastq.gz


ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R1_001.fastq.gz GAST_5_1943-L001_1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R2_001.fastq.gz GAST_5_1943-L001_2.fastq.gz

fastp \
    --in1 GAST_5_1943-L001_1.fastq.gz \
    --in2 GAST_5_1943-L001_2.fastq.gz \
    --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
    --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
    --json GAST_5_1943-L001.fastp.json \
    --html GAST_5_1943-L001.fastp.html \
    --thread 12 \
    --reads_to_process 1000 \
    --detect_adapter_for_pe \


    fastp \
        --in1 GAST_5_1943-L001_1.fastq.gz \
        --in2 GAST_5_1943-L001_2.fastq.gz \
        --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
        --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
        --json GAST_5_1943-L001.fastp.json \
        --html GAST_5_1943-L001.fastp.html \
        --thread 12 \
        --reads_to_process 1000 \
        --trim_front1=20 \
        --trim_tail1=20 \
        --detect_adapter_for_pe \
        -p \
        --adapter_fasta adapters.fasta

>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>polyG

fastp \
    --in1 GAST_5_1943-L001_1.fastq.gz \
    --in2 GAST_5_1943-L001_2.fastq.gz \
    --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
    --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
    --json GAST_5_1943-L001.fastp.json \
    --html GAST_5_1943-L001.fastp.html \
    --thread 12 \
    --reads_to_process 1000 \
    --trim_front1=20 \
    --trim_tail1=20 \
    --detect_adapter_for_pe \
    -p \
    --adapter_fasta adapters.fasta


    fastp \
        --in1 GAST_5_1943-L001_1.fastq.gz \
        --in2 GAST_5_1943-L001_2.fastq.gz \
        --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
        --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
        --json GAST_5_1943-L001.fastp.json \
        --html GAST_5_1943-L001.fastp.html \
        --thread 12 \
        --reads_to_process 1000 \
        --trim_front1=20 \
        --trim_tail1=20 \
        --detect_adapter_for_pe \
        -p \
        -r \
        --adapter_fasta adapters.fasta


#### mt genome assembly

module load GetOrganelle/1.7.7.0
module load GetOrganelleDB

 get_organelle_from_reads.py -1 Arabidopsis_simulated.1.fq.gz -2 Arabidopsis_simulated.2.fq.gz -t 1 -o Arabidopsis_simulated.plastome -F embplant_pt -R 10

get_organelle_from_reads.py -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_2.fastp.fastq.gz -t 1 -o SMAL_2_1967.mt -F animal_mt -R 10

-1 ../00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_1.fastp.fastq.gz -2 ../00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_2.fastp.fastq.gz


#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J mtgenome_assembly
#SBATCH --mail-user=daria.shipilina@gmail.com

module load bioinfo-tools GetOrganelle/1.7.7.0
module load GetOrganelleDB

get_organelle_from_reads.py -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_2.fastp.fastq.gz -t 4 -o SMAL_2_1967.mt -F animal_mt -R 10




###### Snakemake for organelle
### working in separate file

mamba env create --name mtgenome --file mapping.yaml
snakemake
snakemake -np mapped_reads/A.bam
snakemake --cores 1 mapped_reads/A.bam
nano Snakefile
snakemake --cores 1 mapped_reads/{A,B}.bam

conda env update --name snakemake-tutorial --file environment.yaml


#### Trying snakemake for mapping and calling

#Added fastp/0.23.2
conda env update --name snakemake-tutorial --file environment.yaml


rule fastqc:
    input:
        "GAST_2_1941-L001_1.fastp.fastq.gz",
        "GAST_2_1941-L001_2.fastp.fastq.gz"
    output:
        html_report = "GAST_2_1941_output.html",
        zip_archive = "GAST_2_1941_output.zip"
    shell:
        """
        fastqc {input} --threads 4 --outdir . &&
        mv *fastqc.html {output.html_report} &&
        mv *fastqc.zip {output.zip_archive}
        """


rule fastp_preprocess:
    input:
	      r1 = "GAST_5_1943-L001_1.fastq.gz",
        r2 = "GAST_5_1943-L001_2.fastq.gz"
    output:
	      r1_out = "GAST_5_1943-L001_1.fastp.fastq.gz",
        r2_out = "GAST_5_1943-L001_2.fastp.fastq.gz",
        json_out = "GAST_5_1943-L001.fastp.json",
        html_out = "GAST_5_1943-L001.fastp.html"
    params:
	      thread = 12,
        reads_to_process = 1000
    shell:
	"""
	fastp \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --out1 {output.r1_out} \
        --out2 {output.r2_out} \
        --json {output.json_out} \
        --html {output.html_out} \
        --thread {params.thread} \
        --reads_to_process {params.reads_to_process} \
        --detect_adapter_for_pe
        """


r1 = expand("../00_Mapping_Calling_sarek/00_Benchmarking4x/{sample}-L001_1.fastp.fastq.gz", sample=samples),
r2 = expand("../00_Mapping_Calling_sarek/00_Benchmarking4x/{sample}-L001_2.fastp.fastq.gz", sample=samples


rule fastqc:
    input:
        "data/samples/GAST_5_1943-L001_1.fastq.gz",
        "data/samples/GAST_5_1943-L001_2.fastq.gz"
    output:
        html_report = "GAST_5_1943_output.html",
        zip_archive = "GAST_5_1943_output.zip"
    shell:
        """
        fastqc {input} --threads 4 --outdir . &&
        mv *fastqc.html {output.html_report} &&
        mv *fastqc.zip {output.zip_archive}
        """

rule fastp_preprocess:
    input:
	r1 = "data/samples/GAST_5_1943-L001_1.fastq.gz",
        r2 = "data/samples/GAST_5_1943-L001_2.fastq.gz"
    output:
	r1_out = "GAST_5_1943-L001_1.fastp.fastq.gz",
        r2_out = "GAST_5_1943-L001_2.fastp.fastq.gz",
        json_out = "GAST_5_1943-L001.fastp.json",
        html_out = "GAST_5_1943-L001.fastp.html"
    params:
	thread = 12,
        reads_to_process = 1000
    shell:
	"""
	fastp \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --out1 {output.r1_out} \
        --out2 {output.r2_out} \
        --json {output.json_out} \
        --html {output.html_out} \
        --thread {params.thread} \
        --reads_to_process {params.reads_to_process} \
        --detect_adapter_for_pe
        """


snakemake -s fastp.snakemake --verbose --debug-dag --use-conda

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming

2023-06-02 09:47:20 jobid=38364534 jobstate=COMPLETED username=daria account=naiss2023-5-52 nodes=r1225 procs=4 partition=core qos=normal jobname=mtgenome_assembly maxmemory_in_GiB=11.3 maxmemory_node=r1225 timelimit=10:00:00 submit_time=2023-06-02T08:59:48 start_time=2023-06-02T08:59:57 end_time=2023-06-02T09:47:20 runtime=00:47:23 margin=09:12:37 queuetime=00:00:09



mamba install -c bioconda fastqc=0.12.1
mamba install -c bioconda fastp=0.23.4

get_organelle_from_reads.py -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming/GAST_5_1943-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming/GAST_5_1943-L001_2.fastp.fastq.gz -t 4 -o GAST_5_1943.mt -F animal_mt -R 10 --max-reads 40E8  --reduce-reads-for-coverage 50

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/02_mtGenome/SMAL_2_1967.mt/seed/animal_mt.initial.fq

animal_mt.initial.fq



get_organelle_from_reads.py -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming/GAST_5_1943-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming/GAST_5_1943-L001_2.fastp.fastq.gz -s /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/02_mtGenome/SMAL_2_1967.mt/animal_mt.K115.complete.graph1.1.path_sequence.fasta -t 4 -o GAST_5_1943.mt -F animal_mt -s /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/02_mtGenome/SMAL_2_1967.mt/animal_mt.K115.complete.graph1.1.path_sequence.fasta

get_organelle_from_reads.py -s SMAL_2_1967.fasta -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming/GAST_5_1943-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming/GAST_5_1943-L001_2.fastp.fastq.gz -o GAST_5_1943.mt -F animal_mt



### Making contemporary set

### Making sample list
declare -A groups=(["КAZA_2_2006"]="P27562_1068" ["КAZA_1_2006"]="P27562_1038" ["КALM_9_1981"]="P27562_1056" ["КALM_10_1983"]="P27562_1057" ["АLTA_2_2015"]="P27562_1069" ["АLTA_1_2015"]="P27562_1051" ["VORO_1_1998"]="P27562_1071" ["VAST_3_2005"]="P27562_1055" ["VAST_2_1999"]="P27562_1066" ["VAST_1_1983"]="P27562_1065" ["URAL_2_1995"]="P27562_1061" ["URAL_1_1991"]="P27562_1062"  ["STOC_6_1965"]="P27562_1035" ["SMAL_6_2013"]="P27562_1044" ["SMAL_5_1998"]="P27562_1037" ["SMAL_4_1996"]="P27562_1054"  ["SLOV_1_1990"]="P27562_1059" ["RUSS_5_2008"]="P27562_1042" ["RUSS_4_2008"]="P27562_1041" ["RUSS_3_1999"]="P27562_1060" ["RUSS_2_1999"]="P27562_1046" ["RUSS_1_1998"]="P27562_1067" ["POLA_2_2003"]="P27562_1073" ["POLA_1_2003"]="P27562_1072" ["KRAS_4_2002"]="P27562_1075" ["KRAS_3_2002"]="P27562_1074" ["KRAS_2_2002"]="P27562_1070" ["KRAS_1_2002"]="P27562_1052" ["KALM_1_2018"]="P27562_1049" ["JAPA_2_1994"]="P27562_1063" ["JAPA_1_1994"]="P27562_1053" ["CHEH_3_2008"]="P27562_1047" ["CHEH_2_2004"]="P27562_1048" ["CHEH_1_1989"]="P27562_1058" ["BELA_1_2004"]="P27562_1045" ["BAJK_3_2016"]="P27562_1064" ["BAJK_2_2016"]="P27562_1043" ["BAJK_4_2016"]="P27562_1040" ["BAJK_1_2016"]="P27562_1039")


for files in "${!groups[@]}"; do
        read1=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R1')
        echo "ln -sf $read1 $files.R1.fastq.gz"
        read2=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R2')
        echo "ln -sf $read2 $files.R2.fastq.gz"

done



ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1039/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1039_S39_L001_R1_001.fastq.gz BAJK_1_2016.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1039/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1039_S39_L001_R2_001.fastq.gz BAJK_1_2016.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1059/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1059_S59_L001_R1_001.fastq.gz SLOV_1_1990.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1059/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1059_S59_L001_R2_001.fastq.gz SLOV_1_1990.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1061/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1061_S61_L001_R1_001.fastq.gz URAL_2_1995.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1061/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1061_S61_L001_R2_001.fastq.gz URAL_2_1995.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1041/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1041_S41_L001_R1_001.fastq.gz RUSS_4_2008.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1041/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1041_S41_L001_R2_001.fastq.gz RUSS_4_2008.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1053/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1053_S53_L001_R1_001.fastq.gz JAPA_1_1994.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1053/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1053_S53_L001_R2_001.fastq.gz JAPA_1_1994.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1051/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1051_S51_L001_R1_001.fastq.gz АLTA_1_2015.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1051/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1051_S51_L001_R2_001.fastq.gz АLTA_1_2015.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1057/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1057_S57_L001_R1_001.fastq.gz КALM_10_1983.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1057/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1057_S57_L001_R2_001.fastq.gz КALM_10_1983.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1058/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1058_S58_L001_R1_001.fastq.gz CHEH_1_1989.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1058/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1058_S58_L001_R2_001.fastq.gz CHEH_1_1989.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1070/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1070_S70_L001_R1_001.fastq.gz KRAS_2_2002.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1070/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1070_S70_L001_R2_001.fastq.gz KRAS_2_2002.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1067/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1067_S67_L001_R1_001.fastq.gz RUSS_1_1998.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1067/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1067_S67_L001_R2_001.fastq.gz RUSS_1_1998.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1062/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1062_S62_L001_R1_001.fastq.gz URAL_1_1991.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1062/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1062_S62_L001_R2_001.fastq.gz URAL_1_1991.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1035/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1035_S35_L001_R1_001.fastq.gz STOC_6_1965.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1035/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1035_S35_L001_R2_001.fastq.gz STOC_6_1965.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1065/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1065_S65_L001_R1_001.fastq.gz VAST_1_1983.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1065/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1065_S65_L001_R2_001.fastq.gz VAST_1_1983.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1054/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1054_S54_L001_R1_001.fastq.gz SMAL_4_1996.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1054/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1054_S54_L001_R2_001.fastq.gz SMAL_4_1996.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1060/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1060_S60_L001_R1_001.fastq.gz RUSS_3_1999.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1060/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1060_S60_L001_R2_001.fastq.gz RUSS_3_1999.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1056/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1056_S56_L001_R1_001.fastq.gz КALM_9_1981.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1056/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1056_S56_L001_R2_001.fastq.gz КALM_9_1981.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1047/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1047_S47_L001_R1_001.fastq.gz CHEH_3_2008.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1047/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1047_S47_L001_R2_001.fastq.gz CHEH_3_2008.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1068/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1068_S68_L001_R1_001.fastq.gz КAZA_2_2006.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1068/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1068_S68_L001_R2_001.fastq.gz КAZA_2_2006.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1043/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1043_S43_L001_R1_001.fastq.gz BAJK_2_2016.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1043/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1043_S43_L001_R2_001.fastq.gz BAJK_2_2016.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1037/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1037_S37_L001_R1_001.fastq.gz SMAL_5_1998.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1037/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1037_S37_L001_R2_001.fastq.gz SMAL_5_1998.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1055/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1055_S55_L001_R1_001.fastq.gz VAST_3_2005.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1055/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1055_S55_L001_R2_001.fastq.gz VAST_3_2005.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1063/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1063_S63_L001_R1_001.fastq.gz JAPA_2_1994.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1063/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1063_S63_L001_R2_001.fastq.gz JAPA_2_1994.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1046/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1046_S46_L001_R1_001.fastq.gz RUSS_2_1999.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1046/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1046_S46_L001_R2_001.fastq.gz RUSS_2_1999.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1074/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1074_S74_L001_R1_001.fastq.gz KRAS_3_2002.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1074/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1074_S74_L001_R2_001.fastq.gz KRAS_3_2002.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1072/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1072_S72_L001_R1_001.fastq.gz POLA_1_2003.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1072/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1072_S72_L001_R2_001.fastq.gz POLA_1_2003.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R1_001.fastq.gz KALM_1_2018.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R2_001.fastq.gz KALM_1_2018.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1038/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1038_S38_L001_R1_001.fastq.gz КAZA_1_2006.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1038/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1038_S38_L001_R2_001.fastq.gz КAZA_1_2006.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1075/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1075_S75_L001_R1_001.fastq.gz KRAS_4_2002.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1075/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1075_S75_L001_R2_001.fastq.gz KRAS_4_2002.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1044/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1044_S44_L001_R1_001.fastq.gz SMAL_6_2013.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1044/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1044_S44_L001_R2_001.fastq.gz SMAL_6_2013.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1071/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1071_S71_L001_R1_001.fastq.gz VORO_1_1998.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1071/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1071_S71_L001_R2_001.fastq.gz VORO_1_1998.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1064/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1064_S64_L001_R1_001.fastq.gz BAJK_3_2016.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1064/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1064_S64_L001_R2_001.fastq.gz BAJK_3_2016.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1069/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1069_S69_L001_R1_001.fastq.gz АLTA_2_2015.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1069/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1069_S69_L001_R2_001.fastq.gz АLTA_2_2015.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1040/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1040_S40_L001_R1_001.fastq.gz BAJK_4_2016.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1040/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1040_S40_L001_R2_001.fastq.gz BAJK_4_2016.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1048/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1048_S48_L001_R1_001.fastq.gz CHEH_2_2004.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1048/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1048_S48_L001_R2_001.fastq.gz CHEH_2_2004.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1042/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1042_S42_L001_R1_001.fastq.gz RUSS_5_2008.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1042/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1042_S42_L001_R2_001.fastq.gz RUSS_5_2008.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1066/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1066_S66_L001_R1_001.fastq.gz VAST_2_1999.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1066/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1066_S66_L001_R2_001.fastq.gz VAST_2_1999.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1073/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1073_S73_L001_R1_001.fastq.gz POLA_2_2003.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1073/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1073_S73_L001_R2_001.fastq.gz POLA_2_2003.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1052/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1052_S52_L001_R1_001.fastq.gz KRAS_1_2002.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1052/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1052_S52_L001_R2_001.fastq.gz KRAS_1_2002.R2.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1045/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1045_S45_L001_R1_001.fastq.gz BELA_1_2004.R1.fastq.gz
ln -sf /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1045/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1045_S45_L001_R2_001.fastq.gz BELA_1_2004.R2.fastq.gz


#!/bin/bash
#SBATCH -A project_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J snakejob
#SBATCH --error=snakemake.err
#SBATCH --output=snakemake.out

module load bioinfo-tools
module load snakemake

snakemake --cluster "sbatch -A {cluster.project} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -J {cluster.jobname} -o {cluster.output} -e {cluster.error}" --snakefile Snakefile --jobs 999 --use-conda --cluster-config cluster.json


printf "%s\n" P27562_1004_S4_L001_R1_001.fastq.gz P27562_1004_S4_L001_R2_001.fastq.gz | while read f; do [[ $f =~ ^GAST_6_1965-L001.* ]] || ln -s $f GAST_6_1965-L001_$f ; done
    fastqc --quiet --threads 4 GAST_6_1965-L001*

    cat <<-END_VERSIONS > versions.yml
    "NFCORE_SAREK:SAREK:FASTQC":
        fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS



### Running trimming
module load bioinfo-tools
module load snakemake
conda activate snakemake-tutorial

snakemake -s snakefile  # main run

#Adding more inds

#Preparing cluster command
(see cluster.json)

snakemake --cluster "sbatch -A naiss2023-5-52 -p core -n 1 -t 01:00:00 -J cluster.jobname -o cluster.output -e cluster.error" --snakefile snakefile --jobs 999 --cluster-config cluster.json

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J mtgenome_assembly
#SBATCH --mail-user=daria.shipilina@gmail.com

--use-conda


snakemake --snakefile snakefile -p -j 64  --cluster "sbatch -A naiss2023-5-52 -p {cluster.partition} -n {cluster.n} -C {cluster.C} -t {cluster.time} -e {cluster.error} -o {cluster.output}" --cluster-config cluster.json



#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -J mapping
#SBATCH --mail-user=daria.shipilina@gmail.com

module load bioinfo-tools
module load bwa
module load samtools

mkdir -p results/bam/
mkdir -p logs/bwa/
mkdir -p tmp/bwa/POLA_2_2003
bwa mem -R '@RG\tID:POLA_2_2003\tSM:POLA_2_2003\tLB:POLA_2_2003\tPL:Illumina' -t 20 -M /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna ../results/fastp/POLA_2_2003_1.fastp.fastq.gz ../results/fastp/POLA_2_2003_2.fastp.fastq.gz         |samtools sort -m 6G -@20 -T tmp/bwa/POLA_2_2003 - >results/bam/POLA_2_2003.bam.tmp 2>logs/bwa/POLA_2_2003.log         && mv results/bam/POLA_2_2003.bam.tmp results/bam/POLA_2_2003.bam


module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Concervation/00_Mapping_Calling_sarek/00_Benchmarking4x

nextflow run nf-core/sarek --input sample_sheet.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna.gz --skip_tools baserecalibrator --outdir ./results



# Benchmarking step somehow rewrtitten

### Mapping of all the samples

#!/bin/bash -l
module load bioinfo-tools Nextflow nf-core nf-core-pipelines

#Confiruging Nextflow in the working directory
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll
NXF_OPTS='-Xms1g -Xmx4g'

#Providing full path to new data
#gtf="/crex/proj/uppstore2017185/b2014034/nobackup/Dasha/VanessaRNAseq/genome_data/vcard.gtf"
fasta="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna"
nextflow run nf-core/sarek --input contemporary_samples.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --skip_tools baserecalibrator --outdir ./results

[magical_gautier]

#Update sample List
#Update reference

declare -A groups=(["KAZA_2_2006"]="P27562_1068" ["KAZA_1_2006"]="P27562_1038" ["KALM_9_1981"]="P27562_1056" ["KALM_10_1983"]="P27562_1057" ["ALTA_2_2015"]="P27562_1069" ["ALTA_1_2015"]="P27562_1051" ["VORO_1_1998"]="P27562_1071" ["VAST_3_2005"]="P27562_1055" ["VAST_2_1999"]="P27562_1066" ["VAST_1_1983"]="P27562_1065" ["URAL_2_1995"]="P27562_1061" ["URAL_1_1991"]="P27562_1062"  ["STOC_6_1965"]="P27562_1035" ["SMAL_6_2013"]="P27562_1044" ["SMAL_5_1998"]="P27562_1037" ["SMAL_4_1996"]="P27562_1054"  ["SLOV_1_1990"]="P27562_1059" ["RUSS_5_2008"]="P27562_1042" ["RUSS_4_2008"]="P27562_1041" ["RUSS_3_1999"]="P27562_1060" ["RUSS_2_1999"]="P27562_1046" ["RUSS_1_1998"]="P27562_1067" ["POLA_2_2003"]="P27562_1073" ["POLA_1_2003"]="P27562_1072" ["KRAS_4_2002"]="P27562_1075" ["KRAS_3_2002"]="P27562_1074" ["KRAS_2_2002"]="P27562_1070" ["KRAS_1_2002"]="P27562_1052" ["KALM_1_2018"]="P27562_1049" ["JAPA_2_1994"]="P27562_1063" ["JAPA_1_1994"]="P27562_1053" ["CHEH_3_2008"]="P27562_1047" ["CHEH_2_2004"]="P27562_1048" ["CHEH_1_1989"]="P27562_1058" ["BELA_1_2004"]="P27562_1045" ["BAJK_3_2016"]="P27562_1064" ["BAJK_2_2016"]="P27562_1043" ["BAJK_4_2016"]="P27562_1040" ["BAJK_1_2016"]="P27562_1039")


echo 'patient,sample,lane,fastq_1,fastq_2'
for files in "${!groups[@]}"; do
        #echo "$files - ${groups[$files]}";
        read1=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R1')
        read2=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R2')
        echo "$files,$files,L001,$read1,$read2"
done


contemporary_samples.txt

#Update reference

awk '{print $1}' GCA_905220545.2_ilMelAtha1.2_genomic.fna.fai >  GCA_905220545.2_ilMelAtha1.2_genomic.fna.chroms.names
module load seqtk
seqtk subseq GCA_905220545.2_ilMelAtha1.2_genomic.fna GCA_905220545.2_ilMelAtha1.2_genomic.fna.chroms.names > GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna


Error executing process > 'NFCORE_SAREK:SAREK:FASTQC (АLTA_1_2015-L001)'

Caused by:
  Missing output file(s) `*.html` expected by process `NFCORE_SAREK:SAREK:FASTQC (АLTA_1_2015-L001)`

Command executed:

  printf "%s\n" P27562_1051_S51_L001_R1_001.fastq.gz P27562_1051_S51_L001_R2_001.fastq.gz | while read f; do [[ $f =~ ^АLTA_1_2015-L001.* ]] || ln -s $f АLTA_1_2015-L001_$f ; done
  fastqc --quiet --threads 4 АLTA_1_2015-L001*

  cat <<-END_VERSIONS > versions.yml
  "NFCORE_SAREK:SAREK:FASTQC":
      fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
  END_VERSIONS

Command exit status:
  0

Command output:
  (empty)

Command error:
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  INFO:    Environment variable SINGULARITYENV_SNIC_TMP is set, but APPTAINERENV_SNIC_TMP is preferred
  WARNING: Skipping mount /var/apptainer/mnt/session/etc/resolv.conf [files]: /etc/resolv.conf doesn't exist in container
  Skipping '??LTA_1_2015-L001_P27562_1051_S51_L001_R1_001.fastq.gz' which didn't exist, or couldn't be read
  Skipping '??LTA_1_2015-L001_P27562_1051_S51_L001_R2_001.fastq.gz' which didn't exist, or couldn't be read

Work dir:
  /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/work/89/dd505f764287b08302fb689a1e24bf

Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

#restart
nextflow run nf-core/sarek -resume magical_gautier --input contemporary_samples.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --skip_tools baserecalibrator --outdir ./results

#lets call snps with freebayes
mkdir 04_CallingContemporary
module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Concervation/00_Mapping_Calling_sarek/04_CallingContemporary


declare -A groups=(["KAZA_2_2006"]="P27562_1068" ["KAZA_1_2006"]="P27562_1038" ["KALM_9_1981"]="P27562_1056" ["KALM_10_1983"]="P27562_1057" ["ALTA_2_2015"]="P27562_1069" ["ALTA_1_2015"]="P27562_1051" ["VORO_1_1998"]="P27562_1071" ["VAST_3_2005"]="P27562_1055" ["VAST_2_1999"]="P27562_1066" ["VAST_1_1983"]="P27562_1065" ["URAL_2_1995"]="P27562_1061" ["URAL_1_1991"]="P27562_1062"  ["STOC_6_1965"]="P27562_1035" ["SMAL_6_2013"]="P27562_1044" ["SMAL_5_1998"]="P27562_1037" ["SMAL_4_1996"]="P27562_1054"  ["SLOV_1_1990"]="P27562_1059" ["RUSS_5_2008"]="P27562_1042" ["RUSS_4_2008"]="P27562_1041" ["RUSS_3_1999"]="P27562_1060" ["RUSS_2_1999"]="P27562_1046" ["RUSS_1_1998"]="P27562_1067" ["POLA_2_2003"]="P27562_1073" ["POLA_1_2003"]="P27562_1072" ["KRAS_4_2002"]="P27562_1075" ["KRAS_3_2002"]="P27562_1074" ["KRAS_2_2002"]="P27562_1070" ["KRAS_1_2002"]="P27562_1052" ["KALM_1_2018"]="P27562_1049" ["JAPA_2_1994"]="P27562_1063" ["JAPA_1_1994"]="P27562_1053" ["CHEH_3_2008"]="P27562_1047" ["CHEH_2_2004"]="P27562_1048" ["CHEH_1_1989"]="P27562_1058" ["BELA_1_2004"]="P27562_1045" ["BAJK_3_2016"]="P27562_1064" ["BAJK_2_2016"]="P27562_1043" ["BAJK_4_2016"]="P27562_1040" ["BAJK_1_2016"]="P27562_1039")

echo 'patient,sample,cram,crai'
for files in "${!groups[@]}"; do
        echo "$files,$files,../01_MappingAll/results/preprocessing/markduplicates/$files/$files.md.cram,../01_MappingAll/results/preprocessing/markduplicates/$files/$files.md.cram.crai"
done

patient,sample,cram,crai
BAJK_1_2016,BAJK_1_2016,../01_MappingAll/results/preprocessing/markduplicates/BAJK_1_2016/BAJK_1_2016.md.cram,../01_MappingAll/results/preprocessing/markduplicates/BAJK_1_2016/BAJK_1_2016.md.cram.crai
KAZA_1_2006,KAZA_1_2006,../01_MappingAll/results/preprocessing/markduplicates/KAZA_1_2006/KAZA_1_2006.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KAZA_1_2006/KAZA_1_2006.md.cram.crai
SLOV_1_1990,SLOV_1_1990,../01_MappingAll/results/preprocessing/markduplicates/SLOV_1_1990/SLOV_1_1990.md.cram,../01_MappingAll/results/preprocessing/markduplicates/SLOV_1_1990/SLOV_1_1990.md.cram.crai
URAL_2_1995,URAL_2_1995,../01_MappingAll/results/preprocessing/markduplicates/URAL_2_1995/URAL_2_1995.md.cram,../01_MappingAll/results/preprocessing/markduplicates/URAL_2_1995/URAL_2_1995.md.cram.crai
RUSS_4_2008,RUSS_4_2008,../01_MappingAll/results/preprocessing/markduplicates/RUSS_4_2008/RUSS_4_2008.md.cram,../01_MappingAll/results/preprocessing/markduplicates/RUSS_4_2008/RUSS_4_2008.md.cram.crai
JAPA_1_1994,JAPA_1_1994,../01_MappingAll/results/preprocessing/markduplicates/JAPA_1_1994/JAPA_1_1994.md.cram,../01_MappingAll/results/preprocessing/markduplicates/JAPA_1_1994/JAPA_1_1994.md.cram.crai
CHEH_1_1989,CHEH_1_1989,../01_MappingAll/results/preprocessing/markduplicates/CHEH_1_1989/CHEH_1_1989.md.cram,../01_MappingAll/results/preprocessing/markduplicates/CHEH_1_1989/CHEH_1_1989.md.cram.crai
KRAS_2_2002,KRAS_2_2002,../01_MappingAll/results/preprocessing/markduplicates/KRAS_2_2002/KRAS_2_2002.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KRAS_2_2002/KRAS_2_2002.md.cram.crai
RUSS_1_1998,RUSS_1_1998,../01_MappingAll/results/preprocessing/markduplicates/RUSS_1_1998/RUSS_1_1998.md.cram,../01_MappingAll/results/preprocessing/markduplicates/RUSS_1_1998/RUSS_1_1998.md.cram.crai
URAL_1_1991,URAL_1_1991,../01_MappingAll/results/preprocessing/markduplicates/URAL_1_1991/URAL_1_1991.md.cram,../01_MappingAll/results/preprocessing/markduplicates/URAL_1_1991/URAL_1_1991.md.cram.crai
STOC_6_1965,STOC_6_1965,../01_MappingAll/results/preprocessing/markduplicates/STOC_6_1965/STOC_6_1965.md.cram,../01_MappingAll/results/preprocessing/markduplicates/STOC_6_1965/STOC_6_1965.md.cram.crai
VAST_1_1983,VAST_1_1983,../01_MappingAll/results/preprocessing/markduplicates/VAST_1_1983/VAST_1_1983.md.cram,../01_MappingAll/results/preprocessing/markduplicates/VAST_1_1983/VAST_1_1983.md.cram.crai
SMAL_4_1996,SMAL_4_1996,../01_MappingAll/results/preprocessing/markduplicates/SMAL_4_1996/SMAL_4_1996.md.cram,../01_MappingAll/results/preprocessing/markduplicates/SMAL_4_1996/SMAL_4_1996.md.cram.crai
RUSS_3_1999,RUSS_3_1999,../01_MappingAll/results/preprocessing/markduplicates/RUSS_3_1999/RUSS_3_1999.md.cram,../01_MappingAll/results/preprocessing/markduplicates/RUSS_3_1999/RUSS_3_1999.md.cram.crai
CHEH_3_2008,CHEH_3_2008,../01_MappingAll/results/preprocessing/markduplicates/CHEH_3_2008/CHEH_3_2008.md.cram,../01_MappingAll/results/preprocessing/markduplicates/CHEH_3_2008/CHEH_3_2008.md.cram.crai
BAJK_2_2016,BAJK_2_2016,../01_MappingAll/results/preprocessing/markduplicates/BAJK_2_2016/BAJK_2_2016.md.cram,../01_MappingAll/results/preprocessing/markduplicates/BAJK_2_2016/BAJK_2_2016.md.cram.crai
SMAL_5_1998,SMAL_5_1998,../01_MappingAll/results/preprocessing/markduplicates/SMAL_5_1998/SMAL_5_1998.md.cram,../01_MappingAll/results/preprocessing/markduplicates/SMAL_5_1998/SMAL_5_1998.md.cram.crai
VAST_3_2005,VAST_3_2005,../01_MappingAll/results/preprocessing/markduplicates/VAST_3_2005/VAST_3_2005.md.cram,../01_MappingAll/results/preprocessing/markduplicates/VAST_3_2005/VAST_3_2005.md.cram.crai
ALTA_2_2015,ALTA_2_2015,../01_MappingAll/results/preprocessing/markduplicates/ALTA_2_2015/ALTA_2_2015.md.cram,../01_MappingAll/results/preprocessing/markduplicates/ALTA_2_2015/ALTA_2_2015.md.cram.crai
JAPA_2_1994,JAPA_2_1994,../01_MappingAll/results/preprocessing/markduplicates/JAPA_2_1994/JAPA_2_1994.md.cram,../01_MappingAll/results/preprocessing/markduplicates/JAPA_2_1994/JAPA_2_1994.md.cram.crai
KALM_9_1981,KALM_9_1981,../01_MappingAll/results/preprocessing/markduplicates/KALM_9_1981/KALM_9_1981.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KALM_9_1981/KALM_9_1981.md.cram.crai
RUSS_2_1999,RUSS_2_1999,../01_MappingAll/results/preprocessing/markduplicates/RUSS_2_1999/RUSS_2_1999.md.cram,../01_MappingAll/results/preprocessing/markduplicates/RUSS_2_1999/RUSS_2_1999.md.cram.crai
KRAS_3_2002,KRAS_3_2002,../01_MappingAll/results/preprocessing/markduplicates/KRAS_3_2002/KRAS_3_2002.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KRAS_3_2002/KRAS_3_2002.md.cram.crai
POLA_1_2003,POLA_1_2003,../01_MappingAll/results/preprocessing/markduplicates/POLA_1_2003/POLA_1_2003.md.cram,../01_MappingAll/results/preprocessing/markduplicates/POLA_1_2003/POLA_1_2003.md.cram.crai
KALM_1_2018,KALM_1_2018,../01_MappingAll/results/preprocessing/markduplicates/KALM_1_2018/KALM_1_2018.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KALM_1_2018/KALM_1_2018.md.cram.crai
KRAS_4_2002,KRAS_4_2002,../01_MappingAll/results/preprocessing/markduplicates/KRAS_4_2002/KRAS_4_2002.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KRAS_4_2002/KRAS_4_2002.md.cram.crai
SMAL_6_2013,SMAL_6_2013,../01_MappingAll/results/preprocessing/markduplicates/SMAL_6_2013/SMAL_6_2013.md.cram,../01_MappingAll/results/preprocessing/markduplicates/SMAL_6_2013/SMAL_6_2013.md.cram.crai
VORO_1_1998,VORO_1_1998,../01_MappingAll/results/preprocessing/markduplicates/VORO_1_1998/VORO_1_1998.md.cram,../01_MappingAll/results/preprocessing/markduplicates/VORO_1_1998/VORO_1_1998.md.cram.crai
KALM_10_1983,KALM_10_1983,../01_MappingAll/results/preprocessing/markduplicates/KALM_10_1983/KALM_10_1983.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KALM_10_1983/KALM_10_1983.md.cram.crai
BAJK_3_2016,BAJK_3_2016,../01_MappingAll/results/preprocessing/markduplicates/BAJK_3_2016/BAJK_3_2016.md.cram,../01_MappingAll/results/preprocessing/markduplicates/BAJK_3_2016/BAJK_3_2016.md.cram.crai
BAJK_4_2016,BAJK_4_2016,../01_MappingAll/results/preprocessing/markduplicates/BAJK_4_2016/BAJK_4_2016.md.cram,../01_MappingAll/results/preprocessing/markduplicates/BAJK_4_2016/BAJK_4_2016.md.cram.crai
ALTA_1_2015,ALTA_1_2015,../01_MappingAll/results/preprocessing/markduplicates/ALTA_1_2015/ALTA_1_2015.md.cram,../01_MappingAll/results/preprocessing/markduplicates/ALTA_1_2015/ALTA_1_2015.md.cram.crai
CHEH_2_2004,CHEH_2_2004,../01_MappingAll/results/preprocessing/markduplicates/CHEH_2_2004/CHEH_2_2004.md.cram,../01_MappingAll/results/preprocessing/markduplicates/CHEH_2_2004/CHEH_2_2004.md.cram.crai
RUSS_5_2008,RUSS_5_2008,../01_MappingAll/results/preprocessing/markduplicates/RUSS_5_2008/RUSS_5_2008.md.cram,../01_MappingAll/results/preprocessing/markduplicates/RUSS_5_2008/RUSS_5_2008.md.cram.crai
VAST_2_1999,VAST_2_1999,../01_MappingAll/results/preprocessing/markduplicates/VAST_2_1999/VAST_2_1999.md.cram,../01_MappingAll/results/preprocessing/markduplicates/VAST_2_1999/VAST_2_1999.md.cram.crai
POLA_2_2003,POLA_2_2003,../01_MappingAll/results/preprocessing/markduplicates/POLA_2_2003/POLA_2_2003.md.cram,../01_MappingAll/results/preprocessing/markduplicates/POLA_2_2003/POLA_2_2003.md.cram.crai
KRAS_1_2002,KRAS_1_2002,../01_MappingAll/results/preprocessing/markduplicates/KRAS_1_2002/KRAS_1_2002.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KRAS_1_2002/KRAS_1_2002.md.cram.crai
BELA_1_2004,BELA_1_2004,../01_MappingAll/results/preprocessing/markduplicates/BELA_1_2004/BELA_1_2004.md.cram,../01_MappingAll/results/preprocessing/markduplicates/BELA_1_2004/BELA_1_2004.md.cram.crai
KAZA_2_2006,KAZA_2_2006,../01_MappingAll/results/preprocessing/markduplicates/KAZA_2_2006/KAZA_2_2006.md.cram,../01_MappingAll/results/preprocessing/markduplicates/KAZA_2_2006/KAZA_2_2006.md.cram.crai


nextflow run nf-core/sarek -r 3.1.2 --input contemporary_calling_input.csv -profile uppmax --project naiss2023-5-52 --fasta  /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --step variant_calling --tools freebayes --outdir ./results


#Lets map low q
mkdir 03_MappingHistorical


### Making sample list

#!/bin/bash



declare -A groups=(["КALM_8_1969"]="P27562_1029" ["КALM_7_1961"]="P27562_1028" ["КALM_6_1961"]="P27562_1025" ["КALM_5_1961"]="P27562_1021" ["КALM_4_1957"]="P27562_1024" ["КALM_3_1955"]="P27562_1031" ["КALM_2_1955"]="P27562_1030" ["UPPS_3_1969"]="P27562_1022" ["UPPS_2_1958"]="P27562_1023" ["UPPS_1_1951"]="P27562_1026" ["STOC_5_1965"]="P27562_1013" ["STOC_4_1965"]="P27562_1008" ["STOC_3_1965"]="P27562_1007" ["STOC_2_1965"]="P27562_1002" ["STOC_1_1965"]="P27562_1001" ["SMAL_5_1998"]="P27562_1037" ["SMAL_4_1996"]="P27562_1054" ["SMAL_3_1967"]="P27562_1033" ["SMAL_2_1967"]="P27562_1020" ["SMAL_1_1967"]="P27562_1005" ["ITAL_1_1946"]="P27562_1034" ["GAST_9_1965"]="P27562_1014" ["GAST_8_1965"]="P27562_1011" ["GAST_7_1965"]="P27562_1010" ["GAST_6_1965"]="P27562_1004" ["GAST_5_1943"]="P27562_1015" ["GAST_4_1941"]="P27562_1019" ["GAST_3_1941"]="P27562_1017" ["GAST_2_1941"]="P27562_1009" ["GAST_14_1969"]="P27562_1032" ["GAST_13_1969"]="P27562_1027" ["GAST_12_1969"]="P27562_1012" ["GAST_11_1965"]="P27562_1018" ["GAST_10_1965"]="P27562_1016" ["GAST_1_1941"]="P27562_1003" ["DALA_1_1965"]="P27562_1006")


echo 'patient,sample,lane,fastq_1,fastq_2'
for files in "${!groups[@]}"; do
        #echo "$files - ${groups[$files]}";
        read1=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R1')
        read2=$(cat /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/fasta_pathes_britomartis.txt | grep ${groups[$files]} | grep 'R2')
        echo "$files,$files,L001,$read1,$read2"
done


nextflow run nf-core/sarek -r 3.1.2 --input historical_sample_list.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --skip_tools baserecalibrator --outdir ./results



module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Concervation/00_Mapping_Calling_sarek/03_MappingHistorical


historical_sample_list.csv

-resume hopeful_monod

nextflow run nf-core/sarek -r 3.1.2 -resume hopeful_monod --input historical_sample_list.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --skip_tools baserecalibrator --outdir ./results


#lets call snps with freebayes

module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical


  declare -A groups=(["КALM_8_1969"]="P27562_1029" ["КALM_7_1961"]="P27562_1028" ["КALM_6_1961"]="P27562_1025" ["КALM_5_1961"]="P27562_1021" ["КALM_4_1957"]="P27562_1024" ["КALM_3_1955"]="P27562_1031" ["КALM_2_1955"]="P27562_1030" ["UPPS_3_1969"]="P27562_1022" ["UPPS_2_1958"]="P27562_1023" ["UPPS_1_1951"]="P27562_1026" ["STOC_5_1965"]="P27562_1013" ["STOC_4_1965"]="P27562_1008" ["STOC_3_1965"]="P27562_1007" ["STOC_2_1965"]="P27562_1002" ["STOC_1_1965"]="P27562_1001" ["SMAL_5_1998"]="P27562_1037" ["SMAL_4_1996"]="P27562_1054" ["SMAL_3_1967"]="P27562_1033" ["SMAL_2_1967"]="P27562_1020" ["SMAL_1_1967"]="P27562_1005" ["ITAL_1_1946"]="P27562_1034" ["GAST_9_1965"]="P27562_1014" ["GAST_8_1965"]="P27562_1011" ["GAST_7_1965"]="P27562_1010" ["GAST_6_1965"]="P27562_1004" ["GAST_5_1943"]="P27562_1015" ["GAST_4_1941"]="P27562_1019" ["GAST_3_1941"]="P27562_1017" ["GAST_2_1941"]="P27562_1009" ["GAST_14_1969"]="P27562_1032" ["GAST_13_1969"]="P27562_1027" ["GAST_12_1969"]="P27562_1012" ["GAST_11_1965"]="P27562_1018" ["GAST_10_1965"]="P27562_1016" ["GAST_1_1941"]="P27562_1003" ["DALA_1_1965"]="P27562_1006")

  echo 'patient,sample,cram,crai'
  for files in "${!groups[@]}"; do
          echo "$files,$files,../03_MappingHistorical/results/preprocessing/markduplicates/$files/$files.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/$files/$files.md.cram.crai"
  done

nextflow pull nf-core/sarek -r 3.1.2
nextflow run nf-core/sarek -r 3.1.2 --input historical_calling_input.csv -profile uppmax --project naiss2023-5-52 --fasta  /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --step variant_calling --tools freebayes --outdir ./results




patient,sample,cram,crai
GAST_6_1965,GAST_6_1965,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_6_1965/GAST_6_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_6_1965/GAST_6_1965.md.cram.crai
STOC_4_1965,STOC_4_1965,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_4_1965/STOC_4_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_4_1965/STOC_4_1965.md.cram.crai
GAST_14_1969,GAST_14_1969,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_14_1969/GAST_14_1969.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_14_1969/GAST_14_1969.md.cram.crai
UPPS_1_1951,UPPS_1_1951,../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_1_1951/UPPS_1_1951.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_1_1951/UPPS_1_1951.md.cram.crai
KALM_3_1955,KALM_3_1955,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_3_1955/KALM_3_1955.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_3_1955/KALM_3_1955.md.cram.crai
GAST_4_1941,GAST_4_1941,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_4_1941/GAST_4_1941.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_4_1941/GAST_4_1941.md.cram.crai
KALM_5_1961,KALM_5_1961,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_5_1961/KALM_5_1961.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_5_1961/KALM_5_1961.md.cram.crai
GAST_2_1941,GAST_2_1941,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_2_1941/GAST_2_1941.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_2_1941/GAST_2_1941.md.cram.crai
SMAL_4_1996,SMAL_4_1996,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_4_1996/SMAL_4_1996.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_4_1996/SMAL_4_1996.md.cram.crai
UPPS_3_1969,UPPS_3_1969,../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_3_1969/UPPS_3_1969.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_3_1969/UPPS_3_1969.md.cram.crai
KALM_7_1961,KALM_7_1961,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_7_1961/KALM_7_1961.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_7_1961/KALM_7_1961.md.cram.crai
GAST_10_1965,GAST_10_1965,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_10_1965/GAST_10_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_10_1965/GAST_10_1965.md.cram.crai
SMAL_2_1967,SMAL_2_1967,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_2_1967/SMAL_2_1967.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_2_1967/SMAL_2_1967.md.cram.crai
STOC_2_1965,STOC_2_1965,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_2_1965/STOC_2_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_2_1965/STOC_2_1965.md.cram.crai
GAST_7_1965,GAST_7_1965,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_7_1965/GAST_7_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_7_1965/GAST_7_1965.md.cram.crai
SMAL_5_1998,SMAL_5_1998,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_5_1998/SMAL_5_1998.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_5_1998/SMAL_5_1998.md.cram.crai
GAST_11_1965,GAST_11_1965,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_11_1965/GAST_11_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_11_1965/GAST_11_1965.md.cram.crai
GAST_12_1969,GAST_12_1969,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_12_1969/GAST_12_1969.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_12_1969/GAST_12_1969.md.cram.crai
STOC_3_1965,STOC_3_1965,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_3_1965/STOC_3_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_3_1965/STOC_3_1965.md.cram.crai
DALA_1_1965,DALA_1_1965,../03_MappingHistorical/results/preprocessing/markduplicates/DALA_1_1965/DALA_1_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/DALA_1_1965/DALA_1_1965.md.cram.crai
SMAL_3_1967,SMAL_3_1967,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_3_1967/SMAL_3_1967.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_3_1967/SMAL_3_1967.md.cram.crai
STOC_5_1965,STOC_5_1965,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_5_1965/STOC_5_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_5_1965/STOC_5_1965.md.cram.crai
KALM_6_1961,KALM_6_1961,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_6_1961/KALM_6_1961.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_6_1961/KALM_6_1961.md.cram.crai
GAST_13_1969,GAST_13_1969,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_13_1969/GAST_13_1969.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_13_1969/GAST_13_1969.md.cram.crai
GAST_1_1941,GAST_1_1941,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_1_1941/GAST_1_1941.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_1_1941/GAST_1_1941.md.cram.crai
KALM_8_1969,KALM_8_1969,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_8_1969/KALM_8_1969.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_8_1969/KALM_8_1969.md.cram.crai
GAST_8_1965,GAST_8_1965,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_8_1965/GAST_8_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_8_1965/GAST_8_1965.md.cram.crai
KALM_2_1955,KALM_2_1955,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_2_1955/KALM_2_1955.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_2_1955/KALM_2_1955.md.cram.crai
GAST_9_1965,GAST_9_1965,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_9_1965/GAST_9_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_9_1965/GAST_9_1965.md.cram.crai
SMAL_1_1967,SMAL_1_1967,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_1_1967/SMAL_1_1967.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_1_1967/SMAL_1_1967.md.cram.crai
STOC_1_1965,STOC_1_1965,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_1_1965/STOC_1_1965.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/STOC_1_1965/STOC_1_1965.md.cram.crai
KALM_4_1957,KALM_4_1957,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_4_1957/KALM_4_1957.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/KALM_4_1957/KALM_4_1957.md.cram.crai
GAST_3_1941,GAST_3_1941,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_3_1941/GAST_3_1941.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_3_1941/GAST_3_1941.md.cram.crai
ITAL_1_1946,ITAL_1_1946,../03_MappingHistorical/results/preprocessing/markduplicates/ITAL_1_1946/ITAL_1_1946.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/ITAL_1_1946/ITAL_1_1946.md.cram.crai
GAST_5_1943,GAST_5_1943,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_5_1943/GAST_5_1943.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/GAST_5_1943/GAST_5_1943.md.cram.crai
UPPS_2_1958,UPPS_2_1958,../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_2_1958/UPPS_2_1958.md.cram,../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_2_1958/UPPS_2_1958.md.cram.crai


### Let's join Freebayes
### Let's get a list of files to Merge
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes
find . -name "*gz" | awk '{for (i=1; i<=NF; i++) a[NR,i] = $i} NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j];} print str}}'


#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH -J cont_freebayes
#SBATCH --output=cont_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

date
module load bioinfo-tools bcftools
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes
bcftools merge --threads 4 ./SMAL_4_1996/SMAL_4_1996.freebayes.vcf.gz ./KAZA_1_2006/KAZA_1_2006.freebayes.vcf.gz ./POLA_2_2003/POLA_2_2003.freebayes.vcf.gz ./KRAS_4_2002/KRAS_4_2002.freebayes.vcf.gz ./RUSS_3_1999/RUSS_3_1999.freebayes.vcf.gz ./KALM_9_1981/KALM_9_1981.freebayes.vcf.gz ./RUSS_4_2008/RUSS_4_2008.freebayes.vcf.gz ./RUSS_1_1998/RUSS_1_1998.freebayes.vcf.gz ./BELA_1_2004/BELA_1_2004.freebayes.vcf.gz ./RUSS_5_2008/RUSS_5_2008.freebayes.vcf.gz ./KAZA_2_2006/KAZA_2_2006.freebayes.vcf.gz ./URAL_2_1995/URAL_2_1995.freebayes.vcf.gz ./URAL_1_1991/URAL_1_1991.freebayes.vcf.gz ./KALM_1_2018/KALM_1_2018.freebayes.vcf.gz ./KALM_10_1983/KALM_10_1983.freebayes.vcf.gz ./VAST_3_2005/VAST_3_2005.freebayes.vcf.gz ./CHEH_1_1989/CHEH_1_1989.freebayes.vcf.gz ./BAJK_3_2016/BAJK_3_2016.freebayes.vcf.gz ./JAPA_2_1994/JAPA_2_1994.freebayes.vcf.gz ./VORO_1_1998/VORO_1_1998.freebayes.vcf.gz ./SMAL_5_1998/SMAL_5_1998.freebayes.vcf.gz ./RUSS_2_1999/RUSS_2_1999.freebayes.vcf.gz ./KRAS_2_2002/KRAS_2_2002.freebayes.vcf.gz ./ALTA_2_2015/ALTA_2_2015.freebayes.vcf.gz ./SLOV_1_1990/SLOV_1_1990.freebayes.vcf.gz ./POLA_1_2003/POLA_1_2003.freebayes.vcf.gz ./JAPA_1_1994/JAPA_1_1994.freebayes.vcf.gz ./ALTA_1_2015/ALTA_1_2015.freebayes.vcf.gz ./CHEH_3_2008/CHEH_3_2008.freebayes.vcf.gz ./KRAS_1_2002/KRAS_1_2002.freebayes.vcf.gz ./SMAL_6_2013/SMAL_6_2013.freebayes.vcf.gz ./STOC_6_1965/STOC_6_1965.freebayes.vcf.gz ./CHEH_2_2004/CHEH_2_2004.freebayes.vcf.gz ./VAST_2_1999/VAST_2_1999.freebayes.vcf.gz ./BAJK_1_2016/BAJK_1_2016.freebayes.vcf.gz ./BAJK_2_2016/BAJK_2_2016.freebayes.vcf.gz ./BAJK_4_2016/BAJK_4_2016.freebayes.vcf.gz ./VAST_1_1983/VAST_1_1983.freebayes.vcf.gz ./KRAS_3_2002/KRAS_3_2002.freebayes.vcf.gz > merge.freebayes.vcf.gz
date

(base) [daria@rackham3 freebayes]$ sbatch merging.sh
Submitted batch job 38737351

### Let's get a list of historical files to Merge
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes


#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH -J hist_freebayes
#SBATCH --output=hist_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

date
module load bioinfo-tools bcftools
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes
bcftools merge --threads 4 ./GAST_7_1965/GAST_7_1965.freebayes.vcf.gz ./KALM_2_1955/KALM_2_1955.freebayes.vcf.gz ./GAST_10_1965/GAST_10_1965.freebayes.vcf.gz ./GAST_8_1965/GAST_8_1965.freebayes.vcf.gz ./KALM_5_1961/KALM_5_1961.freebayes.vcf.gz ./ITAL_1_1946/ITAL_1_1946.freebayes.vcf.gz ./GAST_11_1965/GAST_11_1965.freebayes.vcf.gz ./UPPS_3_1969/UPPS_3_1969.freebayes.vcf.gz ./GAST_1_1941/GAST_1_1941.freebayes.vcf.gz ./UPPS_2_1958/UPPS_2_1958.freebayes.vcf.gz ./STOC_4_1965/STOC_4_1965.freebayes.vcf.gz ./SMAL_2_1967/SMAL_2_1967.freebayes.vcf.gz ./KALM_3_1955/KALM_3_1955.freebayes.vcf.gz ./GAST_12_1969/GAST_12_1969.freebayes.vcf.gz ./SMAL_1_1967/SMAL_1_1967.freebayes.vcf.gz ./SMAL_3_1967/SMAL_3_1967.freebayes.vcf.gz ./GAST_2_1941/GAST_2_1941.freebayes.vcf.gz ./STOC_3_1965/STOC_3_1965.freebayes.vcf.gz ./STOC_5_1965/STOC_5_1965.freebayes.vcf.gz ./GAST_3_1941/GAST_3_1941.freebayes.vcf.gz ./DALA_1_1965/DALA_1_1965.freebayes.vcf.gz ./GAST_13_1969/GAST_13_1969.freebayes.vcf.gz ./KALM_7_1961/KALM_7_1961.freebayes.vcf.gz ./UPPS_1_1951/UPPS_1_1951.freebayes.vcf.gz ./KALM_4_1957/KALM_4_1957.freebayes.vcf.gz ./GAST_14_1969/GAST_14_1969.freebayes.vcf.gz ./GAST_5_1943/GAST_5_1943.freebayes.vcf.gz ./KALM_6_1961/KALM_6_1961.freebayes.vcf.gz ./KALM_8_1969/KALM_8_1969.freebayes.vcf.gz ./GAST_4_1941/GAST_4_1941.freebayes.vcf.gz ./GAST_6_1965/GAST_6_1965.freebayes.vcf.gz ./STOC_2_1965/STOC_2_1965.freebayes.vcf.gz ./GAST_9_1965/GAST_9_1965.freebayes.vcf.gz ./STOC_1_1965/STOC_1_1965.freebayes.vcf.gz > merge.freebayes.historical.vcf.gz
date

awk '/# DP/ {flag=1;next} /#/ {flag=0} flag {print}' DALA_1_1965.freebayes.bcftools_stats.txt > depth.txt
awk '/# QUAL/ {flag=1;next} /#/ {flag=0} flag {print}' DALA_1_1965.freebayes.bcftools_stats.txt > qual.txt
