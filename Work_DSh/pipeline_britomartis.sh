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
  WARNING: Skipping mount /var/apptainer/mnt/session/etc/resolv.conf [files]: /etc/resolv.conf doesnt exist in container
  Skipping '??LTA_1_2015-L001_P27562_1051_S51_L001_R1_001.fastq.gz' which didn't exist, or couldn't be read
  Skipping '??LTA_1_2015-L001_P27562_1051_S51_L001_R2_001.fastq.gz' which didn't exist, or couldn't be read

Work dir:
  /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/work/89/dd505f764287b08302fb689a1e24bf

Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

#restart
nextflow run nf-core/sarek -resume magical_gautier --input contemporary_samples.csv -profile uppmax --project naiss2023-5-52 --fasta /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna --skip_tools baserecalibrator --outdir ./results



#Separate diary

Contmporary set:
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes/merge.freebayes.vcf.gz

170G	merge.freebayes.vcf.gz

Historical set:

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.vcf.gz

58G	merge.freebayes.historical.vcf.gz


module load bioinfo-tools vcftools bcftools samtools

for i in $(seq 26 55); do vcftools --gzvcf LR9999${i}allsites_hardTEfiltered_noindel.vcf.gz --remove-indv S15_19M002_S15_19M002 --max-maf 0 --min-meanDP 30 --max-meanDP 80 --recode --stdout | bgzip -c > LR9999${i}allsites_hardTEfiltered_noindel_invariant_Qfilter.vcf.gz ; done


vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.vcf.gz --min-meanDP 10 --max-meanDP 50 --remove-indels --max-missing-count 30 --minQ 30 --recode --stdout | bgzip -c > merge.freebayes.historical.Qfilter.vcf.gz

bcftools stats merge.freebayes.historical.vcf.gz | head -n 30 | tail -n 9

#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -J all_freebayes
#SBATCH --output=all_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

date
module load bioinfo-tools bcftools
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Mapping_Calling_sarek/
bgzip -c >
bcftools merge --threads 4 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes/merge.freebayes.vcf.gz /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.vcf.gz > britomartis.merged.freebayes.vcf.gz
date


bgzip -c >  /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes/merge.freebayes.vcf.gz

bgzip -c > /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.vcf.gz


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
module load bioinfo-tools bcftools samtools
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes
bcftools merge --threads 4 ./SMAL_4_1996/SMAL_4_1996.freebayes.vcf.gz ./KAZA_1_2006/KAZA_1_2006.freebayes.vcf.gz ./POLA_2_2003/POLA_2_2003.freebayes.vcf.gz ./KRAS_4_2002/KRAS_4_2002.freebayes.vcf.gz ./RUSS_3_1999/RUSS_3_1999.freebayes.vcf.gz ./KALM_9_1981/KALM_9_1981.freebayes.vcf.gz ./RUSS_4_2008/RUSS_4_2008.freebayes.vcf.gz ./RUSS_1_1998/RUSS_1_1998.freebayes.vcf.gz ./BELA_1_2004/BELA_1_2004.freebayes.vcf.gz ./RUSS_5_2008/RUSS_5_2008.freebayes.vcf.gz ./KAZA_2_2006/KAZA_2_2006.freebayes.vcf.gz ./URAL_2_1995/URAL_2_1995.freebayes.vcf.gz ./URAL_1_1991/URAL_1_1991.freebayes.vcf.gz ./KALM_1_2018/KALM_1_2018.freebayes.vcf.gz ./KALM_10_1983/KALM_10_1983.freebayes.vcf.gz ./VAST_3_2005/VAST_3_2005.freebayes.vcf.gz ./CHEH_1_1989/CHEH_1_1989.freebayes.vcf.gz ./BAJK_3_2016/BAJK_3_2016.freebayes.vcf.gz ./JAPA_2_1994/JAPA_2_1994.freebayes.vcf.gz ./VORO_1_1998/VORO_1_1998.freebayes.vcf.gz ./SMAL_5_1998/SMAL_5_1998.freebayes.vcf.gz ./RUSS_2_1999/RUSS_2_1999.freebayes.vcf.gz ./KRAS_2_2002/KRAS_2_2002.freebayes.vcf.gz ./ALTA_2_2015/ALTA_2_2015.freebayes.vcf.gz ./SLOV_1_1990/SLOV_1_1990.freebayes.vcf.gz ./POLA_1_2003/POLA_1_2003.freebayes.vcf.gz ./JAPA_1_1994/JAPA_1_1994.freebayes.vcf.gz ./ALTA_1_2015/ALTA_1_2015.freebayes.vcf.gz ./CHEH_3_2008/CHEH_3_2008.freebayes.vcf.gz ./KRAS_1_2002/KRAS_1_2002.freebayes.vcf.gz ./SMAL_6_2013/SMAL_6_2013.freebayes.vcf.gz ./STOC_6_1965/STOC_6_1965.freebayes.vcf.gz ./CHEH_2_2004/CHEH_2_2004.freebayes.vcf.gz ./VAST_2_1999/VAST_2_1999.freebayes.vcf.gz ./BAJK_1_2016/BAJK_1_2016.freebayes.vcf.gz ./BAJK_2_2016/BAJK_2_2016.freebayes.vcf.gz ./BAJK_4_2016/BAJK_4_2016.freebayes.vcf.gz ./VAST_1_1983/VAST_1_1983.freebayes.vcf.gz ./KRAS_3_2002/KRAS_3_2002.freebayes.vcf.gz > merge.contemporary.freebayes.vcf.gz
bgzip -c merge.contemporary.freebayes.vcf.gz > merge.contemporary.freebayes.bg.vcf.gz
date


plink --vcf ../merge.freebayes.historical.Qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out britomartis.historical


plink --vcf ../merge.freebayes.historical.Qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.historical.prune.in --make-bed --pca --out britomartis.historical


tabix /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes/merge.contemporary.freebayes.bg.vcf.gz
tabix /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.bg.vcf.gz > britomartis.merged.freebayes.bg.vcf.gz


#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -J all_freebayes
#SBATCH --output=all_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

date
module load bioinfo-tools bcftools samtools
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Mapping_Calling_sarek/06_freebayesRun/
tabix /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes/merge.contemporary.freebayes.bg.vcf.gz
tabix /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.bg.vcf.gz
bcftools merge --threads 8 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes/merge.contemporary.freebayes.bg.vcf.gz /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/05_CallingHistorical/results/variant_calling/freebayes/merge.freebayes.historical.bg.vcf.gz > britomartis.merged.freebayes.vcf.gz
date


bcftools stats merge.contemporary.freebayes.bg.vcf.gz | head -n 30 | tail -n 9

vcftools --gzvcf merge.freebayes.historical.Qfilter.vcf.gz --max-missing-count 30 --minQ 30 --minDP 8 --maxDP 40 --recode --stdout | gzip -c | bcftools stats | head -n 30 | tail -n 9 > stats.out

vcftools --gzvcf britomartis.merged.freebayes.vcf.gz --remove-indels --max-missing-count 10 --minQ 30 --min-meanDP 10 --max-meanDP 40 --recode --stdout | gzip -c > britomartis.merged.freebayes.qfilter.vcf.gz

--min-meanDP 30

3607 britomartis.merged.freebayes.qfilter.vcf.gz


plink --vcf ../britomartis.merged.freebayes.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out britomartis.all
Pruning complete.  1486 of 2203 variants removed.


plink --vcf ../britomartis.merged.freebayes.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out britomartis.all



########## Contmporary only!
module load bioinfo-tools bcftools samtools
bcftools stats merge.contemporary.freebayes.bg.vcf.gz | head -n 30 | tail -n 9

interactive -A naiss2023-5-52  -n 4 -t 10:00:0

vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz --remove-indels --max-missing-count 3 --minQ 30 --min-meanDP 8 --max-meanDP 50 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.qfilter.vcf.gz

bcftools stats merge.contemporary.freebayes.bg.qfilter.vcf.gz | head -n 30 | tail -n 9

SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      81478
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 81478
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   11326


plink --vcf merge.contemporary.freebayes.bg.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out britomartis.all

#### Tailoring filters with subsampled Set
interactive -A naiss2023-5-52  -n 4 -t 10:00:0
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/results/variant_calling/freebayes

bcftools view merge.contemporary.freebayes.bg.vcf.gz | vcfrandomsample -r 0.012 > merge.contemporary.freebayes.bg.subset.vcf


#subsamp.bed
HG992177.1 1 10000000

bcftools filter -T ^subsamp.bed -o merge.contemporary.freebayes.bg.Chr1.vcf.gz merge.contemporary.freebayes.bg.vcf.gz

bcftools stats merge.contemporary.freebayes.bg.Chr1.vcf.gz | head -n 30 | tail -n 9
[E::bgzf_read_block] Failed to read BGZF block data at offset 201326170 expected 7991 bytes; hread returned 404
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      1153951
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 770306
SN      0       number of MNPs: 215555
SN      0       number of indels:       182319
SN      0       number of others:       139468
SN      0       number of multiallelic sites:   258348

#Example:
vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz --remove-indels --max-missing-count 3 --minQ 30 --min-meanDP 8 --max-meanDP 50 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.qfilter.vcf.gz

#Only missingness
vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing-count 3 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.vcf.gz.mxmiss3.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.vcf.gz.mxmiss3.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      630
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 630
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   84

vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing-count 10 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.vcf.gz.mxmiss10.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.vcf.gz.mxmiss10.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      4071
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 4071
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   611

vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing-count 30 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.vcf.gz.mxmiss30.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.vcf.gz.mxmiss30.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      19175
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 19175
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   3137

vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing 1 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.mxmissprc1.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.mxmissprc1.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      196
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 196
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   24

vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing 0.75 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.mxmissprc075.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.mxmissprc075.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      8727
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 8727
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   1360

vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing 0.5 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.mxmissprc05.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.mxmissprc05.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      39
SN      0       number of records:      27557
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 27557
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   4686


### Trying to recall
### Example
cd ~
bcftools mpileup -a AD,DP,SP -Ou -f $REF \
./align/*_sort.bam | bcftools call -f GQ,GP \
-mO z -o ./cichlid.vcf.gz

#### My code/ small recalling
module load bioinfo-tools bcftools samtools vcftools

bcftools mpileup -a AD,DP,SP -Ou -f $REF -b reduced_samle.list | bcftools call -f GQ,GP -mO z -o ./cichlid.vcf.gz


nano reduced_samle.list
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/URAL_2_1995/URAL_2_1995.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/RUSS_4_2008/RUSS_4_2008.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/JAPA_1_1994/JAPA_1_1994.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/CHEH_1_1989/CHEH_1_1989.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/KRAS_2_2002/KRAS_2_2002.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/RUSS_1_1998/RUSS_1_1998.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/URAL_1_1991/URAL_1_1991.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/STOC_6_1965/STOC_6_1965.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/VAST_1_1983/VAST_1_1983.md.cram
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/01_MappingAll/results/preprocessing/markduplicates/SMAL_4_1996/SMAL_4_1996.md.cram

 mkdir 00_mpileupCalling

 bcftools mpileup -Ou -Q 30 -q 30 -B -f /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna -b reduced_samle.list | bcftools call -m -M -O b --threads 4 -o ./reduced_samle_contemp.vcf.gz

 .....testing interactively.....

# Preparing comparable set
vcftools --gzvcf merge.contemporary.freebayes.bg.Chr1.vcf.gz --remove-indels --max-missing 0.5 --recode --stdout  | gzip -c > merge.contemporary.freebayes.bg.Chr1.mxmissprc05.vcf.gz
bcftools stats merge.contemporary.freebayes.bg.Chr1.mxmissprc05.vcf.gz | head -n 30 | tail -n 9


 might be wrong: vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz  --indv URAL_2_1995_URAL_2_1995,RUSS_4_2008_RUSS_4_2008,JAPA_1_1994_JAPA_1_1994,CHEH_1_1989_CHEH_1_1989,KRAS_1_2002_KRAS_1_2002,RUSS_1_1998_RUSS_1_1998,URAL_1_1991_URAL_1_1991,STOC_6_1965_STOC_6_1965,VAST_1_1983_VAST_1_1983,SMAL_4_1996_SMAL_4_1996 --recode --stdout | gzip -c > merge.contemporary.reduced.freebayes.bg.vcf.gz

 vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz --indv URAL_2_1995_URAL_2_1995;RUSS_4_2008_RUSS_4_2008 --recode --stdout

 --keep <filename>

Provide files containing a list of individuals to either include or exclude in subsequent analysis. Each individual ID (as defined in the VCF headerline) should be included on a separate line. If both options are used, then the "--keep" option is executed before the "--remove" option. When multiple files are provided, the union of individuals from all keep files subtracted by the union of individuals from all remove files are kept. No header line is expected.

nano keep.reduced.txt
#
URAL_2_1995_URAL_2_1995
RUSS_4_2008_RUSS_4_2008
JAPA_1_1994_JAPA_1_1994
CHEH_1_1989_CHEH_1_1989
KRAS_1_2002_KRAS_1_2002
RUSS_1_1998_RUSS_1_1998
URAL_1_1991_URAL_1_1991
STOC_6_1965_STOC_6_1965
VAST_1_1983_VAST_1_1983
SMAL_4_1996_SMAL_4_1996

vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz --keep keep.reduced.txt --recode --stdout | gzip -c > merge.contemporary.reduced.freebayes.bg.vcf.gz


###*** --max-indv
bcftools stats merge.contemporary.reduced.freebayes.bg.vcf.gz | head -n 30 | tail -n 9




bcftools stats reduced_samle_contemp.vcf.gz | head -n 30 | tail -n 9

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      10
SN      0       number of records:      96145268
SN      0       number of no-ALTs:      86364234
SN      0       number of SNPs: 9012258
SN      0       number of MNPs: 0
SN      0       number of indels:       768776
SN      0       number of others:       0
SN      0       number of multiallelic sites:   106241

#Only missingness
vcftools --gzvcf reduced_samle_contemp.vcf.gz --remove-indels --max-missing-count 3 --recode --stdout  | gzip -c > reduced_samle_contemp.mxmiss3.vcf.gz
bcftools stats reduced_samle_contemp.mxmiss3.vcf.gz | head -n 30 | tail -n 9


# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      10
SN      0       number of records:      95220215
SN      0       number of no-ALTs:      86207957
SN      0       number of SNPs: 9012258
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   92674

vcftools --gzvcf reduced_samle_contemp.vcf.gz --remove-indels --max-missing-count 3 --recode --stdout  | gzip -c > reduced_samle_contemp.mxmiss3.vcf.gz

#subsamp.bed
HG992177.1 1 10000000

bcftools filter -T subsamp.bed -o reduced_samle_contemp.Chr1.vcf.gz reduced_samle_contemp.vcf.gz

bcftools stats reduced_samle_contemp.Chr1.vcf.gz | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	10
SN	0	number of records:	8447766
SN	0	number of no-ALTs:	7581371
SN	0	number of SNPs:	794417
SN	0	number of MNPs:	0
SN	0	number of indels:	71978
SN	0	number of others:	0
SN	0	number of multiallelic sites:	9464

vcftools --gzvcf reduced_samle_contemp.Chr1.vcf.gz --remove-indels --max-missing-count 3 --recode --stdout  | gzip -c > reduced_samle_contemp.Chr1.mxmiss3.vcf.gz
bcftools stats reduced_samle_contemp.Chr1.mxmiss3.vcf.gz | head -n 30 | tail -n 9

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	10
SN	0	number of records:	8361676
SN	0	number of no-ALTs:	7567259
SN	0	number of SNPs:	794417
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	8187

vcftools --gzvcf reduced_samle_contemp.Chr1.vcf.gz --remove-indels --max-missing 1 --recode --stdout  | gzip -c > reduced_samle_contemp.Chr1.mxmissprc1.vcf.gz
bcftools stats reduced_samle_contemp.Chr1.mxmissprc1.vcf.gz | head -n 30 | tail -n 9


bcftools stats reduced_samle_contemp.Chr1.mxmissprc1.vcf.gz | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	10
SN	0	number of records:	8361676
SN	0	number of no-ALTs:	7567259
SN	0	number of SNPs:	794417
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	8187

vcftools --gzvcf reduced_samle_contemp.Chr1.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --min-meanDP 8 --max-meanDP 50 --recode --stdout  | gzip -c > reduced_samle_contemp.Chr1.qfilter.vcf.gz

0 snps

vcftools --gzvcf reduced_samle_contemp.Chr1.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > reduced_samle_contemp.Chr1.qfilter.vcf.gz

bcftools stats reduced_samle_contemp.Chr1.qfilter.vcf.gz | head -n 30 | tail -n 9

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	10
SN	0	number of records:	665439
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	665439
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0

### Conclusion: consensus caller dont give missing data as ./.
### Allelic inbalance???

### Next: use mpiledup file for Chr1 (first half) to run stats
module load plink
plink --vcf reduced_samle_contemp.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out britomartis.contemp.reduced

plink --vcf reduced_samle_contemp.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.contemp.reduced.prune.in --make-bed --pca --out britomartis.contemp.reduced.pruned

plink --vcf reduced_samle_contemp.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out britomartis.contemp.reduced.noprune

##### Small test

plink --vcf reduced_samle_contemp.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 50 0.5 --out britomartis.contemp.reduced.test
Pruned 465453 variants from chromosome 27, leaving 199986.
Pruning complete.  465453 of 665439 variants removed.

### Call using different algorithm

#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -J all_freebayes
#SBATCH --output=all_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

date
module load bioinfo-tools bcftools samtools
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/00_mpileupCalling
bcftools mpileup -Ou -Q 30 -q 30 -B -f /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna -b reduced_samle.list | bcftools call -m -M -O b --threads 8 -o ./reduced_samle_contemp_mcalling.vcf.gz
date

wc -l reduced_samle_contemp_mcalling.vcf.gz
21734211 reduced_samle_contemp_mcalling.vcf.gz

bcftools filter -T subsamp.bed -o reduced_samle_contemp_mcalling.Chr1.vcf.gz reduced_samle_contemp_mcalling.vcf.gz

bcftools stats reduced_samle_contemp_mcalling.Chr1.vcf.gz | head -n 30 | tail -n 9

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	10
SN	0	number of records:	8447766
SN	0	number of no-ALTs:	7569410
SN	0	number of SNPs:	802171
SN	0	number of MNPs:	0
SN	0	number of indels:	76185
SN	0	number of others:	0
SN	0	number of multiallelic sites:	45255

vcftools --gzvcf reduced_samle_contemp_mcalling.Chr1.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz

bcftools stats reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	10
SN	0	number of records:	521907
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	521907
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0

### Next: use mpiledup file for Chr1 (first half) to run stats
module load plink
plink --vcf reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 40 0.2 --out britomartis.contemp_mcalling.reduced

plink --vcf reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.contemp.reduced.prune.in --make-bed --pca --out britomartis.contemp_mcalling.reduced.pruned

plink --vcf reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out britomartis.contemp_mcalling.reduced.noprune

#### Aug 25 ####
### Trying new subsample

interactive -A naiss2023-5-52 -n 4 -t 10:00:0

##contig=<ID=HG992196.1,length=18366132>

#subsampChr20.bed
HG992196.1  1 18366132

module load bioinfo-tools bcftools samtools vcftools

nohup bcftools filter -T subsampChr20.bed -o reduced_samle_contemp_mcalling.Chr20.vcf.gz reduced_samle_contemp_mcalling.vcf.gz &


###Chrom 10
##contig=<ID=HG992186.1,length=21375889>

#subsampChr10.bed
HG992186.1  1 21375889

nohup bcftools filter -T subsampChr10.bed -t 3 -o reduced_samle_contemp_mcalling.Chr10.vcf.gz reduced_samle_contemp_mcalling.vcf.gz &

# Filtering vcf

vcftools --gzvcf reduced_samle_contemp_mcalling.Chr1.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz


### Aug 28

#Trying to see all the SNPs


###Current corrupted??? file ends with:
HG992194.1	4981543	.	G	.	128.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981544	.	C	.	131.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981545	.	G	.	131.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981546	.	A	.	131.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981547	.	A	.	123.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981548	.	C	.	98.0878	.	DP=6;MQ0F=0;AN=4;DP4=1,1,0,0;MQ=37	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981549	.	T	.	123.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981550	.	T	.	131.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981551	.	A	.	131.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/0
HG992194.1	4981552	.	T	.	131.088	.	DP=6;MQ0F=0;AN=4;DP4=1,2,0,0;MQ=38	GT	./.	0/0	./.	./.	./.	./.	./.	./.	./.	0/00_Mapping_Calling_sarek


#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH -J Filteringcontemp_mcalling.out
#SBATCH --output=cont_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bcftools samtools vcftools
vcftools --gzvcf reduced_samle_contemp_mcalling.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > reduced_samle_contemp_mcalling.qfilter.vcf.gz

module load plink
plink --vcf reduced_samle_contemp_mcalling.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 40 0.2 --out britomartis.contemp_mcalling.full

Pruning complete.  17656071 of 19539605 variants removed.

plink --vcf reduced_samle_contemp_mcalling.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.contemp_mcalling.full.prune.in --make-bed --pca --out britomartis.contemp_mcalling.full.pruned

1883534 variants and 10 people pass filters and QC.

plink --vcf reduced_samle_contemp_mcalling.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out britomartis.contemp_mcalling.full.noprune

19539605 variants and 10 people pass filters and QC.




#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH -J Filteringcontemp_mcalling.out
#SBATCH --output=cont_freebayes.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bcftools samtools vcftools
vcftools --gzvcf reduced_samle_contemp_mcalling.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --min-meanDP 8 --max-meanDP 50 --recode --stdout  | gzip -c > reduced_samle_contemp_mcalling.qfilterDP.vcf.gz



#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J all_sample_mcalling
#SBATCH --output=all_sample_mcalling.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bcftools samtools vcftools
bcftools mpileup -Ou -Q 30 -q 30 -B -f /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna -b full_sample.list | bcftools call -m -M -O b --threads 8 -o britomartis_all_mcaller.vcf.gz

vcftools --gzvcf britomartis_all_mcaller.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > britomartis_all_mcaller.miss3.vcf.gz

vcftools --gzvcf britomartis_all_mcaller.vcf.gz --remove-indels --max-missing-count 20 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > britomartis_all_mcaller.miss20.vcf.gz



### full_sample.list
### Getting sample list from SAREK



awk -F ',' {print $3} full_sample_raw.list

(base) [daria@rackham2 07_mpileupM_joined]$ nano all_sample_mcalling.sh
(base) [daria@rackham2 07_mpileupM_joined]$ sbatch all_sample_mcalling.sh
Subted batch job 40558801


du -h reduced_samle_contemp_mcalling.qfilterDP.vcf.gz
#Nothing here!!! Filter not applicable!!!!

bcftools stats reduced_samle_contemp_mcalling.qfilterDP.vcf.gz | head -n 30 | tail -n 9


 bcftools filter -e 'FORMAT/PL[0] < 50'  reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz


scancel 40558801


bcftools mpileup -a AD,DP,SP -Ou -f $REF -b reduced_samle.list | bcftools call -f GQ,GP -mO z -o ./cichlid.vcf.gz

#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J all_sample_mcalling
#SBATCH --output=all_sample_mcalling.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bcftools samtools vcftools
bcftools mpileup -a AD,DP,SP -Ou -Q 30 -q 30 -B -f /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna -b full_sample.list | bcftools call -f GQ,GP -m -M -O b --threads 8 -o britomartis_all_mcaller.vcf.gz

vcftools --gzvcf britomartis_all_mcaller.vcf.gz --remove-indels --max-missing-count 3 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > britomartis_all_mcaller.miss3.vcf.gz

vcftools --gzvcf britomartis_all_mcaller.vcf.gz --remove-indels --max-missing-count 20 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > britomartis_all_mcaller.miss20.vcf.gz

sbatch all_sample_mcalling.sh
Submitted batch job 40561218

### August 31

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/00_mpileupCalling

zcat britomartis_all_mcaller.vcf.gz | tail -10

#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -J end_check
#SBATCH --output=end_check.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

zcat britomartis_all_mcaller.vcf.gz | tail -10 > current.tail

HG992197.1      2323651 .       A       G       193.215 .       DP=329;VDB=0.517259;SGB=29.0063;RPBZ=-1.34465;MQBZ=-0.868137;MQSBZ=0.320284;BQBZ=-1.3369;SCBZ=-1.19565;MQ0F=0;AC=2;AN=96;DP4=78,93,3,7;MQ=56    GT:PL:DP:SP:AD:GP:GQ    0/0:0,3,46:1:0:1,0:0.979118,0.0208818,1.11337e-08:16    0/0:0,3,60:1:0:1,0:0.979118,0.0208818,4.4324e-10:16     0/0:0,15,205:5:0:5,0:0.998656,0.00134384,1.42962e-24:28 0/0:0,18,255:6:0:6,0:0.999326,0.000673968,1.43058e-29:31        ./.:0,0,0:0:0:0,0:0,0,0:0       0/0:0,18,211:6:0:6,0:0.999326,0.000673968,3.59345e-25:31        0/0:0,39,255:13:0:13,0:0.999995,5.3571e-06,1.43154e-29:52       0/0:0,9,127:3:0:3,0:0.994671,0.00532858,8.98429e-17:22  0/0:0,9,167:3:0:3,0:0.994671,0.00532858,8.98429e-21:22  0/0:0,27,255:9:0:9,0:0.999915,8.48976e-05,1.43142e-29:40        0/0:0,3,46:1:0:1,0:0.979118,0.0208818,1.11337e-08:16    0/0:0,27,255:9:0:9,0:0.999915,8.48976e-05,1.43142e-29:40        ./.:0,0,0:0:0:0,0:0,0,0:0       0/
1-21:42:29

Benchmarking:
25+25+24+23+23+22+21+21+21+21+20+20+20+20+20+20+20
 404 million snps in 2 days


bcftools stats britomartis_all_mcaller.vcf.gz | head -n 30 | tail -n 9

#September 1st

#Example:


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:00:00
#SBATCH -J QCall
#SBATCH --output=QCall.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools vcftools
for i in $(seq 24 55); do
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel.vcf.gz --mac 1 --recode --stdout | gzip -c > LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --freq2 --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant --max-alleles 2
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --depth --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --site-mean-depth --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --site-quality --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --missing-indv --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --missing-site --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --het --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
done


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J SNP_QC
#SBATCH --output=SNP_QC.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools vcftools
vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel.vcf.gz --mac 1 --recode --stdout | gzip -c > LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --freq2 --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant --max-alleles 2
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --depth --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --site-mean-depth --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --site-quality --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --missing-indv --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --missing-site --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
  vcftools --gzvcf /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/LR9999${i}allsites_hardTEfiltered_noindel_variant.vcf.gz --het --out /proj/uppstore2017185/b2014034/nobackup/Dasha/ShortLong_Vanessa/6_SCANS/1_VCFtoolsQC/LR9999${i}_variant
done


#### Testing stringent
module load plink
plink --vcf ../britomartis_all_mcaller.miss3.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 40 0.2 --out britomartis_all_mcaller.miss3

Pruning complete.  112134 of 220333 variants removed.

plink --vcf ../britomartis_all_mcaller.miss3.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis_all_mcaller.miss3.prune.in  --make-bed --pca --out britomartis_all_mcaller.miss3.pruned

108199 variants and 73 people pass filters and QC.

plink --vcf ../britomartis_all_mcaller.miss3.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out britomartis_all_mcaller.miss3.noprune

220333 variants and 73 people pass filters and QC.

bcftools stats britomartis_all_mcaller.miss3.vcf.gz | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	73
SN	0	number of records:	220333
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	220333
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0

##Continue filtering
wc -l britomartis_all_mcaller.miss20.vcf.gz
20026630 britomartis_all_mcaller.miss20.vcf.gz

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined

#### Testing relaxed
module load bioinfo-tools  plink
plink --vcf ../britomartis_all_mcaller.miss20.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 40 0.2 --out britomartis_all_mcaller.miss20

Pruning complete.  4967863 of 8419518 variants removed.


plink --vcf ../britomartis_all_mcaller.miss20.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis_all_mcaller.miss20.prune.in  --make-bed --pca --out britomartis_all_mcaller.miss20.pruned

3451655 variants and 73 people pass filters and QC.

plink --vcf ../britomartis_all_mcaller.miss20.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out britomartis_all_mcaller.miss20.noprune

8419518 variants and 73 people pass filters and QC.

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J SNP_QC
#SBATCH --output=SNP_QC.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools vcftools


vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --freq2 --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --depth --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --site-mean-depth --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --site-quality  --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --missing-indv  --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --missing-site  --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.vcf.gz --het --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller


### Split dataset!
--remove-indv S15_19M002_S15_19M002

module load bcftools

bcftools filter -T ^VcardDToL_kNdSannot_CDS.bed -o LR9999${i}allsites_hardTEfiltered_noindel_variant_Qfilter_nosinglt_CDS.vcf.gz LR9999${i}allsites_hardTEfiltered_noindel_variant_Qfilter_nosinglt.vcf.gz

HG992207.1	W	5.27	37.4
HG992176.1	Z	26.23	34
HG992208.1	MT	0.02	19.7

##contig=<ID=HG992207.1,length=5265952>
##contig=<ID=HG992176.1,length=26233870>
##contig=<ID=HG992208.2,length=15154>

### mtDNA.bed
HG992208.1  1 15154

bcftools filter -T mtDNA.bed -o britomartis_all_mcaller.miss3.mtDNA.vcf.gz britomartis_all_mcaller.miss3.vcf.gz

### sex_chroms.bed
HG992207.1  1 5265952
HG992176.1  1 26233870

-r, --regions
bcftools view -r HG992208.2 britomartis_all_mcaller.miss3.vcf.gz

bcftools view -r HG992207.1,HG992176.1 britomartis_all_mcaller.miss3.vcf.gz -o britomartis_all_mcaller.miss3.ZW.vcf.gz

bcftools view -r HG992177.1,HG992178.1,HG992179.1,HG992180.1,HG992181.1,HG992182.1,HG992183.1,HG992184.1,HG992185.1,HG992186.1,HG992187.1,HG992188.1,HG992189.1,HG992190.1,HG992191.1,HG992192.1,HG992193.1,HG992194.1,HG992195.1,HG992196.1,HG992197.1,HG992198.1,HG992199.1,HG992200.1,HG992201.1,HG992202.1,HG992203.1,HG992204.1,HG992205.1,HG992206.1 britomartis_all_mcaller.miss3.vcf.gz -o britomartis_all_mcaller.miss3.autosomes.vcf.gz


wc -l britomartis_all_mcaller.miss3.autosomes.vcf.gz
585531 britomartis_all_mcaller.miss3.autosomes.vcf.gz
wc -l britomartis_all_mcaller.miss3.ZW.vcf.gz
21299 britomartis_all_mcaller.miss3.ZW.vcf.gz
wc -l britomartis_all_mcaller.miss3.mtDNA.vcf.gz
3617 britomartis_all_mcaller.miss3.mtDNA.vcf.gz


#### mtDNA
module load bioinfo-tools  plink

plink --vcf ../britomartis_all_mcaller.miss3.mtDNA.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out britomartis_all_mcaller.miss3.mtDNA

8419518 variants and 73 people pass filters and QC.


#mtDNA tree

mv britomartis_all_mcaller.miss20.mtDNA.COI.vcf.gz britomartis_all_mcaller.COI.vcf.gz

### CReating alternative rferences

#Dictionary for fasta
gatk CreateSequenceDictionary -R GCA_905220545.2_ilMelAtha1.2_genomic.fna

#Making rfs
gatk FastaAlternateReferenceMaker \
  -R /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna \
  -O britomartisCO1.fasta \
  -V /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.COI.vcf.gz

bcftools index ../britomartis_all_mcaller.COI.vcf.gz
tabix -p vcf ../britomartis_all_mcaller.COI.vcf.gz


### This gave me full gnome and for all the samples together

###
cp /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna
cp /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.dict .

##Let's extract piece of fasta:


nano CO1.bed
# HG992208.2  4154  4810
module load  BEDTools/2.29.2
bedtools getfasta -fi GCA_905220545.2_ilMelAtha1.2_genomic.fna -bed CO1.bed


>HG992208.2:4154-4810
AATAAATGTTGATATAAAATAGGATCACCTCCCCCGGCTGGGTCAAAGAATGAGGTATTAATATTTCGATCGGTAAGAAGTATAGTAATTGCTCCAGCTAATACTGGCAAAGATAATAATAATAAGAGAGCTGTAATACCCACAGCTCAAACAAATAAAGGTATTTGATCAAATGATATATTATTAACACGTATATTAATAATTGTAGTAATAAAATTAATAGCTCCTAAGATTGATGAAATTCCAGCTAAATGTAATGAAAAAATTGCTAAATCAACAGATGATCCTCTATGAGCAATATTAGATGAAAGTGGGGGGTAAACTGTTCATCCTGTTCCTGCTCCATTTTCTACAATTCTTCTGGAAATTAATAAAATTAGTGAGGGGGGGAGTAATCAAAATCTTATATTATTTATTCGAGGGAATGCTATATCAGGAGCTCCTAATATTAAAGGAACTAATCAATTTCCAAATCCTCCAATTATAATAGGTATAACCATAAAAAAAATTATAATAAAAGCATGAGCTGTTACAATAGTATTATAAATTTGATCATCTCCAATTAAAGATCCAGGATTTCCTAATTCAGTTCGAATTAAAAGTCTAAGTGAAGTTCCTAATATACCTGCTCAAATTCCAAAAATAAAATATAAAGT

gatk CreateSequenceDictionary -R GCA_905220545.2_ilMelAtha1.2_CO1.fna
module load samtools
samtools faidx GCA_905220545.2_ilMelAtha1.2_CO1.fna

gatk FastaAlternateReferenceMaker \
  -R GCA_905220545.2_ilMelAtha1.2_CO1.fna \
  -O britomartisCO1_small.fasta \
  -V /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.COI.vcf.gz


NOPE: bcftools view -r HG992208.2:4154-4810 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.COI.vcf.gz


bcftools view -s sample_name input.vcf -o output_sample_name.vcf

for sample in $(bcftools query -l britomartis_all_mcaller.COI.vcf.gz); do
    bcftools view -s $sample britomartis_all_mcaller.COI.vcf.gz -o ${sample}.vcf
done


for sample in $(bcftools query -l britomartis_all_mcaller.COI.vcf.gz); do
      bgzip ${sample}.vcf;
      tabix -p vcf ${sample}.vcf.gz
done


for sample in $(bcftools query -l britomartis_all_mcaller.COI.vcf.gz); do
  gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O britomartisCO1${sample}.fasta -V ${sample}.vcf.gz
done





#!/bin/bash

#SBATCH --job-name=split_vcf
#SBATCH --output=split_vcf_%A_%a.out
#SBATCH --error=split_vcf_%A_%a.err
#SBATCH --array=1-N  # Replace N with the number of samples

# Load necessary modules or activate environments, if needed
# module load bcftools (if using modules)

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.list)

bcftools view -s $SAMPLE input.vcf -o ${SAMPLE}.vcf

# If you want to compress and index the resulting VCFs
# bgzip ${SAMPLE}.vcf
# tabix -p vcf ${SAMPLE}.vcf.gz



gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O britomartisCO1_ALTA_1_2015_ALTA_1_2015.fasta -V ALTA_1_2015_ALTA_1_2015.vcf.gz

sed 's/^>46 \(HG992208.2\):[0-9]*-[0-9]*/>\1/' britomartisCO1_ALTA_1_2015_ALTA_1_2015.fasta > britomartisCO1_ALTA_1_2015_ALTA_1_2015_headfix.fasta

bedtools getfasta -fi britomartisCO1_ALTA_1_2015_ALTA_1_2015_headfix.fasta -bed CO1.bed -fo britomartisCO1_ALTA_1_2015_ALTA_1_2015_headfix.f.fasta


#Try:
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined

module load bioinfo-tools bcftools GATK

gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O britomartisCO1_CHEH_1_1989_CHEH_1_1989.fasta -V CHEH_1_1989_CHEH_1_1989.vcf.gz

sed 's/^>46 \(HG992208.2\):[0-9]*-[0-9]*/>\1/' britomartisCO1_CHEH_1_1989_CHEH_1_1989.fasta > britomartisCO1_CHEH_1_1989_CHEH_1_1989_headfix.fasta

bedtools getfasta -fi britomartisCO1_CHEH_1_1989_CHEH_1_1989_headfix.fasta -bed CO1.bed -fo britomartisCO1_CHEH_1_1989_CHEH_1_1989_headfix.f.fasta


#Making array jobs

bcftools query -l britomartis_all_mcaller.COI.vcf.gz -o samples.txt

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J extractCO1
#SBATCH --output=extractCO1.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-$(cat samples.txt | wc -l)

# Load modules or source software here, if needed
module load bioinfo-tools bcftools GATK BEDTools

# Extract sample name based on the current array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Run the GATK command
gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O britomartisCO1_${SAMPLE}.fasta -V ${SAMPLE}.vcf.gz

# Use sed to adjust the header
sed 's/^>46 \(HG992208.2\):[0-9]*-[0-9]*/>\1/' britomartisCO1_${SAMPLE}.fasta > britomartisCO1_${SAMPLE}_headfix.fasta

# Run bedtools to extract the sequence
bedtools getfasta -fi britomartisCO1_${SAMPLE}_headfix.fasta -bed CO1.bed -fo britomartisCO1_${SAMPLE}_headfix_final.fasta



#!/bin/bash

# Define the output file
output_file="combined.fasta"

# Ensure the output file is empty to start
> "$output_file"

# Loop over all fasta files with the naming pattern
for fasta in britomartisCO1_*_headfix_final.fasta; do
    # Extract ID from the filename
    ID="${fasta#britomartisCO1_}"
    ID="${ID%_headfix_final.fasta}"

    # Take only the first half of the redundant ID
    ID=$(echo "$ID" | awk -F'_' '{print $1"_"$2"_"$3}')

    # Use sed to modify the header and append to the output file
    sed "s/^>HG992208.2:[0-9]*-[0-9]*$/>$ID/" "$fasta" >> "$output_file"
done


for sample in $(bcftools query -l britomartis_all_mcaller.COI.vcf.gz); do
  gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O britomartisCO1${sample}.fasta -V ${sample}.vcf.gz
done


bcftools view -i 'TYPE="snp" || TYPE="indel"' ALTA_1_2015_ALTA_1_2015.vcf.gz -Oz -o outputALTA_1_2015_ALTA_1_2015.vcf.gz
bcftools view -e 'GT="0/0"' outputALTA_1_2015_ALTA_1_2015.vcf.gz -Oz -o outputALTA_1_2015_ALTA_1_2015_alt.vcf.gz

tabix -p vcf outputALTA_1_2015_ALTA_1_2015_alt.vcf.gz

# Run the GATK command
gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O outputALTA_1_2015_ALTA_1_2015_alt.fasta -V outputALTA_1_2015_ALTA_1_2015_alt.vcf.gz

# Use sed to adjust the header
sed 's/^>46 \(HG992208.2\):[0-9]*-[0-9]*/>\1/' outputALTA_1_2015_ALTA_1_2015_alt.fasta > outputALTA_1_2015_ALTA_1_2015_alt_headfix.fasta

# Run bedtools to extract the sequence
bedtools getfasta -fi outputALTA_1_2015_ALTA_1_2015_alt_headfix.fasta -bed CO1.bed -fo outputALTA_1_2015_ALTA_1_2015_alt_headfix_final.fasta





bcftools view -i 'TYPE="snp" || TYPE="indel"' BAJK_2_2016_BAJK_2_2016.vcf.gz -Oz -o outputBAJK_2_2016_BAJK_2_2016.vcf.gz
bcftools view -e 'GT="0/0"' outputBAJK_2_2016_BAJK_2_2016.vcf.gz -Oz -o outputBAJK_2_2016_BAJK_2_2016_alt.vcf.gz

tabix -p vcf outputBAJK_2_2016_BAJK_2_2016_alt.vcf.gz

# Run the GATK command
gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O outputBAJK_2_2016_BAJK_2_2016_alt.fasta -V outputBAJK_2_2016_BAJK_2_2016_alt.vcf.gz

# Use sed to adjust the header
sed 's/^>46 \(HG992208.2\):[0-9]*-[0-9]*/>\1/' outputBAJK_2_2016_BAJK_2_2016_alt.fasta > outputBAJK_2_2016_BAJK_2_2016_alt_headfix.fasta

# Run bedtools to extract the sequence
bedtools getfasta -fi outputBAJK_2_2016_BAJK_2_2016_alt_headfix.fasta -bed CO1.bed -fo outputBAJK_2_2016_BAJK_2_2016_alt_headfix_final.fasta


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J extractCO1
#SBATCH --output=extractCO1.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-73

# Load modules or source software here, if needed
module load bioinfo-tools bcftools GATK BEDTools

# Extract sample name based on the current array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Filter vcf file for invariant sites, enforce removing 0/0 records
bcftools view -i 'TYPE="snp" || TYPE="indel"' ${SAMPLE}.vcf.gz -Oz -o ${SAMPLE}.snponly.vcf.gz
bcftools view -e 'GT="0/0"' ${SAMPLE}.snponly.vcf.gz -Oz -o ${SAMPLE}.no0.vcf.gz

#Index resulting vcf files
tabix -p vcf ${SAMPLE}.no0.vcf.gz

# Run the GATK command
gatk FastaAlternateReferenceMaker -R GCA_905220545.2_ilMelAtha1.2_genomic.fna -O britomartisCO1_${SAMPLE}_f.fasta -V ${SAMPLE}.no0.vcf.gz

# Use sed to adjust the header
sed 's/^>46 \(HG992208.2\):[0-9]*-[0-9]*/>\1/' britomartisCO1_${SAMPLE}_f.fasta > britomartisCO1_${SAMPLE}_f_headfix.fasta

# Run bedtools to extract the sequence
bedtools getfasta -fi britomartisCO1_${SAMPLE}_f_headfix.fasta -bed CO1.bed -fo britomartisCO1_${SAMPLE}_f_headfix_final.fasta

#!/bin/bash

# Define the output file
output_file="combined.fasta"

# Ensure the output file is empty to start
> "$output_file"

# Loop over all fasta files with the naming pattern
for fasta in britomartisCO1_*_f_headfix_final.fasta; do
    # Extract ID from the filename
    ID="${fasta#britomartisCO1_}"
    ID="${ID%_f_headfix_final.fasta}"

    # Take only the first half of the redundant ID
    ID=$(echo "$ID" | awk -F'_' '{print $1"_"$2"_"$3}')

    # Use sed to modify the header and append to the output file
    sed "s/^>HG992208.2:[0-9]*-[0-9]*$/>$ID/" "$fasta" >> "$output_file"
done

#Check quality of those:

#Vastmanland_2_1999 no
#Stockholm_4_1965 really low QC
#Gastrikland_13_1969
#(and perhaps Gastrikland_1_1941

#Check SNP quality overall

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J SNP_QC.miss20
#SBATCH --output=SNP_QC.miss20.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools vcftools


vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --freq2 --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --depth --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --site-mean-depth --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --site-quality  --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --missing-indv  --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --missing-site  --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20

vcftools --gzvcf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.miss20.vcf.gz --het --out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/britomartis_all_mcaller.miss20


for sample in $(bcftools query -l britomartis_all_mcaller.COI.vcf.gz); do
  grep './.' ${sample}.no0.vcf.gz
done

for sample in $(bcftools query -l britomartis_all_mcaller.COI.vcf.gz); do zgrep -H '\./\.' ${sample}.no0.vcf.gz ; done

#Quality control performed in python3

#Samples to filter

KALM_4_1957 - diamina
ALTA_2_2015 - diamina
BAJK_4_2016 - diamina
ITAL_1_1946 - diamina
VORO_1_1998 - athalia
VAST_3_2005 - athalia
URAL_1_1991 - athalia
RUSS_2_1999 - athalia
SMAL_6_2013 - athalia
KRAS_4_2002 - athalia
KRAS_2_2002 - athalia
CHEH_3_2008 - athalia
CHEH_2_2004 - athalia
POLA_1_2003 - athalia

#Quality red flags, deph lower then 2.5
ind depth
KALM_7_1961 2.25957
STOC_3_1965 0.960228
GAST_13_1969 1.40341
GAST_8_1965 1.10823
GAST_9_1965 1.31725
SMAL_1_1967 1.97408
STOC_1_1965 0.520048

##############################################
##############################################
###### NEW ERA ##############################

#Workplan
# 1. Remove extra
### a. Analyse with plink
### b. Analyse with ANGSD
### c. Plot PCA
### d. PLot bsaic stats incl. missingness/heterozygosity etc

# 2. Subsampe big Files
### a. Analyse from BAM with ANGSD
### b. Plot

# 3. Assemble mtDNA


#1. Remove low quality and invariant
#Provide files containing a list of individuals to either include or exclude in subsequent analysis. Each individual ID (as defined in the VCF headerline) should be included on a separate line. If both options are used, then the "--keep" option is executed before the "--remove" option. When multiple files are provided, the union of individuals from all keep files subtracted by the union of individuals from all remove files are kept. No header line is expected.

nano keep.reduced.txt
#
URAL_2_1995_URAL_2_1995
RUSS_4_2008_RUSS_4_2008
JAPA_1_1994_JAPA_1_1994
CHEH_1_1989_CHEH_1_1989
KRAS_1_2002_KRAS_1_2002
RUSS_1_1998_RUSS_1_1998
URAL_1_1991_URAL_1_1991
STOC_6_1965_STOC_6_1965
VAST_1_1983_VAST_1_1983
SMAL_4_1996_SMAL_4_1996

vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz --keep keep.reduced.txt --recode --stdout | gzip -c > merge.contemporary.reduced.freebayes.bg.vcf.gz

module load bioinfo-tools vcftools bcftools
bcftools query -l britomartis_all_mcaller.COI.vcf.gz

nano keep.HQ.indviduals.list

ALTA_1_2015_ALTA_1_2015
BAJK_1_2016_BAJK_1_2016
BAJK_2_2016_BAJK_2_2016
BAJK_3_2016_BAJK_3_2016
BELA_1_2004_BELA_1_2004
CHEH_1_1989_CHEH_1_1989
DALA_1_1965_DALA_1_1965
GAST_10_1965_GAST_10_1965
GAST_11_1965_GAST_11_1965
GAST_1_1941_GAST_1_1941
GAST_12_1969_GAST_12_1969
GAST_14_1969_GAST_14_1969
GAST_2_1941_GAST_2_1941
GAST_3_1941_GAST_3_1941
GAST_4_1941_GAST_4_1941
GAST_5_1943_GAST_5_1943
GAST_6_1965_GAST_6_1965
GAST_7_1965_GAST_7_1965
JAPA_1_1994_JAPA_1_1994
JAPA_2_1994_JAPA_2_1994
KALM_10_1983_KALM_10_1983
KALM_1_2018_KALM_1_2018
KALM_2_1955_KALM_2_1955
KALM_3_1955_KALM_3_1955
KALM_5_1961_KALM_5_1961
KALM_6_1961_KALM_6_1961
KALM_8_1969_KALM_8_1969
KALM_9_1981_KALM_9_1981
KAZA_1_2006_KAZA_1_2006
KAZA_2_2006_KAZA_2_2006
KRAS_1_2002_KRAS_1_2002
KRAS_3_2002_KRAS_3_2002
POLA_2_2003_POLA_2_2003
RUSS_1_1998_RUSS_1_1998
RUSS_3_1999_RUSS_3_1999
RUSS_4_2008_RUSS_4_2008
RUSS_5_2008_RUSS_5_2008
SLOV_1_1990_SLOV_1_1990
SMAL_2_1967_SMAL_2_1967
SMAL_3_1967_SMAL_3_1967
SMAL_4_1996_SMAL_4_1996
SMAL_5_1998_SMAL_5_1998
STOC_2_1965_STOC_2_1965
STOC_4_1965_STOC_4_1965
STOC_5_1965_STOC_5_1965
STOC_6_1965_STOC_6_1965
UPPS_1_1951_UPPS_1_1951
UPPS_2_1958_UPPS_2_1958
UPPS_3_1969_UPPS_3_1969
URAL_2_1995_URAL_2_1995
VAST_1_1983_VAST_1_1983
VAST_2_1999_VAST_2_1999

(base) [daria@rackham1 07_mpileupM_joined]$ wc -l keep.HQ.indviduals.list
52 keep.HQ.indviduals.list

vcftools --gzvcf merge.contemporary.freebayes.bg.vcf.gz --keep keep.HQ.indviduals.list --recode --stdout | gzip -c > merge.contemporary.reduced.freebayes.bg.vcf.gz

vcftools --gzvcf britomartis_all_mcaller.vcf.gz --keep keep.HQ.indviduals.list --max-missing-count 0 --remove-indels --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > britomartis_all_mcaller.miss3.vcf.gz

bcftools view -i 'FORMAT/GQ > 30' input.vcf -Ov -o filtered.vcf


#!/bin/sh
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1-00:00:00
#SBATCH -J all_sample_filter
#SBATCH --output=all_sample_filter.out
#SBATCH --mail-user=daria.shipilina@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bcftools samtools vcftools

vcftools --gzvcf britomartis_all_mcaller.vcf.gz --keep keep.HQ.indviduals.list --max-missing-count 0 --remove-indels --min-alleles 2 --max-alleles 2 --minQ 30 --recode --stdout  | gzip -c > britomartis_all_mcaller.HQ.miss0.vcf.gz

tabix britomartis_all_mcaller.HQ.miss0.vcf.gz

bcftools view -i 'FORMAT/GQ >= 30' britomartis_all_mcaller.HQ.miss0.vcf.gz -Oz -o britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz

FORMAT/DP[0-20] > 30


# 3. Assemble mtDNAm

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J mtgenome_assembly
#SBATCH --mail-user=daria.shipilina@gmail.com

module load bioinfo-tools GetOrganelle/1.7.7.0
module load GetOrganelleDB

get_organelle_from_reads.py -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_2.fastp.fastq.gz -t 4 -o SMAL_2_1967.mt -F animal_mt -R 10


/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J mtgenome_assembly
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --array=1-3  # Replace N with the number of lines in samples.txt

module load bioinfo-tools GetOrganelle/1.7.7.0
module load GetOrganelleDB

# Extract the paths from the samples.txt file using the array task ID as the line number
READ1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fasta_pathes_britomartis.txt | cut -f1)
READ2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fasta_pathes_britomartis.txt | cut -f2)

# Set the output directory based on the sample name
OUTPUT_DIR=$(basename $(dirname $READ1))

get_organelle_from_reads.py -1 $READ1 -2 $READ2 -t 4 -o ${OUTPUT_DIR}.mt -F animal_mt -R 10


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J mtgenome_assembly
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --array=1-N  # Replace N with the number of lines in samples.csv (excluding the header)

module load bioinfo-tools GetOrganelle/1.7.7.0
module load GetOrganelleDB

# Skip the header (using tail) and then extract the paths using awk and the SLURM_ARRAY_TASK_ID as the record number
READ1=$(tail -n +2 sample_sheetHQ.csv | awk -F, -v line=${SLURM_ARRAY_TASK_ID} 'NR==line{print $4}')
READ2=$(tail -n +2 sample_sheetHQ.csv | awk -F, -v line=${SLURM_ARRAY_TASK_ID} 'NR==line{print $5}')

# Set the output directory based on the sample name
OUTPUT_DIR=$(tail -n +2 sample_sheetHQ.csv | awk -F, -v line=${SLURM_ARRAY_TASK_ID} 'NR==line{print $2}')

get_organelle_from_reads.py -1 $READ1 -2 $READ2 -t 4 -o ${OUTPUT_DIR}.mt -F animal_mt -R 10

#NOT ENOUGH MEMORY


# 2. Subsampe big Files
### a. Analyse from BAM with ANGSD
### b. Plot

#DownsampleSam function from GATK

Downsample file, keeping about 10% of the reads

 java -jar picard.jar DownsampleSam \
       I=input.bam \
       O=downsampled.bam \
       P=0.1

gatk DownsampleSam I=MUZAV26_94.md.cram  O=downsampled.bam P=0.1

samtools view -b MUZAV26_94.md.cram | gatk DownsampleSam -I /dev/stdin -O downsampled.bam -P 0.1

###################
##  Direction 1  ##
###################

bcftools stats britomartis_all_mcaller.HQ.miss0.vcf.gz | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	52
SN	0	number of records:	584552
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	584552
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0

bcftools stats britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	52
SN	0	number of records:	584552
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	584552
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0

#### Filter fixed alternative

bcftools view -i 'GT!="1/1"' britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz

bcftools view -i 'FORMAT/DP[0-20] > 30' britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz

FORMAT/DP[0-20] > 30

vcftools --gzvcf britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz --minGQ 30


--recode --stdout | bgzip -c > LR9999${i}allsites_hardTEfiltered_noindel_invariant_Qfilter.vcf.gz

'MEAN(FMT/DP)>10'


bcftools view -i 'N_PASS(FORMAT/GQ >= 30) > 51' britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz -Oz -o britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.vcf.gz

bcftools stats britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.vcf.gz | head -n 30 | tail -n 9

tabix britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.vcf.gz

bcftools view -r HG992177.1,HG992178.1,HG992179.1,HG992180.1,HG992181.1,HG992182.1,HG992183.1,HG992184.1,HG992185.1,HG992186.1,HG992187.1,HG992188.1,HG992189.1,HG992190.1,HG992191.1,HG992192.1,HG992193.1,HG992194.1,HG992195.1,HG992196.1,HG992197.1,HG992198.1,HG992199.1,HG992200.1,HG992201.1,HG992202.1,HG992203.1,HG992204.1,HG992205.1,HG992206.1 britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.vcf.gz -Oz -o britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.autosomes.vcf.gz

bcftools stats britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.autosomes.vcf.gz  | head -n 30 | tail -n 9

bcftools view -r HG992208.2 britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.vcf.gz -Oz -o britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.mtDNA.vcf.gz

bcftools stats britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.mtDNA.vcf.gz  | head -n 30 | tail -n 9

###### PCA ######

plink --vcf ../britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out britomartis.all
Pruning complete.  1486 of 2203 variants removed.


plink --vcf ../britomartis.merged.freebayes.qfilter.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out britomartis.all

plink --vcf ../britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.autosomes
19453 variants and 52 people pass filters and QC.

plink --vcf ../britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.mtDNA.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.mtDNA

module load plink


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J SNP_QC.miss20
#SBATCH --output=SNP_QC.miss20.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools vcftools

# Define the input VCF file and output directory
input_vcf="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.vcf.gz"
output_dir="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/07_mpileupM_joined/QC/"

# Create a loop to iterate over each VCF operation with a sample name
vcf_operations=("freq2" "depth" "site-mean-depth" "site-quality" "missing-indv" "missing-site" "het")

for operation in "${vcf_operations[@]}"; do
    output_file="${output_dir}britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.${operation}"
    vcftools --gzvcf "$input_vcf" --"$operation" --out "$output_file"
done



########## Removing

vcftools --gzvcf britomartis_all_mcaller.HQ.miss0.GQfilteredNPASS.autosomes.vcf.gz --keep keep.HQextra.indviduals.list --recode --stdout  | gzip -c > britomartis_all_mcaller.HQ.miss0.GQextra.filteredNPASS.autosomes.vcf.gz

plink --vcf ../britomartis_all_mcaller.HQ.miss0.GQextra.filteredNPASS.autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out britomartis_all_mcaller.HQ.miss0.GQextra.filteredNPASS.autosomes



####### Trying ANGSD

$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL -r 11 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 60 -setMaxDepth 400 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null


ANGSD/0.940-stable

PCAngsd


module load bioinfo-tools ANGSD PCAngsd

./angsd -GL 2 -doGlf 2 -b bam.filelist -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1

angsd -GL 2 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam bam.filelist

bcftools view -i 'N_PASS(FORMAT/GQ >= 30) > 51' britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz

vcftools --gzvcf britomartis_all_mcaller.HQ.miss0.GQfiltered.vcf.gz --keep keep.HQextra.indviduals.list --maf 0.05 --recode --stdout  | gzip -c > britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.vcf.gz

bcftools stats britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.vcf.gz | head -n 30 | tail -n 9

plink --vcf ../britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF

plink --vcf ../britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.prune.in --make-bed --pca --out britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF

britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.eigenval
britomartis_all_mcaller.HQ.miss0.GQfiltered.MAF.eigenvec

http://www.popgen.dk/software/index.php/PCAngsd

Test set of
/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/04_CallingContemporary/00_mpileupCalling


##### TESTING ANGST #####

module load bioinfo-tools ANGSD PCAngsd

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD


angsd -GL 2 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam reduced_sample_ANGSTtest1.list

pcangsd -beagle genolike.beagle.gz -out output -threads 64

### Try angsd for mtDNA

### extract mtDNA


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J mtDNAbam
#SBATCH --output=mtDNAbam.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-3

# Load the samtools module (modify this if needed)
module load bioinfo-tools samtools

# Set the output directory
output_dir=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams

# Read the input file paths from reducedset.info line by line
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/reduced_sample_ANGSTtest1.list)

# Extract the sample name from the input file path
sample_name=$(basename "${input_file}" .md.cram)

# Run samtools view and save the output to a BAM file
samtools view "${input_file}" HG992208.2 -o "${output_dir}/${sample_name}.cram"
tabix "${output_dir}/${sample_name}.cram"


input_file=$(sed -n 1p /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/reduced_sample_ANGSTtest1.list)

samtools view "${input_file}" HG992208.2


### Create bam list

#!/bin/bash

# File containing all file paths
all_files_file="full_sample.list"

# File containing samples to keep
samples_to_keep_file="keep.HQextra.indviduals.list"

# Output file to store the filtered file paths
output_file="filtered_paths_bam.txt"

while IFS= read -r sample_name; do
    # Iterate through each line in the list of all files
    while IFS= read -r file_path; do
        # Extract the sample name from the file path
        file_sample_name=$(basename "$file_path" .md.cram)

        # Check if the sample name matches the current sample to keep
        if [ "$sample_name" == "$file_sample_name" ] || [ "$sample_name" == "${file_sample_name}_${file_sample_name}" ]; then
            # Sample is in the list, so add the file path to the output
            echo "$file_path" >> "$output_file"
        fi
    done < "$all_files_file"
done < "$samples_to_keep_file"

echo "Filtered files are saved in $output_file"

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J mtDNAbam
#SBATCH --output=mtDNAbam.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-3

# Load the samtools module (modify this if needed)
module load bioinfo-tools samtools

# Set the output directory
output_dir=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams

# Read the input file paths from reducedset.info line by line
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/filtered_paths_bam.txt)

# Extract the sample name from the input file path
sample_name=$(basename "${input_file}" .md.cram)

# Run samtools view and save the output to a BAM file
samtools view "${input_file}" HG992208.2 -o "${output_dir}/${sample_name}.cram"
tabix "${output_dir}/${sample_name}.cram"


./BAJK_2_2016.cram
./URAL_2_1995.cram
./RUSS_1_1998.cram
./KAZA_1_2006.cram
./SMAL_2_1967.cram
./UPPS_3_1969.cram
./CHEH_1_1989.cram
./GAST_4_1941.cram
./KALM_3_1955.cram
./UPPS_1_1951.cram
./KRAS_1_2002.cram
./STOC_4_1965.cram
./UPPS_2_1958.cram
./DALA_1_1965.cram
./GAST_10_1965.cram
./VAST_1_1983.cram
./GAST_1_1941.cram
./KALM_10_1983.cram
./GAST_5_1943.cram
./SLOV_1_1990.cram
./GAST_14_1969.cram
./GAST_12_1969.cram
./GAST_3_1941.cram
./STOC_5_1965.cram
./BAJK_1_2016.cram
./KALM_1_2018.cram
./GAST_6_1965.cram
./RUSS_3_1999.cram
./JAPA_1_1994.cram
./KALM_2_1955.cram
./SMAL_3_1967.cram
./SMAL_5_1998.cram
./STOC_2_1965.cram
./KALM_8_1969.cram
./KAZA_2_2006.cram
./KALM_6_1961.cram
./GAST_11_1965.cram
./RUSS_5_2008.cram
./JAPA_2_1994.cram
./POLA_2_2003.cram
./KRAS_3_2002.cram
./SMAL_4_1996.cram
./GAST_2_1941.cram
./KALM_9_1981.cram
./BELA_1_2004.cram
./GAST_7_1965.cram
./RUSS_4_2008.cram
./ALTA_1_2015.cram
./STOC_6_1965.cram
./KALM_5_1961.cram


module load bioinfo-tools ANGSD PCAngsd samtools

angsd -GL 2 -out mtDNA_assmann -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam sample_mtDNA.list

pcangsd -b mtDNA_assmann.beagle.gz -o mtDNA_assmann -threads 64

###Moving to python

### Trying to evenout coverage

samtools depth ALTA_1_2015.cram	| awk '{sum+=$3} END {print "Mean Coverage:", sum/NR}'

### Plan Oct.30
#1. Try subsampling
#2. Run whole genomes
#3. Subtract repeats
#4. Run relatedness
#5. pi whole and within populations, 4-fold
#6. Taj D?

#1. Try subsampling
#!/bin/bash

# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*.cram; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%.cram}"

    # Run samtools depth and calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Print the result
    echo "File: $file_name_no_ext - Mean Coverage: $mean_coverage"
done

File: ALTA_1_2015 - Mean Coverage: 25930.4
File: BAJK_1_2016 - Mean Coverage: 344.696
File: BAJK_2_2016 - Mean Coverage: 2946.04
File: BELA_1_2004 - Mean Coverage: 2473.57
File: CHEH_1_1989 - Mean Coverage: 8737.7
File: DALA_1_1965 - Mean Coverage: 351.153
File: GAST_10_1965 - Mean Coverage: 203.476
File: GAST_11_1965 - Mean Coverage: 411.907
File: GAST_1_1941 - Mean Coverage: 602.453
File: GAST_12_1969 - Mean Coverage: 735.909
File: GAST_14_1969 - Mean Coverage: 298.635
File: GAST_2_1941 - Mean Coverage: 125.333
File: GAST_3_1941 - Mean Coverage: 422.677
File: GAST_4_1941 - Mean Coverage: 203.288
File: GAST_5_1943 - Mean Coverage: 547.975
File: GAST_6_1965 - Mean Coverage: 160.235
File: GAST_7_1965 - Mean Coverage: 298.282
File: JAPA_1_1994 - Mean Coverage: 1572.59
File: JAPA_2_1994 - Mean Coverage: 3239.56
File: KALM_10_1983 - Mean Coverage: 6271.84
File: KALM_1_2018 - Mean Coverage: 26890.8
File: KALM_2_1955 - Mean Coverage: 258.309
File: KALM_3_1955 - Mean Coverage: 394.039
File: KALM_5_1961 - Mean Coverage: 121.643
File: KALM_6_1961 - Mean Coverage: 284.727
File: KALM_8_1969 - Mean Coverage: 259.423
File: KALM_9_1981 - Mean Coverage: 6428.52
File: KAZA_1_2006 - Mean Coverage: 4034.71
File: KAZA_2_2006 - Mean Coverage: 2179.59
File: KRAS_1_2002 - Mean Coverage: 1130.79
File: KRAS_3_2002 - Mean Coverage: 1873.8
File: POLA_2_2003 - Mean Coverage: 1173.58
File: RUSS_1_1998 - Mean Coverage: 2037.39
File: RUSS_3_1999 - Mean Coverage: 6646.78
File: RUSS_4_2008 - Mean Coverage: 6443.06
File: RUSS_5_2008 - Mean Coverage: 887.994
File: SLOV_1_1990 - Mean Coverage: 1636.64
File: SMAL_2_1967 - Mean Coverage: 171.244
File: SMAL_3_1967 - Mean Coverage: 346.479
File: SMAL_4_1996 - Mean Coverage: 4432.66
File: SMAL_5_1998 - Mean Coverage: 5142.24
File: STOC_2_1965 - Mean Coverage: 28.9061
File: STOC_4_1965 - Mean Coverage: 7.51993
File: STOC_5_1965 - Mean Coverage: 15.1845
File: STOC_6_1965 - Mean Coverage: 740.16
File: UPPS_1_1951 - Mean Coverage: 309.581
File: UPPS_2_1958 - Mean Coverage: 185.247
File: UPPS_3_1969 - Mean Coverage: 241.404
File: URAL_2_1995 - Mean Coverage: 10654.1
File: VAST_1_1983 - Mean Coverage: 5219.36

module load GATK

#!/bin/bash

# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Define the target coverage
target_coverage=120

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*.cram; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%.cram}"

    # Run samtools depth to calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Calculate P based on target coverage and mean coverage
    P=$(awk "BEGIN {print $target_coverage / $mean_coverage}")

    # Run gatk DownsampleSam with the calculated P value
    gatk DownsampleSam I="$cram_file" O="${file_name_no_ext}_downsampled.bam" P="$P"
done

# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*_downsampled.bam; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%_downsampled.bam}"

    # Run samtools depth and calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Print the result
    echo "File: $file_name_no_ext - Mean Coverage: $mean_coverage"
done

samtools depth KALM_10_1983_downsampled.bam	| awk '{sum+=$3} END {print "Mean Coverage:", sum/NR}'

samtools view -o ALTA_1_2015.bam ALTA_1_2015.cram

gatk DownsampleSam -I /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams/GAST_2_1941.bam  -O GAST_2_1941_downsampled.bam -P 0.1


samtools depth GAST_2_1941_downsampled.bam | awk '{sum+=$3} END {print "Mean Coverage:", sum/NR}'


samtools view -b -o GAST_2_1941.bam GAST_2_1941.cram




# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Define the target coverage
target_coverage=120

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*.cram; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%.cram}"
    output_bam="${file_name_no_ext}.bam"

    # Run samtools depth to calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Calculate P based on target coverage and mean coverage
    P=$(awk "BEGIN {print $target_coverage / $mean_coverage}")

    # Run gatk DownsampleSam with the calculated P value
    samtools view -b -o "$output_bam" "$cram_file"
    gatk DownsampleSam I="$output_bam" O="${file_name_no_ext}_downsampled.bam" P="$P"
done

# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*_downsampled.bam; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%_downsampled.bam}"

    # Run samtools depth and calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Print the result
    echo "File: $file_name_no_ext - Mean Coverage: $mean_coverage"
done


File: ALTA_1_2015 - Mean Coverage: 118.561
File: BAJK_1_2016 - Mean Coverage: 120.225
File: BAJK_2_2016 - Mean Coverage: 117.391
File: BELA_1_2004 - Mean Coverage: 120.574
File: CHEH_1_1989 - Mean Coverage: 117.255
File: DALA_1_1965 - Mean Coverage: 119.262
File: GAST_10_1965 - Mean Coverage: 119.513
File: GAST_11_1965 - Mean Coverage: 122.198
File: GAST_1_1941 - Mean Coverage: 120.158
File: GAST_12_1969 - Mean Coverage: 118.82
File: GAST_14_1969 - Mean Coverage: 118.789
File: GAST_2_1941 - Mean Coverage: 120.093
File: GAST_3_1941 - Mean Coverage: 119.445
File: GAST_4_1941 - Mean Coverage: 120.871
File: GAST_5_1943 - Mean Coverage: 119.709
File: GAST_6_1965 - Mean Coverage: 119.796
File: GAST_7_1965 - Mean Coverage: 120.239
File: JAPA_1_1994 - Mean Coverage: 118.382
File: JAPA_2_1994 - Mean Coverage: 118.499
File: KALM_10_1983 - Mean Coverage: 119.082
File: KALM_1_2018 - Mean Coverage: 118.512
File: KALM_2_1955 - Mean Coverage: 120.191
File: KALM_3_1955 - Mean Coverage: 120.599
File: KALM_5_1961 - Mean Coverage: 120.045
File: KALM_6_1961 - Mean Coverage: 119.297
File: KALM_8_1969 - Mean Coverage: 118.83
File: KALM_9_1981 - Mean Coverage: 121.514
File: KAZA_1_2006 - Mean Coverage: 120.168
File: KAZA_2_2006 - Mean Coverage: 120.49
File: KRAS_1_2002 - Mean Coverage: 120.863
File: KRAS_3_2002 - Mean Coverage: 117.622
File: POLA_2_2003 - Mean Coverage: 122.002
File: RUSS_1_1998 - Mean Coverage: 122.27
File: RUSS_3_1999 - Mean Coverage: 118.258
File: RUSS_4_2008 - Mean Coverage: 120.459
File: RUSS_5_2008 - Mean Coverage: 120.551
File: SLOV_1_1990 - Mean Coverage: 118.125
File: SMAL_2_1967 - Mean Coverage: 119.942
File: SMAL_3_1967 - Mean Coverage: 119.777
File: SMAL_4_1996 - Mean Coverage: 121.546
File: SMAL_5_1998 - Mean Coverage: 120.729
File: STOC_6_1965 - Mean Coverage: 118.046
File: UPPS_1_1951 - Mean Coverage: 118.786
File: UPPS_2_1958 - Mean Coverage: 119.458
File: UPPS_3_1969 - Mean Coverage: 119.829
File: URAL_2_1995 - Mean Coverage: 121.602
File: VAST_1_1983 - Mean Coverage: 121.229


SMAL_5_1998_downsampled.bam

# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*_downsampled.bam; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%_downsampled.bam}"

    # Run samtools depth and calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Print the result
    echo "$cram_file"
done



angsd -GL 2 -out mtDNA_assmann_downsamp -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam sample_mtDNA_downsamp.list

pcangsd -b mtDNA_assmann_downsamp.beagle.gz -o mtDNA_assmann_downsamp -threads 64


#!/bin/bash

# Define the file containing the list of BAM files
file_list="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/filtered_paths_bam.txt"  # Replace with the actual path to your file

# Iterate through BAM files listed in the file
while IFS= read -r cram_file; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%.cram}"

    # Run samtools depth and calculate mean coverage
    mean_coverage=$(samtools depth "$cram_file" | awk '{sum+=$3} END {print sum/NR}')

    # Print the result
    echo "File: $file_name_no_ext - Mean Coverage: $mean_coverage"
done < "$file_list"

###############
#2. Run whole genomes
###############

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -J mtDNAbam
#SBATCH --output=nomtDNAbam.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-50

# Load the samtools module (modify this if needed)
module load bioinfo-tools samtools

# Set the output directory
output_dir=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams

# Read the input file paths from reducedset.info line by line
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/filtered_paths_bam.txt)

# Extract the sample name from the input file path
sample_name=$(basename "${input_file}" .md.cram)

# Run samtools view and save the output to a BAM file
samtools view -h "${input_file}" HG992176.1 HG992177.1 HG992178.1 HG992179.1 HG992180.1 HG992181.1 HG992182.1 HG992183.1 HG992184.1 HG992185.1 HG992186.1 HG992187.1 HG992188.1 HG992189.1 HG992190.1 HG992191.1 HG992192.1 HG992193.1 HG992194.1 HG992195.1 HG992196.1 HG992197.1 HG992198.1 HG992199.1 HG992200.1 HG992201.1 HG992202.1 HG992203.1 HG992204.1 HG992205.1 HG992206.1 HG992207.1 -o "${output_dir}/${sample_name}.cram"
tabix "${output_dir}/${sample_name}.cram"



samtools view -h ../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_3_1969/UPPS_3_1969.md.cram HG992176.1 HG992177.1 HG992178.1 HG992179.1 HG992180.1 HG992181.1 HG992182.1 HG992183.1 HG992184.1 HG992185.1 HG992186.1 HG992187.1 HG992188.1 HG992189.1 HG992190.1 HG992191.1 HG992192.1 HG992193.1 HG992194.1 HG992195.1 HG992196.1 HG992197.1 HG992198.1 HG992199.1 HG992200.1 HG992201.1 HG992202.1 HG992203.1 HG992204.1 HG992205.1 HG992206.1 HG992207.1 -o autos.cram

for cram_file in /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams/*.cram; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%.cram}"

    # Print the result
    echo "$file_name_no_ext"
done

angsd -GL 2 -out chroms_assmann -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam sample_chroms.list


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 09:00:00
#SBATCH -J angsd
#SBATCH --output=angsd.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams
angsd -GL 2 -out chroms_assmann -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam sample_chroms.list



#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 09:00:00
#SBATCH -J angsd
#SBATCH --output=angsd.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams
angsd -b sample_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out chroms_assmann_saf \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
        -GL 1 -doSaf 1


pcangsd -b chroms_assmann.beagle.gz -o chroms_assmann -t 64

###############
#3. RUN PCAngsd

###############
module load bioinfo-tools ANGSD PCAngsd
pcangsd -b mtDNA_assmann_downsamp.beagle.gz -o mtDNA_assmann_downsamp --inbreedSamples

# Define the folder containing the .cram files
folder_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/mtDNAbams"

# Iterate through .cram files in the folder
for cram_file in "$folder_path"/*_downsampled.bam; do
    # Extract the file name without the path and extension
    file_name=$(basename -- "$cram_file")
    file_name_no_ext="${file_name%_downsampled.bam}"

    # Print the result
    echo "$file_name_no_ext"
done


gzip: chroms_assmann.beagle.gz: unexpected end of file
Loaded 3280971 sites and 47 individuals.
Estimating minor allele frequencies.
EM (MAF) converged at iteration: 41
Number of sites after MAF filtering (0.05): 1623233


tree?



### Figuring out pop stats ###
https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial4_summary_stats/markdowns/summary_stats.md


angsd -b $DIR/$POP'_bams.txt' -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
                -GL 1 -doSaf 1

realSFS print chroms_assmann_saf.saf.idx | less -S

realSFS chroms_assmann_saf.saf.idx > chroms_assmann_saf.sfs


realSFS saf2theta chroms_assmann_saf.saf.idx -sfs chroms_assmann_saf.sfs -outname chroms_assmann_saf

awk -F' ' '{print NF; exit}' chroms_assmann_saf.sfs




############# 10.11.2023 ###############
############# Select all Swedish samples ###############
# 0. Filter repeats
# 1. New angsd
# 2. PCA, theta, taj D, relatedness - Swedish

### Pipeline development

# Set pathes
DATA=$BASEDIR/bam
REF=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna
ANC=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna

module load bioinfo-tools ANGSD PCAngsd NgsRelate
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams

#Call GL and MAF simultaneously
#angsd -bam sample_chroms.list -out chroms_assmann -GL 2 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1

angsd -bam sample_mtDNA.list -out chroms_assmann -ref $REF -anc $ANC -GL 2 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -doSaf 1 -rf non_TEintervals.txt

angsd -bam sweden_chroms.list -out sweden_chroms -ref $REF -anc $ANC -GL 2 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -doSaf 1 -rf non_TEintervals.txt -nThreads 16

-nThreads

-minInd 22 quals to 75% of individuals

? -doGlf 2/5

#cd $BASEDIR
#$ANGSD -b $BASEDIR/sample_lists/PANY_bams.txt -ref $REF -out $BASEDIR/results/PANY \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
#        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
#        -GL 2 -doGlf 2

#Call SAF
#angsd -b sample_chroms.list -ref $REF -anc $ANC -out chroms_assmann_saf \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
#        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
#        -GL 1 -doSaf 1 -rf

#### optional
# Call gnotypes

#cd $BASEDIR
#$ANGSD -b $BASEDIR/sample_lists/PANY_bams.txt -ref $REF -out $BASEDIR/results/PANY \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
#        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
#        -GL 2 -doGlf 4

#-rf for excluding Repeats

# Call SNPs (geno)

#$ANGSD -glf $BASEDIR/results/PANY.glf.gz -fai $REF.fai -nInd 15 -out $BASEDIR/results/PANY \
#    -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1

    # iterate over some cutoffs (you can change these)
#    for PV in 0.05 1e-2 1e-4 1e-6
#    do
#            if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
#            $ANGSD -glf ? -nInd 15 -fai ? -out $BASEDIR/results/PANY.$PV \
#                    -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
#                    -SNP_pval $PV &> /dev/null
#            echo $PV `zcat $BASEDIR/results/PANY.$PV.mafs.gz | tail -n+2 | wc -l`
#    done

# # cd $BASEDIR
# $ANGSD -b $BASEDIR'/sample_lists/ALL_bams.txt' -anc $REF -out $BASEDIR'/results/MME_SNPs' \
# -minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
# -GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1

## Call SAF FOLDED!!!

#angsd -b $DIR/$POP'_bams.txt' -ref $REF -anc $ANC -out Results/$POP \
#                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
#                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
#                -GL 1 -doSaf 1

# stats
# PCA in PCangsd
pcangsd -b chroms_assmann.beagle.gz -o chroms_assmann -t 64
# Fis
pcangsd -b mtDNA_assmann_downsamp.beagle.gz -o mtDNA_assmann_downsamp --inbreedSamples
# TRee
pcangsd -b mtDNA_assmann_downsamp.beagle.gz -o mtDNA_assmann_downsamp --tree -tree_samples XXXXX
# SFS can be filtered here!

realSFS print chroms_assmann_saf.saf.idx | less -S

realSFS chroms_assmann_saf.saf.idx -r HG992176.1:10-100000 -fold 1
realSFS chroms_assmann_saf.saf.idx > chroms_assmann_saf.sfs #(UNFOLDED)
realSFS chroms_assmann_saf.saf.idx -r HG992176.1:10-100000 -fold 1 > chroms_assmann_ch1.folded.sfs
$REALSFS $BASEDIR/results/$POP.folded.saf.idx -fold 1 > $BASEDIR/results/$POP.folded.sfs #(FOLDED)
realSFS saf2theta chroms_assmann_saf.saf.idx -r HG992176.1:10-100000 -sfs chroms_assmann_ch1.folded.sfs -outname chroms_assmann_theta
#remove
thetaStat do_stat chroms_assmann_theta.thetas.idx

#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)   Chr     WinCenter       tW      tP      tF      tH      tL      Tajima  fuf     fud     fayh    zeng    nSites
(0,15876)(3029,99785)(0,99785)  HG992176.1      49892   121.389684      186.970604      159.190415      136.322426      161.646515      1.845950        0.442196        -0.939451       0.272558        0.200808        15876

#remove
#Intraspecific variation#

# reletedness
NGSrelate
./ngsrelate -g angsdput.glf.gz -n 100 -f freq  -O newres

### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
./angsd -b filelist -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate
./ngsrelate  -g angsdput.glf.gz -n 100 -f freq  -O newres


# Output files for python
$BASEDIR/results/JIGA.sfs
cov


#### python code for TE identification

import re

def find_non_masked_regions(fasta_file):
    with open(fasta_file, 'r') as file:
        # Skip the header line
        next(file)
        sequence = file.read().replace('\n', '')

    # Regular expression to find non-masked (capital letter) regions
    non_masked_regions = [(match.start(), match.end()) for match in re.finditer(r'[A-Z]+', sequence)]

    return non_masked_regions

# Replace with your FASTA file path
fasta_file ="/crex/proj/uppstore2017185/cb2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna"
regions = find_non_masked_regions(fasta_file)

# Print the non-masked regions
for start, end in regions:
    print(f'chr:start-end -> {start}-{end}')




    import re

    def find_non_masked_regions(fasta_file):
        regions = []
        with open(fasta_file, 'r') as file:
            sequence = ''
            chrm_id = ''
            for line in file:
                if line.startswith('>'):
                    # Process the previous chromosome
                    if sequence:
                        non_masked_regions = [(match.start(), match.end()) for match in re.finditer(r'[A-Z]+', sequence)]
                        for start, end in non_masked_regions:
                            regions.append(f'{chrm_id}:{start}-{end}')
                    # Update the chromosome ID
                    chrm_id = line.split()[0][1:]  # Remove '>' and split by space, taking the first part
                    sequence = ''
                else:
                    sequence += line.strip()

            # Process the last chromosome
            if sequence:
                non_masked_regions = [(match.start(), match.end()) for match in re.finditer(r'[A-Z]+', sequence)]
                for start, end in non_masked_regions:
                    regions.append(f'{chrm_id}:{start}-{end}')

        return regions

    # Replace with your FASTA file path
    fasta_file = 'path_to_your_fasta_file.fasta'
    regions = find_non_masked_regions(fasta_file)

    # Print the non-masked regions with chromosome IDs
    for region in regions:
        print(region)


grep -v "CAJNA" non_TEintervals.txt > non_TEintervals_chroms.txt

15467.020064 164.788129 10.326439 2.573697 1.640540 1.847954 2.799214 4.721829 7.846736 11.949729 16.107125 19.067260 20.024360 19.040159 16.798602 14.088673 11.469279 9.209387 7.372525 5.919084 4.775954 3.871481 3.148335 2.565112 2.093427 1.713789 1.412001 1.176698 0.998120 0.867786 0.778663 0.725503 0.705160 0.716861 0.762452 0.846628 0.977080 1.164348 1.421004 1.759633 2.189098 2.709005 3.303215 3.934630 4.544537 5.059395 5.405474 2.763829 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000


angsd -GL 2 -out chroms_assmannGLF -nThreads 10 -doGlf 5 -doMajorMinor 1 -SNP_pval 1e-6 -bam sample_mtDNA.list


#Script for angsd

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 22:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams
REF=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna
ANC=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna

angsd -bam sweden_chroms.list -out sweden_chroms_wg -ref $REF -anc $ANC -GL 2 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -doSaf 1 -rf non_TEintervals_chroms.txt -nThreads 16


angsd -b sample_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out chroms_sweden_saf \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
        -GL 1 -doSaf 1


        angsd -b sweden_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out sweden_saf \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 \
                -GL 2 -doSaf 1


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams

angsd -b sweden_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out sweden_2 -nThreads 8 -uniqueOnly 1 -remove_bads 1 -doGlf 2 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 3-00:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams

angsd -b sweden_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test1 -nThreads 8 -uniqueOnly 1 -remove_bads 1 -doGlf 2 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -rf non_TEintervals_chroms.txt -doMajorMinor 1 -doMaf 1


samtools sort SMAL_5_1998.cram -o SMAL_5_1998.sorted.cram
samtools index SMAL_5_1998.sorted.cram

samtools sort STOC_6_1965.cram -o STOC_6_1965.sorted.cram
samtools index STOC_6_1965.sorted.cram


angsd -b sweden_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out sweden_tst -nThreads 16 -uniqueOnly 1 -remove_bads 1 -doGlf 2 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1

angsd -b sweden_chroms.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out sweden_tst -nThreads 16 -doGlf 2 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1

####re sorting

samtools sort STOC_6_1965.cram -o STOC_6_1965.sorted.cram
samtools index STOC_6_1965.sorted.cram


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J mtDNAbam2
#SBATCH --output=mtDNAbam2.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-3

# Load the samtools module (modify this if needed)
module load bioinfo-tools samtools

# Set the output directory
output_dir=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams2

# Read the input file paths from reducedset.info line by line
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams2/sweden_chroms.list)

# Extract the sample name from the input file path
sample_name=$(basename "${input_file}" .md.cram)

# Run samtools view and save the output to a BAM file
samtools sort "${input_file}" -o "${output_dir}/${sample_name}.sorted.cram"
tabix "${output_dir}/${sample_name}.sorted.cram"
samtools index "${output_dir}/${sample_name}.sorted.cram"







#!/bin/bash
#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J mtDNAbam
#SBATCH --output=mtDNAbam.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-3

# Load the samtools module (modify this if needed)
module load bioinfo-tools samtools

# Directory containing your CRAM files
CRAM_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams2"

# Directory to store sorted CRAM files
SORTED_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/nomtDNAcrams2"

# Check if the sorted directory exists, if not, create it
mkdir -p "$SORTED_DIR"

# Loop over each CRAM file in the directory
for cram_file in "$CRAM_DIR"/*.cram; do
    # Extract the base name of the file
    base_name=$(basename "$cram_file" .cram)

    # Define the output file path
    sorted_file="$SORTED_DIR/${base_name}_sorted.cram"

    # Sort the CRAM file and output the sorted version
    samtools sort "$cram_file" -o "$sorted_file"

    # Index the sorted CRAM file
    samtools index "$sorted_file"
done

# echo "Sorting complete."


DALA_1_1965.cram.sorted.cram
GAST_10_1965.cram.sorted.cram
GAST_11_1965.cram.sorted.cram

nano small.list


angsd -b small.list -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out sweden_tst -nThreads 8 -doGlf 2 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1



samtools view -h ../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_3_1969/UPPS_3_1969.md.cram HG992176.1 HG992177.1 HG992178.1 HG992179.1 HG992180.1 HG992181.1 HG992182.1 HG992183.1 HG992184.1 HG992185.1 HG992186.1 HG992187.1 HG992188.1 HG992189.1 HG992190.1 HG992191.1 HG992192.1 HG992193.1 HG992194.1 HG992195.1 HG992196.1 HG992197.1 HG992198.1 HG992199.1 HG992200.1 HG992201.1 HG992202.1 HG992203.1 HG992204.1 HG992205.1 HG992206.1 HG992207.1 -o autos.cram


angsd -b filtered_paths_bam.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test4 -nThreads 8 -doGlf 2 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50

angsd -b filtered_paths_bam.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test6 -nThreads 8 -doGlf 2 -minMapQ 10 -minQ 10 -minInd 40 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 1 -C 50

angsd -b filtered_paths_bam.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test8 -nThreads 8 -doGlf 2 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 10 -C 50 -minmaf 0.05

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test9 -nThreads 8 -doGlf 2 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doCounts 1 -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minmaf 0.05


angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test11 -nThreads 8 -doGlf 2 -minMapQ 20 -minQ 20 -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test12 -nThreads 8 -doGlf 2  -minInd 22 -setMinDepth 5 -setMaxDepth 30 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test12 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test13 -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test14 -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10


angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test15 -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf Zonly.list

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out test15 -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf non_TEintervals_Zchrom.txt


@SQ	SN:HG992177.1	LN:25133022	M5:7fa9b50c82266c962fcd75b3ab08b30c	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992178.1	LN:24877037	M5:27fc270259ff632e4114c8ee4d60f384	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992179.1	LN:23648884	M5:a7a4e8ed6853fe54a3ff85ca5f3cb7b8	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992180.1	LN:22925277	M5:7fc9614edb8918ddf026d8c2119d4883	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992181.1	LN:22895925	M5:e02dd9a1b420843eaddbc9b90858895f	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992182.1	LN:22794856	M5:cd5c4fd792d9457a0fe3cf3be6f11208	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992183.1	LN:21872383	M5:1cdd7e22917795ba24790d65d3d99c5f	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992184.1	LN:21419982	M5:be2188058d87d8265e2a3eb2e6d3d91b	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992185.1	LN:21389510	M5:4c7b8257935f70c1aeaa94adecf30d0b	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992186.1	LN:21375889	M5:0e8e7f54ce0625f02b1d85f1769bb86d	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992187.1	LN:21234190	M5:f1269e6e3b70a3bafe1a7debd483d8f5	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992188.1	LN:20508949	M5:8e58e479979be88eadab8a952ee96ebe	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992189.1	LN:20295254	M5:1ad15d58348d00e9c8806a3b1b58f7f8	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992190.1	LN:20209291	M5:f95845f9d27a5df6de0ca31100bc07d9	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992191.1	LN:19986373	M5:cecb217ccbb3d97b50680374e635665f	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992192.1	LN:19818840	M5:40197286d65834b4bbffc54664657741	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992193.1	LN:19640724	M5:87d42ddd8e2e20b55d02f34c5c81fee4	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992194.1	LN:19553016	M5:4e70f810261f69e1b51ff7d0cf063683	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992195.1	LN:18440301	M5:9fcb2abc5ed1c5fc4b0eb29c81a81648	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992196.1	LN:18366132	M5:25965eba3b45ccaba37ce480485380b9	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992197.1	LN:17060131	M5:51aaea1eba59d8d0e3e8dbc5218e21d8	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992198.1	LN:16621731	M5:621d368e2b1cf954d92e778a81789c85	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992199.1	LN:15218959	M5:276b31dccd5a7848c78793414456976d	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992200.1	LN:15145429	M5:7e282a201e17f87afa0bd0d4387d2e59	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992201.1	LN:14968207	M5:1deb4206c49c7dc1a5fe832329e79259	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992202.1	LN:14402919	M5:b0928eb4c70d977865034d047700251a	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992203.1	LN:13271753	M5:066862bbb0df85b62c4dd0716e5ed7b3	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992204.1	LN:12763386	M5:99ad602e114ce45f628d0dcd07c9a161	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992205.1	LN:11974655	M5:d6fffbfb33363fda2f413d9b578994be	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992206.1	LN:10898131	M5:dbc784c16f8035f1c2646e3fe1532dfb	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992207.1	LN:5265952	M5:9a5f6fa634bcb24209afe3e2b217b6bd	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna
@SQ	SN:HG992176.1	LN:26233870	M5:9e4d67f3fd26622b14b38dc272cc36f1	UR:/scratch/38675654/nxf.1NNvBKX4OK/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.fna



chromo  position        major   minor   ref     anc     knownEM nInd
HG992176.1      159212  G       A       T       T       0.315583        11
HG992176.1      159228  G       T       G       G       0.074429        11
HG992176.1      159240  T       C       T       T       0.176845        11
HG992176.1      213355  G       T       G       G       0.073204        13
HG992176.1      703170  T       G       T       T       0.063740        10
HG992176.1      1147217 A       T       A       A       0.127632        10
HG992176.1      1147218 A       G       A       A       0.063367        10
HG992176.1      1620893 T       A       T       T       0.093377        10
HG992176.1      1746046 T       A       T       T       0.074459        10
HG992176.1      2292842 A       T       A       A       0.068497        10
HG992176.1      2292847 T       C       T       T       0.073135        10
HG992176.1      2580077 C       T       C       C       0.161427        10
HG992176.1      3031387 T       G       T       T       0.082412        12
HG992176.1      3031401 T       A       T       T       0.261492        14
HG992176.1      3031404 G       C       G       G       0.065952        16
HG992176.1      3031406 C       G       C       C       0.148772        16
HG992176.1      3031413 T       C       C       C       0.314479        12
HG992176.1      3099658 G       T       G       G       0.075900        10
HG992176.1      3536838 A       C       A       A       0.073715        10
HG992176.1      5076861 A       G       A       A       0.062530        10
HG992176.1      5076876 T       A       T       T       0.063995        10
HG992176.1      5939552 T       C       T       T       0.139345        10
HG992176.1      5939555 C       T       C       C       0.141604        10
HG992176.1      7306887 A       T       A       A       0.056553        14
HG992176.1      7454322 C       T       C       C       0.074877        10
HG992176.1      7673135 T       C       T       T       0.065881        10
HG992176.1      7673159 T       A       T       T       0.073161        10
HG992176.1      8094648 T       A       T       T       0.060473        10
HG992176.1      8094649 A       G       A       A       0.061028        10
HG992176.1      8423143 C       T       C       C       0.068921        10
HG992176.1      8423146 C       A       C       C       0.306269        10
HG992176.1      8423158 A       G       A       A       0.062721        10
HG992176.1      8423167 G       A       G       G       0.066819        10
HG992176.1      8423172 C       G       C       C       0.065584        10
HG992176.1      8423174 T       C       T       T       0.066301        10
HG992176.1      8423176 G       A       G       G       0.065224        10
HG992176.1      8423178 T       A       T       T       0.065069        10
HG992176.1      8423181 G       A       G       G       0.062353        10
HG992176.1      8423187 G       T       G       G       0.060473        10
HG992176.1      8423194 C       A       C       C       0.064400        10
HG992176.1      8423204 C       T       C       C       0.063034        10
HG992176.1      8423208 T       A       T       T       0.343978        10
HG992176.1      8423209 A       T       A       A       0.062953        10
HG992176.1      9257559 G       T       G       G       0.063761        11
HG992176.1      9454038 A       T       A       A       0.067501        10
HG992176.1      9514081 A       G       A       A       0.052580        10
HG992176.1      9514082 A       G       A       A       0.246404        10



#!/bin/bash

# Path to the file containing the intervals
FILE_PATH="non_TEintervals_Zchrom.txt"
total_length=0
while read -r line; do
    start=$(echo $line | cut -d ':' -f 2 | cut -d '-' -f 1)
    end=$(echo $line | cut -d ':' -f 2 | cut -d '-' -f 2)
    length=$((end - start + 1))
    total_length=$((total_length + length))
done < "$FILE_PATH"
echo "Total length of regions: $total_length"



#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out angsd_new/sweden -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf coordinates_noZ_nomt.list

angsd -b filtered_paths_bam.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out angsd_new/all -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf coordinates_noZ_nomt.list



## Prepare a geno file by subsampling one SNP in every 50 SNPs in the beagle file
zcat sweden.beagle.gz  | awk 'NR % 50 == 0' | cut -f 4- | gzip  > sweden.LDprep.beagle.gz

## Prepare a pos file by subsampling one SNP in every 50 SNPs in the mafs filre
zcat sweden.mafs.gz | cut -f 1,2 |  awk 'NR % 50 == 0' | sed 's/:/_/g'| gzip > sweden.LDprep.pos.gz

 module load bioinfo-tools ngsLD


 $NGSLD/ngsLD \
--geno $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.beagle.gz \
--pos $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.pos.gz \
--probs \
--n_ind 60 \
--n_sites 1134 \
--max_kb_dist 0 \
--n_threads 1 \
--out $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld


ngsLD \
--geno sweden.LDprep.test.beagle.gz \
--pos sweden.LDprep.pos.gz \
--probs \
--n_ind 29 \
--n_sites 109 \
--max_kb_dist 0 \
--n_threads 1 \
--out sweden.LDprep.ld

ngsLD \
--geno sweden.LDprep.beagle \
--pos sweden.LDprep.pos \
--probs \
--n_ind 29 \
--n_sites 109 \
--max_kb_dist 10000 \
--n_threads 1 \
--out sweden.LDprep.10k.ld

cp sweden.beagle.gz sweden.test.beagle.gz
gunzip sweden.test.beagle.gz

zcat sweden.test.beagle.gz  | awk 'NR % 50 == 0' | cut -f 4- | gzip  > sweden.LDprep.test.beagle.gz

less sweden.LDprep.test.beagle.gz



ngsLD \
--geno sweden.LDprep.test.beagle.gz \
--pos sweden.LDprep.pos.gz \
--probs \
--n_ind 29 \
--n_sites 109 \
--max_kb_dist 10000 \
--n_threads 1 \
--out sweden.10k.LDprep.ld


perl /sw/bioinfo/ngsLD/1.1.1/rackham/scripts/prune_graph.pl \
--in_file $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld \
--max_kb_dist 2000 \
--min_weight 0.5 \
--out $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled_unlinked.id

q
ngsLD \
--geno sweden.LDprep.test.beagle.gz \
--pos sweden.LDprep.pos.gz \
--probs \
--n_ind 29 \
--n_sites 109 \
--max_kb_dist 1 \
--n_threads 1 \
--out sweden.LDprep.1kb.ld


python prune_ngsLD.py
usage: prune_ngsLD.py --input sweden.LDprep.ld --output sweden.LDprep.id --max_dist 2000 --min_weight 0.5 --print_excl excel

python prune_ngsLD.py --input sweden.LDprep.ld --output sweden.LDprep.id --max_dist 2000 --min_weight 0.5

conda activate snakemake-tutorial

python ngsLD/scripts/prune_ngsLD.py --input sweden.LDprep.ld --output sweden.LDprep.id --max_dist 2000 --min_weight 0.5

./ngsLD --geno sweden.LDprep.test.beagle.gz --pos sweden.LDprep.pos.gz --probs --n_ind 29 --n_sites 109 --max_kb_dist 0 --n_threads 1 --out sweden.LDprep.ld


#!/bin/bash

# Path to your input file
FILE="sweden.teest.mafs.gz"

# Skip the header and read the file line by line
awk 'NR>1 {
    # Extract the position from the current line
    current_position = $2

    # For the first line, just store the position
    if (NR == 2) {
        previous_position = current_position
    } else {
        # Calculate the difference
        difference = current_position - previous_position

        # Check if the difference is greater than 500
        if (difference > 500) {
            print "Line " NR-1 " (position " previous_position ") and Line " NR " (position " current_position ") are more than 500 bp apart."
        }

        # Update the previous_position for the next iteration
        previous_position = current_position
    }
}' "$FILE"


#!/bin/bash

# Path to your gzipped input file
FILE="sweden.teest.mafs.gz"

# Initialize the previous position variable
previous_position=0

# Read the gzipped file line by line, starting from the second line
zcat "$FILE" | tail -n +2 | while read line; do
    # Extract the position from the current line
    current_position=$(echo $line | awk '{print $2}')

    # For the first line processed, just store the position
    if [ $previous_position -eq 0 ]; then
        previous_position=$current_position
        continue
    fi

    # Calculate the difference
    difference=$((current_position - previous_position))

    # Check if the difference is greater than 500
    if [ $difference -gt 500 ]; then
        echo "Positions $previous_position and $current_position are more than 500 bp apart."
    fi

    # Update the previous_position for the next iteration
    previous_position=$current_position
done



#!/bin/bash

# Path to your gzipped input file
FILE="sweden.teest.mafs.gz"

# Initialize the previous position variable
previous_position=0

# Read the gzipped file line by line, starting from the second line
zcat "$FILE" | tail -n +2 | while read line; do
    # Extract the position from the current line
    current_position=$(echo $line | awk '{print $2}')

    # For the first line processed, just store the position
    if [ $previous_position -eq 0 ]; then
        previous_position=$current_position
        continue
    fi

    # Calculate the difference
    difference=$((current_position - previous_position))

    # Check if the difference is less than 100
    if [ $difference -lt 100 ]; then
        echo "Positions $previous_position and $current_position are less than 100 bp apart."
    fi

    # Update the previous_position for the next iteration
    previous_position=$current_position
done


#!/bin/bash

# Path to your gzipped input file
FILE="sweden.teest.mafs.gz"

# Initialize variables
previous_position=0
line_number=0

# Read the gzipped file line by line, starting from the second line
zcat "$FILE" | tail -n +2 | while read line; do
    # Increment line number
    ((line_number++))

    # Extract the position from the current line
    current_position=$(echo $line | awk '{print $2}')

    # Skip the first line
    if [ $previous_position -eq 0 ]; then
        previous_position=$current_position
        previous_line_number=$line_number
        continue
    fi

    # Calculate the difference
    difference=$((current_position - previous_position))

    # Check if the difference is less than 100
    if [ $difference -lt 100 ]; then
        # Update the previous_position and continue to the next line
        previous_position=$current_position
        previous_line_number=$line_number
    else
        # Print the message and update the previous_position
        echo "Positions at lines $previous_line_number and $line_number are more than 100 bp apart."
        previous_position=$current_position
        previous_line_number=$line_number
    fi
done



#!/bin/bash

# Path to your gzipped input file
FILE="sweden.teest.mafs.gz"

# Initialize variables
previous_position=0
line_number=0

# Read the gzipped file line by line, starting from the second line
zcat "$FILE" | tail -n +2 | while read line; do
    # Increment line number
    ((line_number++))

    # Extract the position from the current line
    current_position=$(echo $line | awk '{print $2}')

    # Skip the first line
    if [ $previous_position -eq 0 ]; then
        previous_position=$current_position
        previous_line_number=$line_number
        continue
    fi

    # Calculate the difference
    difference=$((current_position - previous_position))

    # If the difference is less than 100, update previous_position and continue to next line
    while [ $difference -lt 100 ]; do
        # Read the next line
        if ! read line; then
            # Exit if no more lines
            break 2
        fi
        ((line_number++))
        current_position=$(echo $line | awk '{print $2}')
        difference=$((current_position - previous_position))
    done

    # Print the message for the pair that is 100 bp or more apart
    echo "Positions at lines $previous_line_number and $line_number are 100 bp or more apart."
    previous_position=$current_position
    previous_line_number=$line_number
done

pcangsd -b sweden2.beagle.gz -o sweden2.test -t 64

../01_MappingAll/results/preprocessing/markduplicates/KALM_10_1983/KALM_10_1983.md.cram
../01_MappingAll/results/preprocessing/markduplicates/KALM_1_2018/KALM_1_2018.md.cram 1,2,3,
../01_MappingAll/results/preprocessing/markduplicates/KALM_9_1981/KALM_9_1981.md.cram 4,5,6
../01_MappingAll/results/preprocessing/markduplicates/SMAL_4_1996/SMAL_4_1996.md.cram 7,8,9
../01_MappingAll/results/preprocessing/markduplicates/SMAL_5_1998/SMAL_5_1998.md.cram 10,11,12
../01_MappingAll/results/preprocessing/markduplicates/STOC_6_1965/STOC_6_1965.md.cram 13,14,15
../01_MappingAll/results/preprocessing/markduplicates/VAST_1_1983/VAST_1_1983.md.cram 16,17,18
../03_MappingHistorical/results/preprocessing/markduplicates/DALA_1_1965/DALA_1_1965.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_10_1965/GAST_10_1965.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_11_1965/GAST_11_1965.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_1_1941/GAST_1_1941.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_12_1969/GAST_12_1969.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_14_1969/GAST_14_1969.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_2_1941/GAST_2_1941.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_3_1941/GAST_3_1941.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_4_1941/GAST_4_1941.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_5_1943/GAST_5_1943.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_6_1965/GAST_6_1965.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/GAST_7_1965/GAST_7_1965.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/KALM_2_1955/KALM_2_1955.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/KALM_3_1955/KALM_3_1955.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/KALM_5_1961/KALM_5_1961.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/KALM_6_1961/KALM_6_1961.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/KALM_8_1969/KALM_8_1969.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_2_1967/SMAL_2_1967.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/SMAL_3_1967/SMAL_3_1967.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_1_1951/UPPS_1_1951.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_2_1958/UPPS_2_1958.md.cram
../03_MappingHistorical/results/preprocessing/markduplicates/UPPS_3_1969/UPPS_3_1969.md.cram


realSFS print sweden.saf.idx -fold 1 # | less -S
realSFS sweden.saf.idx -fold 1 > sweden.sfs #(UNFOLDED
realSFS saf2theta sweden.saf.idx -sfs sweden.sfs -outname sweden_theta
#remove
thetaStat do_stat chroms_assmann_theta.thetas.idx


pcangsd -b all.beagle.gz -o all.assmann -t 64


# Remove sample vastmanland

Example: zcat original_beagle.beagle.gz | cut -f10,11,12,25,26,27 --complement |  gzip > edited_beagle.beagle.gz

zcat sweden2.beagle.gz | cut -f22,23,24 --complement |  gzip > sweden2.edited_beagle.noVAST.beagle.gz

gunzip -c sweden2.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Ind%d: %d out of %d\n", (i-1)/3, sum, maxCount}}'

Cols 4-6 (Ind1): 880 out of 5937
Cols 7-9 (Ind2): 855 out of 5937
Cols 10-12 (Ind3): 590 out of 5937
Cols 13-15 (Ind4): 569 out of 5937
Cols 16-18 (Ind5): 2557 out of 5937
Cols 19-21 (Ind6): 1288 out of 5937
Cols 22-24 (Ind7): 460 out of 5937
Cols 25-27 (Ind8): 5817 out of 5937
Cols 28-30 (Ind9): 5892 out of 5937
Cols 31-33 (Ind10): 5533 out of 5937
Cols 34-36 (Ind11): 5870 out of 5937
Cols 37-39 (Ind12): 3113 out of 5937
Cols 40-42 (Ind13): 4326 out of 5937
Cols 43-45 (Ind14): 5891 out of 5937
Cols 46-48 (Ind15): 5087 out of 5937
Cols 49-51 (Ind16): 5045 out of 5937
Cols 52-54 (Ind17): 5668 out of 5937
Cols 55-57 (Ind18): 5749 out of 5937
Cols 58-60 (Ind19): 5885 out of 5937
Cols 61-63 (Ind20): 3448 out of 5937
Cols 64-66 (Ind21): 5590 out of 5937
Cols 67-69 (Ind22): 3999 out of 5937
Cols 70-72 (Ind23): 5820 out of 5937
Cols 73-75 (Ind24): 3023 out of 5937
Cols 76-78 (Ind25): 3931 out of 5937
Cols 79-81 (Ind26): 5125 out of 5937
Cols 82-84 (Ind27): 5892 out of 5937
Cols 85-87 (Ind28): 5260 out of 5937
Cols 88-90 (Ind29): 5681 out of 5937

gunzip -c sweden2.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; print "Row " NR ": " count " occurrences"}'

gunzip -c sweden2.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Ind%d: %d out of %d\n", (i-4)/3, sum, maxCount}}'


gunzip -c sweden2.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-1)/3, sum, maxCount}}'

gunzip -c sweden.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-1)/3, sum, maxCount}}'

gunzip -c sweden.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-4)/3, sum, maxCount}}'

Cols 4-6 (Ind0): 2540 out of 16410
Cols 7-9 (Ind1): 2439 out of 16410
Cols 10-12 (Ind2): 1724 out of 16410
Cols 13-15 (Ind3): 1697 out of 16410
Cols 16-18 (Ind4): 7274 out of 16410
Cols 19-21 (Ind5): 3561 out of 16410
Cols 22-24 (Ind6): 1331 out of 16410
Cols 25-27 (Ind7): 16112 out of 16410
Cols 28-30 (Ind8): 16276 out of 16410
Cols 31-33 (Ind9): 15147 out of 16410
Cols 34-36 (Ind10): 16201 out of 16410
Cols 37-39 (Ind11): 8433 out of 16410
Cols 40-42 (Ind12): 11793 out of 16410
Cols 43-45 (Ind13): 16291 out of 16410
Cols 46-48 (Ind14): 13994 out of 16410
Cols 49-51 (Ind15): 13992 out of 16410
Cols 52-54 (Ind16): 15603 out of 16410
Cols 55-57 (Ind17): 15863 out of 16410
Cols 58-60 (Ind18): 16279 out of 16410
Cols 61-63 (Ind19): 9745 out of 16410
Cols 64-66 (Ind20): 15416 out of 16410
Cols 67-69 (Ind21): 10541 out of 16410
Cols 70-72 (Ind22): 16093 out of 16410
Cols 73-75 (Ind23): 8215 out of 16410
Cols 76-78 (Ind24): 10825 out of 16410
Cols 79-81 (Ind25): 14125 out of 16410
Cols 82-84 (Ind26): 16272 out of 16410
Cols 85-87 (Ind27): 14593 out of 16410
Cols 88-90 (Ind28): 15525 out of 16410

zcat sweden2.beagle.gz | cut -f25-36,43-45,52-60,64-72,70-90 --complement | gzip > edited_beagle.beagle.gz

gunzip -c edited_beagle.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-4)/3, sum, maxCount}}'

zcat sweden2.beagle.gz | cut -f25-36,43-45,52-60,64-72,76-90 --complement | gzip > edited_beagle2.beagle.gz

gunzip -c edited_beagle.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; print "Row " NR ": " count " occurrences"}'

gunzip -c edited_beagle.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; print "Row " NR ": " count / 3 " occurrences"}'

zcat sweden.beagle.gz | cut -f25-36,43-45,52-60,64-72,76-90 --complement | gzip > edited_beagle2.beagle.gz

gunzip -c edited_beagle2.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-4)/3, sum, maxCount}}'

zcat edited_beagle2.beagle.gz | cut -f28-36 --complement | gzip > edited_beagle3.beagle.gz
gunzip -c edited_beagle3.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; print "Row " NR ": " count / 3 " occurrences"}'

gunzip -c edited_beagle3.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; if (count / 3 > 4) linesWithHighCount++; print "Row " NR ": " count / 3 " occurrences"} END {print "Number of lines where count/3 > 4: " linesWithHighCount}'

gunzip -c edited_beagle3.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; if (count / 3 > 4) print "Row " NR " meets criteria"}'

gunzip -c edited_beagle3.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; if (count / 3 <= 4) print $0}' | gzip > edited_beagle3_filtered.beagle.gz

gunzip -c edited_beagle3_filtered.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; if (count / 3 > 4) linesWithHighCount++; print "Row " NR ": " count / 3 " occurrences"} END {print "Number of lines where count/3 > 4: " linesWithHighCount}'



module load bioinfo-tools ANGSD PCAngsd

pcangsd -b edited_beagle3_filtered.beagle.gz -o edited_beagle3_filtered.assmann -t 64

Ind0
Ind1
Ind2
Ind3
Ind4
Ind5
Ind6 remove!!!
Ind11
Ind19
Ind23

zcat edited_beagle2.beagle.gz | cut -f28-36 --complement | gzip > edited_beagle3.beagle.gz

gunzip -c edited_beagle3_filtered.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-4)/3, sum, maxCount}}'

zcat edited_beagle3_filtered.beagle.gz | cut -f22-24 --complement | gzip > edited_beagle3_filtered.beagle.noVAST.gz

pcangsd -b edited_beagle3_filtered.beagle.noVAST.gz -o edited_beagle3_filtered.beagle.noVAST.assmann -t 64

Loaded 4911 sites and 9 individuals.
Estimating minor allele frequencies.
EM (MAF) converged at iteration: 11
Number of sites after MAF filtering (0.05): 3875

edited_beagle3_filtered.beagle.noVAST.assmann.cov


gunzip -c all.beagle.gz | awk '{totalLines++; for (i=4; i<=NF; i++) if ($i == "0.333333") count[i]++} END {maxCount = totalLines * 3; for (i=4; i<=NF; i+=3) {sum=count[i]+count[i+1]+count[i+2]; printf "Cols %d-%d (Ind%d): %d out of %d\n", i, i+2, (i-4)/3, sum, maxCount}}'

gunzip -c all.beagle.gz | awk '{count=0; for (i=4; i<=NF; i++) if ($i == "0.333333") count++; if (count / 3 > 4) linesWithHighCount++; print "Row " NR ": " count / 3 " occurrences"} END {print "Number of lines where count/3 > 4: " linesWithHighCount}'

Cols 4-6 (Ind0): 3870132 out of 5097084
Cols 7-9 (Ind1): 4831501 out of 5097084
Cols 10-12 (Ind2): 3778823 out of 5097084
Cols 13-15 (Ind3): 4067540 out of 5097084
Cols 16-18 (Ind4): 1815833 out of 5097084
Cols 19-21 (Ind5): 1824075 out of 5097084
Cols 22-24 (Ind6): 4402292 out of 5097084
Cols 25-27 (Ind7): 3136955 out of 5097084
Cols 28-30 (Ind8): 2819603 out of 5097084
Cols 31-33 (Ind9): 2453630 out of 5097084
Cols 34-36 (Ind10): 4786209 out of 5097084
Cols 37-39 (Ind11): 2598301 out of 5097084
Cols 40-42 (Ind12): 3733653 out of 5097084
Cols 43-45 (Ind13): 1952516 out of 5097084
Cols 46-48 (Ind14): 2578372 out of 5097084
Cols 49-51 (Ind15): 2180108 out of 5097084
Cols 52-54 (Ind16): 2197163 out of 5097084
Cols 55-57 (Ind17): 3677403 out of 5097084
Cols 58-60 (Ind18): 3741993 out of 5097084
Cols 61-63 (Ind19): 1796207 out of 5097084
Cols 64-66 (Ind20): 2930181 out of 5097084
Cols 67-69 (Ind21): 4701370 out of 5097084
Cols 70-72 (Ind22): 3684989 out of 5097084
Cols 73-75 (Ind23): 3597823 out of 5097084
Cols 76-78 (Ind24): 1924572 out of 5097084
Cols 79-81 (Ind25): 5096552 out of 5097084
Cols 82-84 (Ind26): 5096890 out of 5097084
Cols 85-87 (Ind27): 5088540 out of 5097084
Cols 88-90 (Ind28): 5096862 out of 5097084
Cols 91-93 (Ind29): 4912544 out of 5097084
Cols 94-96 (Ind30): 5032624 out of 5097084
Cols 97-99 (Ind31): 5096847 out of 5097084
Cols 100-102 (Ind32): 5072431 out of 5097084
Cols 103-105 (Ind33): 5071766 out of 5097084
Cols 106-108 (Ind34): 5093450 out of 5097084
Cols 109-111 (Ind35): 5095530 out of 5097084
Cols 112-114 (Ind36): 5096836 out of 5097084
Cols 115-117 (Ind37): 4957762 out of 5097084
Cols 118-120 (Ind38): 5088810 out of 5097084
Cols 121-123 (Ind39): 5007696 out of 5097084
Cols 124-126 (Ind40): 5096546 out of 5097084
Cols 127-129 (Ind41): 4935454 out of 5097084
Cols 130-132 (Ind42): 5005724 out of 5097084
Cols 133-135 (Ind43): 5067015 out of 5097084
Cols 136-138 (Ind44): 5021640 out of 5097084
Cols 139-141 (Ind45): 5038827 out of 5097084
Cols 142-144 (Ind46): 5015808 out of 5097084
Cols 145-147 (Ind47): 5096918 out of 5097084
Cols 148-150 (Ind48): 5077061 out of 5097084
Cols 151-153 (Ind49): 5092003 out of 5097084


Cols 4-6 (Ind0): 3870132 out of 5097084
Cols 7-9 (Ind1): 4831501 out of 5097084
Cols 10-12 (Ind2): 3778823 out of 5097084
Cols 13-15 (Ind3): 4067540 out of 5097084
Cols 16-18 (Ind4): 1815833 out of 5097084
Cols 19-21 (Ind5): 1824075 out of 5097084
Cols 22-24 (Ind6): 4402292 out of 5097084
Cols 25-27 (Ind7): 3136955 out of 5097084
Cols 28-30 (Ind8): 2819603 out of 5097084
Cols 31-33 (Ind9): 2453630 out of 5097084
Cols 34-36 (Ind10): 4786209 out of 5097084
Cols 37-39 (Ind11): 2598301 out of 5097084
Cols 40-42 (Ind12): 3733653 out of 5097084
Cols 43-45 (Ind13): 1952516 out of 5097084
Cols 46-48 (Ind14): 2578372 out of 5097084
Cols 49-51 (Ind15): 2180108 out of 5097084
Cols 52-54 (Ind16): 2197163 out of 5097084
Cols 55-57 (Ind17): 3677403 out of 5097084
Cols 58-60 (Ind18): 3741993 out of 5097084
Cols 61-63 (Ind19): 1796207 out of 5097084
Cols 64-66 (Ind20): 2930181 out of 5097084
Cols 67-69 (Ind21): 4701370 out of 5097084
Cols 70-72 (Ind22): 3684989 out of 5097084
Cols 73-75 (Ind23): 3597823 out of 5097084
Cols 76-78 (Ind24): 1924572 out of 5097084
Cols 91-93 (Ind29): 4912544 out of 5097084
Cols 115-117 (Ind37): 4957762 out of 5097084
Cols 127-129 (Ind41): 4935454 out of 5097084


######################
# BAM QUALITY CHECKS #
######################

HG992177.1_13411770     1       0
HG992177.1_13939774     2       0
HG992177.1_14019323     3       1


HG992177.1  13411760  14019340

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools ANGSD PCAngsd
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD

angsd -b filtered_paths_bam2.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out angsd_new/sweden -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf coordinates_noZ_nomt.list

angsd -b filtered_paths_bam.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna -out angsd_new/all -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf coordinates_noZ_nomt.list
angsd_sweden.sh (END)

samtools view "${input_file}" HG992208.2


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J mtDNAbam
#SBATCH --output=mtDNAbam.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-3

# Load the samtools module (modify this if needed)
module load bioinfo-tools samtools

# Set the output directory
output_dir=/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/visual_test

# Read the input file paths from reducedset.info line by line
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/filtered_paths_bam2.txt )

# Extract the sample name from the input file path
sample_name=$(basename "${input_file}" .md.cram)

# Run samtools view and save the output to a BAM file
samtools view "${input_file}" HG992177.1:13411760-14019340 -o "${output_dir}/${sample_name}.bam"
tabix "${output_dir}/${sample_name}.bam"


  samtools view ../../01_MappingAll/results/preprocessing/markduplicates/KALM_10_1983/KALM_10_1983.md.cram HG992177.1:13411760-14019340 -o KALM_10_1983.sect.bam

../01_MappingAll/results/preprocessing/markduplicates/KALM_10_1983/KALM_10_1983.md.cram



#!/bin/bash

# Set the output directory
output_dir="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/visual_test"

# Path to the file containing input file paths
input_file_path="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/filtered_paths_bam2.txt"

# Iterate over each line in the input file
while IFS= read -r line; do
    input_file="$line"

    # Extract the sample name from the input file path
    sample_name=$(basename "${input_file}" .md.cram)

    # Run samtools view and save the output to a BAM file
    samtools view "${input_file}" HG992177.1:13411760-14019340 -o "${output_dir}/${sample_name}.sect.bam"
    tabix "${output_dir}/${sample_name}.sect.bam"
done < "$input_file_path"



AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTGGCCTGTGTAG

TTCAAGCAGAAGACGGCATACGAGATATGATTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

GGATCTTGTAAGCGCTTTAACAAAACTATGCAGCTTATTTTCTTTAATACACGTTTTACCAGATCGGAAGAGCACACGTCTGAACT-

TGCAGTATTTCTTGCTCTATCAAAAAACTGTGTAAAACAATCACAACCGGA-CACTGCTGGTAAACTTCATTATCTAATTATCTTTTATTTTTTTTGTCTTCT

AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAATCATATCTCGTATGCC--


AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAATCATATCTCGTATGCCGTCTTC


AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCATCCGGATCTCGTATGCCGTCTTCTGCTTGAAAAA-

AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCATCCGGATCTCGTATGCCGTCTTCTGCTTGAAAAT-

AGATCGGAAGAGCGTCGTGTAGGGAAAGA

TTAAATAGAAAGGTTACAGAATAATTCGAGCTGCAGCAATTATCACTGAATACAGCTTGCTGCTTTAGCGATGCC


/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/


### Mapping from before
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


fastp -i {input.R1} -I {input.R2} -p -c --merge --merged_out={output.merged} \
-o {output.R1_un} -O {output.R2_un}  -h {output.html} -j {output.json} -R '{params.report}' \
 -w {threads} -l {params.readlength} 2> {log}

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
     --adapter_fasta adapters.fasta \
     --dedup


--correction --detect_adapter_for_pe --cut_front 3 --cut_tail 3 \
 --qualified_quality-phred 20 --average_qual 20 --length_required 35 \
 --unqualified_percent_limit 10 --n_base_limit 0


 # global trimming options
 -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
 -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
 -b, --max_len1                       if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
 -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
 -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
 -B, --max_len2                       if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])

#Testing: #1
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
     --adapter_fasta adapters.fasta \
     --dedup

Filtering result:
reads passed filter: 1640
reads failed due to low quality: 34
reads failed due to too many N: 0
reads failed due to too short: 326
reads with adapter trimmed: 1916
bases trimmed due to adapters: 124066

nohup fastp \
    --in1 GAST_5_1943-L001_1.fastq.gz \
    --in2 GAST_5_1943-L001_2.fastq.gz \
    --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
    --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
    --json GAST_5_1943-L001.fastp.json \
    --html GAST_5_1943-L001.fastp.html \
    --thread 12 \
    --trim_front1=20 \
    --trim_tail1=20 \
    --detect_adapter_for_pe \
    -p \
    --adapter_fasta adapters.fasta \
    --dedup  \
    --correction --cut_front 3 --cut_tail 3 \
     --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
     --unqualified_percent_limit 10 --n_base_limit 0 &

     Filtering result:
     reads passed filter: 622
     reads failed due to low quality: 232
     reads failed due to too many N: 0
     reads failed due to too short: 1146
     reads with adapter trimmed: 1907
     bases trimmed due to adapters: 120193
     reads corrected by overlap analysis: 27
     bases corrected by overlap analysis: 56

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
         --adapter_fasta adapters.fasta \
         --dedup  \
         --correction --cut_front 5 --cut_tail 5 \
          --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
          --unqualified_percent_limit 10 --n_base_limit 0

#### Mapping test

# snakemake rules
rule map_historical:
    """Map trimmed and merged reads from historical samples to reference. BWA aln for short Illumina reads, parameters according to Palkopoulou et al. 2015"""
    input:
        ref=config["ref_path"],
        index=rules.bwa_index_reference.output,
        fastq_hist=rules.fastp_historical.output.merged,
    output:
        sai=temp("results/historical/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sai"),
    log:
        "results/logs/2_mapping/historical/" + REF_NAME + "/{sample}_{index}_{lane}_map_historical.log",
    threads: 8
    singularity:
        "docker://biocontainers/bwa:v0.7.17-3-deb_cv1"
    shell:
        """
        bwa aln -l 16500 -n 0.01 -o 2 -t  {input.ref} {input.fastq_hist} > {output.sai} 2> {log}
        """

bwa aln -l 16500 -n 0.01 -o 2 -t 8 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna  GAST_5_1943-L001_2.fastp.merged.fastq.gz > GAST_5_1943-L001.sai


fastp -i {input.R1} -I {input.R2} -p -c --merge --merged_out={output.merged} -o {output.R1_un} -O {output.R2_un} \
-h {output.html} -j {output.json} -R '{params.report}' -w {threads} -l {params.readlength} 2> {log}

nohup fastp \
    --in1 GAST_5_1943-L001_1.fastq.gz \
    --in2 GAST_5_1943-L001_2.fastq.gz \
    --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
    --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
    --json GAST_5_1943-L001.fastp.json \
    --html GAST_5_1943-L001.fastp.html \
    --thread 12 \
    --trim_front1=20 \
    --trim_tail1=20 \
    --detect_adapter_for_pe \
    -p \
    --adapter_fasta adapters.fasta \
    --dedup  \
    --correction --cut_front 5 --cut_tail 5 \
    --merge --merged_out=GAST_5_1943-L001_2.fastp.merged.fastq.gz \
     --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
     --unqualified_percent_limit 10 --n_base_limit 0 &

     Filtering result:
     reads passed filter: 33442576
     reads failed due to low quality: 8720074
     reads failed due to too many N: 6110
     reads failed due to too short: 60912010
     reads with adapter trimmed: 98238796
     bases trimmed due to adapters: 6229092140
     reads corrected by overlap analysis: 902705
     bases corrected by overlap analysis: 1856620

bwa samse /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna GAST_5_1943-L001.sai GAST_5_1943-L001_2.fastp.merged.fastq.gz | \
     samtools sort -@ {threads} - > GAST_5_1943-L001_2.fastp.merged.bam


     bwa samse /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna GAST_5_1943-L001.sai GAST_5_1943-L001_2.fastp.merged.fastq.gz | \
          samtools sort - > GAST_5_1943-L001_2.fastp.merged.bam

bwa samse ref.fa reads.sai reads.fq > aln-se.sam

#### Testing GenErode

conda env create -n generode -f environment.yml


rule fastp_historical:
    """Remove adapters, quality trim (phred 15) and merge overlapping paired-end reads in historical samples"""
    """fastp automatically detects adapter sequences for removal"""
    """NextSeq and NovaSeq samples are automatically detected and poly-G tails are removed"""
    """Minimum read length specified in config file"""
    input:
        R1=rules.fastq_historical_symbolic_links.output.fastq_r1,
        R2=rules.fastq_historical_symbolic_links.output.fastq_r2,
    output:
        R1_un=temp("results/historical/trimming/{sample}_{index}_{lane}_R1_unmerged.fastq.gz"),
        R2_un=temp("results/historical/trimming/{sample}_{index}_{lane}_R2_unmerged.fastq.gz"),
        merged="results/historical/trimming/{sample}_{index}_{lane}_trimmed_merged.fastq.gz",
        html="results/historical/trimming/stats/{sample}_{index}_{lane}_fastp_report.html",
        json=temp("results/modern/trimming/stats/{sample}_{index}_{lane}_fastp_report.json"),
    params:
        readlength=config["hist_readlength"],
        report="fastp report for {sample}_{index}_{lane}",
    log:
        "results/logs/1.1_fastq_processing/historical/{sample}_{index}_{lane}_fastp_historical.log",
    threads: 4
    singularity:
        "docker://quay.io/biocontainers/fastp:0.22.0--h2e03b76_0"
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -p -c --merge --merged_out={output.merged} -o {output.R1_un} -O {output.R2_un} \
        -h {output.html} -j {output.json} -R '{params.report}' -w {threads} -l {params.readlength} 2> {log}
        """




#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-00:00:00
#SBATCH -J remapping
#SBATCH --output=remapping.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/02_Benchmarking_trimming
bwa aln -l 16500 -n 0.01 -o 2 -t 16 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna  GAST_5_1943-L001_2.fastp.merged.fastq.gz > GAST_5_1943-L001.sai
bwa samse /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna GAST_5_1943-L001.sai GAST_5_1943-L001_2.fastp.merged.fastq.gz | samtools sort -@ 16 - > GAST_5_1943-L001_2.fastp.merged.bam


#### Trying GenErode

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode

module load Snakemake

########### Example ##########
#historical
samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file
S01_07_L2 BHJ5YVABCD.L2.07 illumina /proj/DeepSeq/P12345_1007_S7_L002_R1_001.fastq.gz /proj/DeepSeq/P12345_1007_S7_L002_R2_001.fastq.gz
#contemporary
samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file
S47_01_L7 AHCHL7XYZA.L7.01 illumina /proj/modern/P1234_101_S47_L007_R1_001.fastq.gz /proj/modern/P1234_101_S47_L007_R2_001.fastq.gz
########### Example ##########

@A00689:768:HKFY7DSX5:1:1101:4698:1047 2:N:0:CGACCTG+TACTCGC
@A00689:768:HKFY7DSX5:1:1101:14696:1016 2:N:0:AACCTGC+CGCAAGG


samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file
GAST51943_01_L1 AHKFY7DSX5.L001.S15 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R2_001.fastq.gz

samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file
KALM12018_01_L1 AHKFY7DSX5.L001.S49 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R2_001.fastq.gz

nano historical_samples_paths_D.txt
nano modern_samples_paths_D.txt

P27562_1001_S1_L001_R1_001_AHKFY7DSX5

GAST_5_1943,GAST_5_1943,L001,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R1_001.fastq.gz,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R2_001.fastq.gz
KALM_1_2018,KALM_1_2018,L001,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R1_001.fastq.gz,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R2_001.fastq.gz


KALM_4_1957 - diamina
ALTA_2_2015 - diamina
BAJK_4_2016 - diamina ?
ITAL_1_1946 - diamina
VORO_1_1998 - athalia
VAST_3_2005 - athalia
URAL_1_1991 - athalia
RUSS_2_1999 - athalia
SMAL_6_2013 - athalia
KRAS_4_2002 - athalia
KRAS_2_2002 - athalia
CHEH_3_2008 - athalia
CHEH_2_2004 - athalia
POLA_1_2003 - athalia

BAJK_3_2016
VAST_2_1999

#Quality red flags, deph lower then 2.5
ind depth
KALM_7_1961 2.25957
STOC_3_1965 0.960228
GAST_13_1969 1.40341
GAST_8_1965 1.10823
GAST_9_1965 1.31725
SMAL_1_1967 1.97408
STOC_1_1965 0.520048


#### 3) Edit the file "config/config.yaml"

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna

Fixing headers

sed -e 's/\(^>[^ ]*\) .*/\1/' GCA_905220545.2_ilMelAtha1.2_genomic.fa > GCA_905220545.2_ilMelAtha1.2_genomic.simple_headers.fa

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa


#Running Snakemake

module load fastp


snakemake -j 100 --use-singularity --cluster-config config/slurm/cluster.yaml --cluster "sbatch -A naiss2023-5-52 -p core --ntasks 17 --cpus-per-task {cluster.cpus-per-task} -t {cluster.time}" -npr &> YYMMDD_dry_run.out

snakemake -j 100 --use-singularity --cluster-config config/slurm/cluster.yaml --cluster "sbatch -A naiss2023-5-52 -p core --ntasks 17 -t 1-00:00:00" -npr &> 231206_1_dry_run.out

snakemake -j 100 --use-singularity --cluster-config config/slurm/cluster.yaml --cluster "sbatch -A naiss2023-5-52"

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config


screen -S britomartis


screen -r britomartis


######## Making sample sample lists

BAJK_1_2016,BAJK_1_2016,L001,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1039/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1039_S39_L001_R1_001.fastq.gz,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1039/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1039_S39_L001_R2_001.fastq.gz
GAST_6_1965,GAST_6_1965,L001,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1004/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1004_S4_L001_R1_001.fastq.gz,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1004/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1004_S4_L001_R2_001.fastq.gz
КALM_11_1997,КALM_11_1997,L001,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1036/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1036_S36_L001_R1_001.fastq.gz,/proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1036/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1036_S36_L001_R2_001.fastq.gz

KALM12018_01_L1 AHKFY7DSX5.L001.S49 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1049/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1049_S49_L001_R2_001.fastq.gz

awk -F "_"  'print $1,$2,$3,'_01_L1'' smaplelist.tsv

awk -F'[_|,]' -v OFS='' '{print $1, $2, $3, "_01_L1 AHKFY7DSX5.L001.", $18," illumina "}' smaplelist.tsv > sample_section.list
awk -F',' -v OFS=' ' '{print $4,$5}' smaplelist.tsv > sample_section2.list

paste -d '' sample_section.list sample_section2.list

CHEH22004_01_L1 AHKFY7DSX5.L001.S48 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1048/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1048_S48_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1048/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1048_S48_L001_R2_001.fastq.gz
RUSS52008_01_L1 AHKFY7DSX5.L001.S42 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1042/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1042_S42_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1042/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1042_S42_L001_R2_001.fastq.gz
VAST21999_01_L1 AHKFY7DSX5.L001.S66 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1066/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1066_S66_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1066/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1066_S66_L001_R2_001.fastq.gz
GAST51943_01_L1 AHKFY7DSX5.L001.S15 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1015/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1015_S15_L001_R2_001.fastq.gz
POLA22003_01_L1 AHKFY7DSX5.L001.S73 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1073/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1073_S73_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1073/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1073_S73_L001_R2_001.fastq.gz
KRAS12002_01_L1 AHKFY7DSX5.L001.S52 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1052/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1052_S52_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1052/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1052_S52_L001_R2_001.fastq.gz
BELA12004_01_L1 AHKFY7DSX5.L001.S45 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1045/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1045_S45_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1045/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1045_S45_L001_R2_001.fastq.gz
UPPS21958_01_L1 AHKFY7DSX5.L001.S23 illumina /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1023/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1023_S23_L001_R1_001.fastq.gz /proj/uppstore2017185/b2014034/private/raw_data/Assmann/DataDelivery_2023-03-07_15-07-05_ngisthlm00193/files/P27562/P27562_1023/02-FASTQ/230303_A00689_0768_AHKFY7DSX5/P27562_1023_S23_L001_R2_001.fastq.gz


snakemake -j 100 --use-singularity --cluster-config config/slurm/cluster.yaml --rerun-incomplete --cluster "sbatch -A naiss2023-5-52 -p core --ntasks 16 -t 1-00:00:00" -npr &> 231206_3_dry_run.out

snakemake -j 100 --use-singularity --cluster-config config/slurm/cluster.yaml --rerun-incomplete --cluster "sbatch -A naiss2023-5-52"


samtools view KALM12018_01_L1.sorted.bam  HG992177.1:13411760-14019340 -o KALM12018_01_L1.sorted.sect.cram


snakemake -j 100 --use-singularity --cluster-config config/slurm/cluster.yaml --rerun-incomplete --cluster "sbatch -A naiss2023-5-52" -np


# Modifying GenErode
# Changed trimming settings
### trying to run


# Back to ANGSD

less

# Historical

# contemporary


/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/data/raw_reads_symlinks/historical

fastp \
    --in1 GAST_5_1943-L001_1.fastq.gz \
    --in2 GAST_5_1943-L001_2.fastq.gz \
    --out1 GAST_5_1943-L001_1.fastp.fastq.gz \
    --out2 GAST_5_1943-L001_2.fastp.fastq.gz \
    --json GAST_5_1943-L001.fastp.json \
    --html GAST_5_1943-L001.fastp.html \
    --thread 12 \
    --trim_front1=20 \
    --trim_tail1=20 \
    --detect_adapter_for_pe \
    -p \
    --adapter_fasta adapters.fasta \
    --dedup  \
    --correction --cut_front 5 --cut_tail 5 \
    --merge --merged_out=GAST_5_1943-L001_2.fastp.merged.fastq.gz \
     --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
     --unqualified_percent_limit 10 --n_base_limit 0




 #!/bin/bash

 #SBATCH --job-name=fastp_processing
 #SBATCH --output=fastp_%A_%a.out
 #SBATCH --error=fastp_%A_%a.err
 #SBATCH --array=1-<number_of_lines>/2
 #SBATCH --time=02:00:00
 #SBATCH --cpus-per-task=12
 #SBATCH --mem=4G

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J hist_fastp
#SBATCH --output=historical_trimming.out
#SBATCH --array=1-33
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual

module load bioinfo-tools fastp

# Reading file names from an input file, where each file is listed on a new line
FILE_LIST="path.list"

# Calculate the line numbers for R1 and R2 files
LINE_NUM=$((SLURM_ARRAY_TASK_ID * 2 - 1))
NEXT_LINE_NUM=$((LINE_NUM + 1))

# Extract the file names
FILE_R1=$(sed -n "${LINE_NUM}p" $FILE_LIST)
FILE_R2=$(sed -n "${NEXT_LINE_NUM}p" $FILE_LIST)

# Extracting the common prefix from the first filename
PREFIX=$(echo $FILE_R1 | cut -d'_' -f 1,2,3,4,5,6)

# Run the fastp command
fastp \
   --in1 $FILE_R1 \
   --in2 $FILE_R2 \
   --out1 ${PREFIX}_1.fastp.fastq.gz \
   --out2 ${PREFIX}_2.fastp.fastq.gz \
   --json ${PREFIX}.fastp.json \
   --html ${PREFIX}.fastp.html \
   --thread 10 \
   --trim_front1=10 \
   --trim_tail1=10 \
   --detect_adapter_for_pe \
   -p \
   --dedup \
   --correction --cut_front 5 --cut_tail 5 \
   --merge --merged_out=${PREFIX}_2.fastp.merged.fastq.gz \
   --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
   --unqualified_percent_limit 10 --n_base_limit 0


awk 'print {'/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/data/raw_reads_symlinks/historical',$1}' historical_symlinks.list


/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/trimming

#Controllinf this trimming

ls -1 *merged* > file_list.txt

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:20:00
#SBATCH -J hist_fastqc
#SBATCH --output=historical_trimmingQC.out
#SBATCH --array=1-33
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/trimming

module load bioinfo-tools fastp FastQC

# Reading file names from an input file, where each file is listed on a new line
FILE_LIST="file_list.txt"

# Calculate the line numbers for R1 and R2 files
#LINE_NUM=$((SLURM_ARRAY_TASK_ID * 2 - 1))

# Extract the file names
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

# Extracting the common prefix from the first filename
PREFIX=$(echo $FILE_R1 | cut -d'_' -f 1,2,3,4,5,6)

# Run the fastp command
fastqc $FILE_1


#Mapping

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J remapping
#SBATCH --output=remapping.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-33

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/trimming

module load bioinfo-tools bwa samtools

# Reading file names from an input file, where each file is listed on a new line
FILE_LIST="file_list.txt"

# Extract the file names
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

# Extracting the common prefix from the first filename
PREFIX=$(echo $FILE_R1 | cut -d'_' -f 1)

#Run alignment
bwa aln -l 16500 -n 0.01 -o 2 -t 20 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa  $FILE_1 > $PREFIX.sai
bwa samse /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa $PREFIX.sai $FILE_1 | samtools sort -@ 20 - > $FILE_1.fastp.sorted.bam


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J remapping
#SBATCH --output=remapping.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-33

# Load required modules
module load bwa
module load samtools

# Define variables (replace these with actual paths or values)
REF_PATH="/path/to/reference"
INDEX="/path/to/index"
FASTQ_MOD_R1="/path/to/fastq_mod_R1"
FASTQ_MOD_R2="/path/to/fastq_mod_R2"
RG_PATH="/path/to/readgroup_file"
OUTPUT_BAM="/path/to/output.bam"
LOG_FILE="/path/to/logfile.log"

# Run the command
bwa mem -M -t 8 -R $(cat ${RG_PATH}) ${REF_PATH} ${FASTQ_MOD_R1} ${FASTQ_MOD_R2} | \
samtools sort -@ 8 -o ${OUTPUT_BAM} - 2> ${LOG_FILE}




#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J remapping
#SBATCH --output=remapping.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-33

# Stop on error
set -e

# Directories and file paths
PROJECT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"
TRIM_DIR="${PROJECT_DIR}/05_manual/historical/trimming"
REF_DIR="${PROJECT_DIR}/04_GenErode/GenErode/config/reference"
REF_GENOME="${REF_DIR}/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
FILE_LIST="${TRIM_DIR}/file_list.txt"

# Load modules
module load bioinfo-tools bwa samtools

# Navigate to the trimming directory
cd "$TRIM_DIR"

# Extract the file names
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")

# Extracting the common prefix from the first filename
PREFIX=$(echo "$FILE_1" | cut -d'_' -f 1)

# Run alignment
bwa aln -l 16500 -n 0.01 -o 2 -t 20 "$REF_GENOME" "$FILE_1" > "${PREFIX}.sai"
bwa samse "$REF_GENOME" "${PREFIX}.sai" "${PREFIX}" | samtools sort -@ 20 -o "${FILE_1}.fastp.sorted.bam"

# Clean up intermediate files
rm "${PREFIX}.sai"

# Add optional commands for error output, if any specific error logging is required


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J remapping
#SBATCH --output=remapping.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-33

echo "Starting Dry Run of SLURM Job"

# Directories and file paths
echo "Setting up directories and file paths..."
PROJECT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"
TRIM_DIR="${PROJECT_DIR}/05_manual/historical/trimming"
REF_DIR="${PROJECT_DIR}/04_GenErode/GenErode/config/reference"
REF_GENOME="${REF_DIR}/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
FILE_LIST="${TRIM_DIR}/file_list.txt"

# Simulating module load
echo "Loading modules: bioinfo-tools bwa samtools"

# Simulating navigation to the trimming directory
echo "Navigating to directory: $TRIM_DIR"

# Extract the file names (simulated)
echo "Extracting file name from list..."
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")
echo "FILE_1: $FILE_1"

# Extracting the common prefix from the first filename
PREFIX=$(echo "$FILE_1" | cut -d'_' -f 1)
echo "PREFIX: $PREFIX"

# Display the bwa alignment command
echo "bwa aln -l 16500 -n 0.01 -o 2 -t 20 $REF_GENOME $FILE_1 > ${PREFIX}.sai"

# Display the bwa samse and samtools sort commands
echo "bwa samse $REF_GENOME ${PREFIX}.sai $FILE_1 | samtools sort -@ 20 -o ${PREFIX}.sorted.bam"

# Simulating cleanup of intermediate files
echo "Cleaning up intermediate files: rm ${PREFIX}.sai"

echo "Dry Run Completed"




#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J remapping
#SBATCH --output=remapping.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-33

# Stop on error
set -e

# Directories and file paths
PROJECT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"
TRIM_DIR="${PROJECT_DIR}/05_manual/historical/trimming"
REF_DIR="${PROJECT_DIR}/04_GenErode/GenErode/config/reference"
REF_GENOME="${REF_DIR}/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
FILE_LIST="${TRIM_DIR}/file_list.txt"

# Load modules
module load bioinfo-tools bwa samtools

# Navigate to the trimming directory
cd "$TRIM_DIR"

# Extract the file names
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")

# Extracting the common prefix from the first filename
PREFIX=$(echo "$FILE_1" | cut -d'_' -f 1)

# Run alignment
bwa aln -l 16500 -n 0.01 -o 2 -t 20 "$REF_GENOME" "$FILE_1" > "${PREFIX}.sai"
bwa samse "$REF_GENOME" "${PREFIX}.sai" "$FILE_1" | samtools sort -@ 20 -o "${PREFIX}.fastp.sorted.bam"

# Clean up intermediate files
rm "${PREFIX}.sai"

# Add optional commands for error output, if any specific error logging is required


#Remapping modern

# 1. Remove really bad samples

BELA_1_2004
KALM_10_1983
KALM111997
SMAL_5_19
RUSS_6_2008
JAPA_2_1994

### 2. Trim Agin

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 01:00:00
#SBATCH -J hist_fastp
#SBATCH --output=modern_super_trimming.out
#SBATCH --array=1-33
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern

module load bioinfo-tools fastp FastQC

# Reading file names from an input file, where each file is listed on a new line
FILE_LIST="path.list"

# Calculate the line numbers for R1 and R2 files
LINE_NUM=$((SLURM_ARRAY_TASK_ID * 2 - 1))
NEXT_LINE_NUM=$((LINE_NUM + 1))

# Extract the file names
FILE_R1=$(sed -n "${LINE_NUM}p" $FILE_LIST)
FILE_R2=$(sed -n "${NEXT_LINE_NUM}p" $FILE_LIST)

# Extracting the common prefix from the first filename
PREFIX=$(echo $FILE_R1 | cut -d'_' -f 1)

# Run the fastp command
fastp \
   --in1 $FILE_R1 \
   --in2 $FILE_R2 \
   --out1 ${PREFIX}_1.superfastp.fastq.gz \
   --out2 ${PREFIX}_2.superfastp.fastq.gz \
   --json ${PREFIX}.superfastp.json \
   --html ${PREFIX}.superfastp.html \
   --thread 10 \
   --detect_adapter_for_pe \
   -p \
   --dedup \
   --correction --cut_front 5 --cut_tail 5 \
   --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
   --unqualified_percent_limit 10 --n_base_limit 0

# Run the fastqc command
fastqc ${PREFIX}_1.superfastp.fastq.gz
fastqc ${PREFIX}_2.superfastp.fastq.gz





###ls -1 *merged* > file_list.txt

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 01:00:00
#SBATCH -J modern_super
#SBATCH --output=modern_super_trimming.out
#SBATCH --array=1-21
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Base directory
BASE_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern"

# Navigate to the base directory
cd "$BASE_DIR"

# Load required modules
module load bioinfo-tools fastp FastQC

# File list
FILE_LIST="$BASE_DIR/path.list"

# Ensure the file list exists
if [ ! -f "$FILE_LIST" ]; then
    echo "File list not found: $FILE_LIST"
    exit 1
fi

# Calculate line numbers for R1 and R2 files
LINE_NUM=$((SLURM_ARRAY_TASK_ID * 2 - 1))
NEXT_LINE_NUM=$((LINE_NUM + 1))

# Extract file names
FILE_R1=$(sed -n "${LINE_NUM}p" "$FILE_LIST")
FILE_R2=$(sed -n "${NEXT_LINE_NUM}p" "$FILE_LIST")

# Ensure the input files exist
if [ ! -f "$FILE_R1" ] || [ ! -f "$FILE_R2" ]; then
    echo "Input files not found: $FILE_R1, $FILE_R2"
    exit 1
fi

# Extract common prefix from the first filename
PREFIX=$(echo "$FILE_R1" | cut -d'_' -f 1)

# Run the fastp command
fastp \
   --in1 "$FILE_R1" \
   --in2 "$FILE_R2" \
   --out1 "${PREFIX}_1.superfastp.fastq.gz" \
   --out2 "${PREFIX}_2.superfastp.fastq.gz" \
   --json "${PREFIX}.superfastp.json" \
   --html "${PREFIX}.superfastp.html" \
   --thread 10 \
   --detect_adapter_for_pe \
   -p \
   --dedup \
   --correction --cut_front 5 --cut_tail 5 \
   --qualified_quality_phred 20 --average_qual 30 --length_required 35 \
   --unqualified_percent_limit 10 --n_base_limit 0

# Run the fastqc command
fastqc "${PREFIX}_1.superfastp.fastq.gz"
fastqc "${PREFIX}_2.superfastp.fastq.gz"


# /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern







#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J remappingm
#SBATCH --output=remappingm.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-21

# Stop on error
set -e

# Directories and file paths
PROJECT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"
TRIM_DIR="${PROJECT_DIR}/05_manual/historical/trimming"
REF_DIR="${PROJECT_DIR}/04_GenErode/GenErode/config/reference"
REF_GENOME="${REF_DIR}/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
FILE_LIST="${TRIM_DIR}/file_list.txt"
OUTPUT_DIR="${PROJECT_DIR}/05_manual/modern/mapping"

# Load modules
module load bioinfo-tools bwa samtools picard

# Navigate to the trimming directory
cd "$TRIM_DIR"

# Calculate line numbers for R1 and R2 files
LINE_NUM=$((SLURM_ARRAY_TASK_ID * 2 - 1))
NEXT_LINE_NUM=$((LINE_NUM + 1))

# Extract file names
FILE_R1=$(sed -n "${LINE_NUM}p" "$FILE_LIST")
FILE_R2=$(sed -n "${NEXT_LINE_NUM}p" "$FILE_LIST")

# Extracting the common prefix from the first filename
PREFIX=$(echo "$FILE_R1" | cut -d'_' -f 1)

# Run alignment
bwa mem -M -t 20 "$REF_GENOME" "$FILE_R1" "$FILE_R2" | \
samtools sort -@ 20 - > "${OUTPUT_DIR}/${PREFIX}.fastp.sorted.bam"
java -jar $PICARD_ROOT/picard.jar MarkDuplicates -Xmx8g \
   INPUT="${OUTPUT_DIR}/${PREFIX}.fastp.sorted.bam" \
   OUTPUT="${OUTPUT_DIR}/${PREFIX}.fastp.sorted.rmdup.bam" \
   METRICS_FILE="${OUTPUT_DIR}/${PREFIX}.merged.rmdup_metrics.txt"









########
## Remove duplicates
########

# historical

rule rmdup_historical_bams:
    """Remove PCR duplicates from historical samples using Pontus Skoglund's custom script for duplicate removal (checks both ends of a read)"""
    """The script was modified so that also unmapped reads are printed to the output bam file so that they are not lost"""
    input:
        merged=rules.merge_historical_bams_per_index.output.merged,
        index="results/historical/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai",
    output:
        rmdup=temp("results/historical/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam"),
    threads: 6
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/{sample}_{index}_rmdup_historical_bams.log",
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        """
        samtools view -@ {threads} -h {input.merged} | python3 workflow/scripts/samremovedup.py  | samtools view -b -o {output.rmdup} 2> {log}
        """


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 01:00:00
#SBATCH -J hist_rmdup
#SBATCH --output=hist_rmdup.out
#SBATCH --array=1-33
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load required modules
module load bioinfo-tools samtools

# Base directories
TRIM_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/trimming"
OUTPUT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/mapping"
cd "${TRIM_DIR}"

# Path to the Python script for duplicate removal
PYTHON_SCRIPT="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/workflow/scripts/samremovedup.py"

# File list
FILE_LIST="${TRIM_DIR}/bam_list.txt"

# Extract the file names
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")

# Extracting the common prefix from the filename
PREFIX=$(basename "$FILE_1" | cut -d'.' -f 1)

# Define output BAM file
OUTPUT_BAM="${OUTPUT_DIR}/${PREFIX}.fastp.sorted.rmdup.bam"

# Run the command
samtools index "$FILE_1"
samtools view -@ 6 -h "$FILE_1" | python3 "$PYTHON_SCRIPT" | samtools view -b -o "$OUTPUT_BAM"
samtools index "$OUTPUT_BAM"


ls -1 *sorted.bam* > bam_list.txt

# modern
rule rmdup_modern_bams:
    """Mark duplicates in modern samples"""
    input:
        merged=rules.merge_modern_bams_per_index.output.merged,
        index="results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai",
    output:
        rmdup=temp("results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam"),
        metrix=temp("results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup_metrics.txt"),
    threads: 2
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/{sample}_{index}_rmdup_modern_bams.log",
    singularity:
        "docker://quay.io/biocontainers/picard:2.26.6--hdfd78af_0"
    shell:
        ""
        mem=$(((6 * {threads}) - 2))
        picard MarkDuplicates -Xmx${{mem}}g INPUT={input.merged} OUTPUT={output.rmdup} METRICS_FILE={output.metrix} 2> {log}
        ""

        #!/bin/bash
        #SBATCH -A naiss2023-5-52
        #SBATCH -p core
        #SBATCH -n 6
        #SBATCH -t 01:00:00
        #SBATCH -J hist_rmdup
        #SBATCH --output=hist_rmdup.out
        #SBATCH --array=1-33
        #SBATCH --mail-user=daria.shipilina@gmail.com
        #SBATCH --mail-type=ALL

        echo "Starting Dry Run"

        # Simulating module load
        echo "Loading modules: bioinfo-tools, samtools"

        # Base directories
        TRIM_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/trimming"
        OUTPUT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/mapping"
        echo "Changing directory to TRIM_DIR: ${TRIM_DIR}"

        # Path to the Python script for duplicate removal
        PYTHON_SCRIPT="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/workflow/scripts/samremovedup.py"

        # File list
        FILE_LIST="${TRIM_DIR}/bam_list.txt"
        echo "Using file list: ${FILE_LIST}"

        # Extract the file names
        echo "Extracting file name based on SLURM_ARRAY_TASK_ID"
        FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")
        echo "FILE_1: $FILE_1"

        # Extracting the common prefix from the filename
        echo "Extracting prefix from FILE_1"
        PREFIX=$(basename "$FILE_1" | cut -d'.' -f 1)
        echo "PREFIX: $PREFIX"

        # Define output BAM file
        OUTPUT_BAM="${OUTPUT_DIR}/${PREFIX}.fastp.sorted.rmdup.bam"
        echo "Output BAM file: $OUTPUT_BAM"

        # Display the samtools index command
        echo "samtools index command: samtools index \"$FILE_1\""

        # Display the samtools view and python script pipeline
        echo "samtools view and Python script pipeline:"
        echo "samtools view -@ 6 -h \"$FILE_1\" | python3 \"$PYTHON_SCRIPT\" | samtools view -b -o \"$OUTPUT_BAM\""

        # Display the final samtools index command
        echo "Final samtools index command: samtools index \"$OUTPUT_BAM\""

        echo "Dry Run Completed"




rmdup_historical_bams:
time: 3-00:00:00
cpus-per-task: 6
rmdup_modern_bams:
time: 3-00:00:00
cpus-per-task: 2


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 01:00:00
#SBATCH -J hist_rmdup
#SBATCH --output=hist_rmdup.out
#SBATCH --array=1-21
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load required modules
module load bioinfo-tools samtools picard

# Base directories
TRIM_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern/trimming"
OUTPUT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern/mapping"
cd "${TRIM_DIR}"

# File list
FILE_LIST="${TRIM_DIR}/bam_list.txt"

# Extract the file names
FILE_1=$(sed -n 1p "$FILE_LIST")

# Extracting the common prefix from the filename
PREFIX=$(basename "$FILE_1" | cut -d'.' -f 1)

# Calculate memory allocation for Picard
THREADS=2
MEM=$((6 * THREADS - 2))

# Run Picard MarkDuplicates command
java -jar picard.jar MarkDuplicates -Xmx${MEM}g \
   INPUT="$FILE_1" \
   OUTPUT="${FILE_1}.rmdup.bam" \
   METRICS_FILE="${METRICS_FILE}.merged.rmdup_metrics.txt"



   # Base directories
   TRIM_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/trimming"
   OUTPUT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/mapping"
   cd "${TRIM_DIR}"

   # Path to the Python script for duplicate removal
   PYTHON_SCRIPT="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/workflow/scripts/samremovedup.py"

   # File list
   FILE_LIST="${TRIM_DIR}/bam_list.txt"


########
## BAM manual QC
########
module load samtools
samtools index DALA11965.fastp.sorted.bam
samtools view DALA11965.fastp.sorted.bam  HG992177.1:13411760-14019340 -o bam_visual_control/DALA11965.fastp.sorted.sect.bam

samtools index GAST51943.fastp.sorted.bam
samtools view GAST51943.fastp.sorted.bam  HG992177.1:13411760-14019340 -o bam_visual_control/GAST51943.fastp.sect.bam

samtools index SMAL21967.fastp.sorted.bam
samtools view SMAL21967.fastp.sorted.bam  HG992177.1:13411760-14019340 -o bam_visual_control/SMAL21967.fastp.sorted.sect.bam

########
## ANGSD
########


echo "${PREFIX}"


# Modifying GenErode
# Changed trimming settings
### trying to run


# Back to ANGSD

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools ANGSD
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/

angsd -b rerun_Dec7/filtered_paths_historical_sweden.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna \
-anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna \
-out rerun_Dec7/historical_sweden -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 \
 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -skipTriallelic 1 \
 -SNP_pval 1e-3 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf coordinates_noZ_nomt.list

# Historical

# contemporary
#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J angsd
#SBATCH --output=angsd_sweden.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools ANGSD
cd /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/

angsd -b rerun_Dec7/filtered_paths_modern_world.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna \
-anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/reference_genome_M.athalia/GCA_905220545.2_ilMelAtha1.2_genomic.fna \
-out rerun_Dec7/modern_world -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 \
 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -skipTriallelic 1 \
 -SNP_pval 1e-3 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 10 -rf coordinates_noZ_nomt.list



 [detached from 31550.pts-8.rackham1]


 ###################
 ##### Retrimmed ANGSD
##################

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J angsd_hist_75
#SBATCH --output=angsd_hist_75.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Base directory
BASE_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"

# Minimum number of individuals
MIN_IND=25

# Load modules
module load bioinfo-tools ANGSD

# Change to ANGSD directory
ANGSD_DIR="${BASE_DIR}/05_manual/historical/calling_angsd"
cd "$ANGSD_DIR"

# Define paths
REF_GENOME="${BASE_DIR}/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
FILTERED_PATHS="${BASE_DIR}/05_manual/historical/calling_angsd/bam_paths_historical_sweden.txt"
COORDINATES_LIST="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list"

# Run ANGSD
angsd -b "${BASE_DIR}/05_manual/historical/calling_angsd/bam_paths_historical_sweden.txt" -ref "$REF_GENOME" \
-anc "$REF_GENOME" \
-out "${ANGSD_DIR}/ansgd_hist_${MIN_IND}" -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2  -setMinDepth 5 -setMaxDepth 30 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -skipTriallelic 1 \
-SNP_pval 1e-3 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd "$MIN_IND" -rf "$COORDINATES_LIST"




ls -1 *sorted.rmdup.bam* > bam_list.txt

awk -F"/" '{print "/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/mapping/"$3}' bam_paths_historical_sweden.txt
awk -F"/" '{print $3}' bam_paths_historical_sweden.txt

awk -F"/" '{print "/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/mapping/"$11}' bam_paths_historical_sweden.txt


angsd -b /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd/bam_paths_historical_sweden_.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd/ansgd_hist_25 -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2 -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -skipTriallelic 1 -SNP_pval 1e-3 -doMaf 1 -doCounts 1 -GL 2 -doSaf 1 -minmaf 0.05 -minInd 25 -rf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list



########
## PCA
########

/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/rerun_Dec7


pcangsd -beagle ansgd_hist_16.beagle.gz -out ansgd_hist_16 -threads 64


########
## vcf approach
########


angsd -b /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd/bam_paths_historical_sweden_.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -out /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd/ansgd_hist_16 -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2 -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -gl 1 -dopost 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -skipTriallelic 1 -SNP_pval 1e-6 -doMaf 1 -doCounts 1 -doSaf 1 -dobcf 1 -rf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list

bcftools stats ansgd_hist_16.bcf | head -n 30 | tail -n 9

(base) [daria@r123 calling_angsd]$ bcftools stats ansgd_hist_16.bcf | head -n 30 | tail -n 9
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	32
SN	0	number of records:	1567755
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	1567755
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0


bcftools view -Oz v -o ansgd_hist_16.vcf ../ansgd_hist_16.bcf


vcftools --bcf ../ansgd_hist_16.bcf --max-missing-count 3 --recode --stdout | gzip -c > ansgd_hist_16.qfilter.vcf.gz


angsd -b bam_paths_modern.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -anc /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -out ansgd_modern -minMapQ 20 -minQ 20 -nThreads 8 -doGlf 2 -setMinDepth 5 -setMaxDepth 30 -uniqueOnly 1 -remove_bads 1 -gl 1 -dopost 1 -only_proper_pairs 1 -trim 0 -C 50 -doMajorMinor 1 -skipTriallelic 1 -SNP_pval 1e-6 -doMaf 1 -doCounts 1 -doSaf 1 -dobcf 1 -rf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list

plink --vcf ansgd_hist_16.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --out britomartis.historical

plink --tped outnames.tped  --double-id --allow-extra-chr --set-missing-var-ids @:# --out britomartis.historical


angsd -bam /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd/bam_paths_historical_sweden_.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -out outnames -doPlink 2 -doGeno -4 -doMajorMinor 1 -doCounts 1 -doMaf 2 -postCutoff 0.99 -uniqueOnly 1 -remove_bads 1 -gl 1 -dopost 1 -only_proper_pairs 1 -trim 0 -C 50 -skipTriallelic 1 -SNP_pval 1e-6 -geno_minDepth 4 -rf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list


plink --tped outnames.tped --tfam outnames.tfam  --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out britomartis.historical


plink --tped outnames.tped --tfam outnames.tfam --mind 0.05 --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.historical.prune.in --make-bed --pca --out britomartis.historical


angsd -bam /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd/bam_paths_historical_sweden_.txt -ref /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa -out outnames -doPlink 2 -doGeno -4 -doMajorMinor 1 -doCounts 1 -doMaf 2 -postCutoff 0.99 -uniqueOnly 1 -remove_bads 1 -gl 1 -dopost 1 -only_proper_pairs 1 -trim 0 -C 50 -skipTriallelic 1 -SNP_pval 1e-6 -geno_minDepth 4 -rf /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list


#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH -J angsd_hist_plink
#SBATCH --output=angsd_hist_plink.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load necessary modules
module load bioinfo-tools ANGSD

# Define directories and file paths for easier management
PROJECT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"
BAM_PATHS="${PROJECT_DIR}/05_manual/historical/calling_angsd/bam_paths_historical_sweden_.txt"
REF_GENOME="${PROJECT_DIR}/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
COORDINATES_LIST="${PROJECT_DIR}/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list"
OUTPUT_DIR="${PROJECT_DIR}/05_manual/historical/calling_angsd" # Make sure this directory exists or create it before running the script

# Navigate to output directory
cd "$OUTPUT_DIR"

# Run ANGSD
angsd -bam "$BAM_PATHS" \
      -ref "$REF_GENOME" \
      -out angsd_hist_plink \
      -P 20 \
      -doPlink 2 \
      -doGeno -4 \
      -doMajorMinor 1 \
      -doCounts 1 \
      -doMaf 2 \
      -postCutoff 0.99 \
      -uniqueOnly 1 \
      -remove_bads 1 \
      -gl 1 \
      -dopost 1 \
      -only_proper_pairs 1 \
      -trim 0 \
      -C 50 \
      -skipTriallelic 1 \
      -SNP_pval 1e-6 \
      -geno_minDepth 4 \
      -rf "$COORDINATES_LIST"



#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 01:00:00
#SBATCH -J modern_rmdup
#SBATCH --output=modern_rmdup.out
#SBATCH --array=1-21
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load required modules
module load bioinfo-tools samtools picard

# Base directories
OUTPUT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern/mapping"
cd "${OUTPUT_DIR}"

# File list
FILE_LIST="${OUTPUT_DIR}/bam_paths_modern_minus1.txt"

# Extract the file names
FILE_1=$(sed -n "${SLURM_ARRAY_TASK_ID}p"  "$FILE_LIST")

# Extracting the common prefix from the filename
PREFIX=$(basename "$FILE_1" | cut -d'.' -f 1)

# Run Picard MarkDuplicates command
java -jar $PICARD_ROOT/picard.jar MarkDuplicates \
   INPUT="$PREFIX".fastp.sorted.bam \
   OUTPUT="$PREFIX".fastp.sorted.rmdup.bam \
   METRICS_FILE="$PREFIX".merged.rmdup_metrics.txt

samtools index "$PREFIX".fastp.sorted.rmdup.bam


 ls -1 *sorted.rmdup.bam$* | awk '{print "/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/modern/mapping/"$11}'

#!/bin/bash
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH -J angsd_modern_plink
#SBATCH --output=angsd_modern_plink.out
#SBATCH --mail-user=daria.shipilina@gmail.com
#SBATCH --mail-type=ALL

# Load necessary modules
module load bioinfo-tools ANGSD

# Define directories and file paths for easier management
PROJECT_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation"
BAM_PATHS="${PROJECT_DIR}/05_manual/modern/calling_angsd/bam_paths_modern.txt"
REF_GENOME="${PROJECT_DIR}/04_GenErode/GenErode/config/reference/GCA_905220545.2_ilMelAtha1.2_genomic.chroms.simple_headers.fa"
COORDINATES_LIST="${PROJECT_DIR}/00_Mapping_Calling_sarek/08_ANGSD/coordinates_noZ_nomt.list"
OUTPUT_DIR="${PROJECT_DIR}/05_manual/modern/calling_angsd" # Make sure this directory exists or create it before running the script

# Navigate to output directory
cd "$OUTPUT_DIR"

# Run ANGSD
angsd -bam "$BAM_PATHS" \
      -ref "$REF_GENOME" \
      -out angsd_hist_plink \
      -P 20 \
      -doPlink 2 \
      -doGeno -4 \
      -doMajorMinor 1 \
      -doCounts 1 \
      -doMaf 2 \
      -postCutoff 0.99 \
      -uniqueOnly 1 \
      -remove_bads 1 \
      -gl 1 \
      -dopost 1 \
      -only_proper_pairs 1 \
      -trim 0 \
      -C 50 \
      -skipTriallelic 1 \
      -SNP_pval 1e-6 \
      -geno_minDepth 4 \
      -rf "$COORDINATES_LIST"


-> Done reading data waiting for calculations to finish
       -> Done waiting for threads
       -> Output filenames:
               ->"angsd_hist_plink.arg"
               ->"angsd_hist_plink.mafs.gz"
               ->"angsd_hist_plink.tfam"
               ->"angsd_hist_plink.tped"
       -> Tue Dec 12 19:39:41 2023
       -> Arguments and parameters for all analysis are located in .arg file
       -> Total number of sites analyzed: 351317124
       -> Number of sites retained after filtering: 653249
       [ALL done] cpu-time used =  49440.23 sec
       [ALL done] walltime used =  6904.00 sec

########
## plink development
########

##### 1. Prune
plink --tped outnames.tped --tfam outnames.tfam --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 30 0.1 --out britomartis.hist

# 2. Check missing genotypes
plink --tped outnames.tped --tfam outnames.tfam --double-id --allow-extra-chr --set-missing-var-ids @:# --missing --out miss_stat

plink --tped outnames.tped --tfam outnames.tfam --geno 0.3 --maf 0.05 --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.hist.prune.in --make-bed --pca --out britomartis.PCAhist

# 3.



##### plink analysis historical

BASE_DIR="/crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/05_manual/historical/calling_angsd"
TPED_FILE="${BASE_DIR}/angsd_hist_plink.tped"
TFAM_FILE="${BASE_DIR}/angsd_hist_plink.tfam"

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --missing --out miss_stat_historical

Total genotyping rate is 0.075041

FID  IID MISS_PHENO   N_MISS   N_GENO   F_MISS
  1    1          Y   641996   653249   0.9828
  2    1          Y   650713   653249   0.9961
  3    1          Y   630056   653249   0.9645
  4    1          Y   648347   653249   0.9925
  5    1          Y   456373   653249   0.6986
  6    1          Y   651870   653249   0.9979
  7    1          Y   552817   653249   0.8463
  8    1          Y   650098   653249   0.9952
  9    1          Y   578185   653249   0.8851
 10    1          Y   602971   653249    0.923
 11    1          Y   608874   653249   0.9321
 12    1          Y   629880   653249   0.9642
 13    1          Y   650265   653249   0.9954
 14    1          Y   652207   653249   0.9984
 15    1          Y   651406   653249   0.9972
 16    1          Y   487451   653249   0.7462
 17    1          Y   617299   653249    0.945
 18    1          Y   618878   653249   0.9474
 19    1          Y   555393   653249   0.8502
 20    1          Y   642271   653249   0.9832
 21    1          Y   626884   653249   0.9596
 22    1          Y   405832   653249   0.6213
 23    1          Y   642805   653249    0.984
 24    1          Y   479258   653249   0.7337
 25    1          Y   626700   653249   0.9594
 26    1          Y   652802   653249   0.9993
 27    1          Y   615543   653249   0.9423
 28    1          Y   649904   653249   0.9949
 29    1          Y   582761   653249   0.8921
 30    1          Y   651365   653249   0.9971
 31    1          Y   605262   653249   0.9265
 32    1          Y   618848   653249   0.9473

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.3 --out britomartis.hist

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.4 --maf 0.05 --double-id \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.hist.prune.in \
    --make-bed --out britomartis.hist

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.4 --maf 0.05 --double-id --pca \
    --allow-extra-chr --mind 0.9800 --set-missing-var-ids @:# --extract britomartis.hist.prune.in  \
    --make-bed --out britomartisSTRING.hist

1013 variants and 22 people pass filters and QC.


plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.hist.prune.in  \
    --make-bed --out britomartisALLIND.hist

132 variants and 32 people pass filters and QC.


#!/bin/bash

# Path to the input file
FILE_PATH="bam_paths_historical_sweden.txt"
# Read each line in the file
while read -r line; do
    # Extract the filename without path and extension
    filename=$(basename "$line" | cut -d'.' -f 1)

    # Extract the first four characters and the last four characters
    first_four=$(echo "$filename" | cut -c 1-4)
    last_four=$(echo "$filename" | rev | cut -c 1-4 | rev)

    # Get the middle part of the filename
    middle_part=$(echo "$filename" | cut -c 5- | rev | cut -c 5- | rev)

    # Output the formatted string
    echo "${first_four}_${middle_part}_${last_four}"
done < "$FILE_PATH"



##### plink analysis modern

BASE_DIR="modern/calling_angsd"
TPED_FILE="${BASE_DIR}/angsd_modern_plink.tped"
TFAM_FILE="${BASE_DIR}/angsd_modern_plink.tfam"

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --missing --out miss_stat_modern

10253654 variants loaded from .bim file.
Total genotyping rate is 0.541703.

FID  IID MISS_PHENO   N_MISS   N_GENO   F_MISS
  1    1          Y  4774887 10253654   0.4657
  2    1          Y  7866506 10253654   0.7672
  3    1          Y  7156705 10253654    0.698
  4    1          Y  5878027 10253654   0.5733
  5    1          Y  2982585 10253654   0.2909
  6    1          Y  4112024 10253654    0.401
  7    1          Y  4851065 10253654   0.4731
  8    1          Y  3760292 10253654   0.3667
  9    1          Y  5683210 10253654   0.5543
 10    1          Y  4780817 10253654   0.4663
 11    1          Y  4261909 10253654   0.4156
 12    1          Y  3708229 10253654   0.3616
 13    1          Y  4375313 10253654   0.4267
 14    1          Y  4208238 10253654   0.4104
 15    1          Y  3500252 10253654   0.3414
 16    1          Y  5051678 10253654   0.4927
 17    1          Y  5456788 10253654   0.5322
 18    1          Y  3361741 10253654   0.3279
 19    1          Y  3551320 10253654   0.3463
 20    1          Y  5231210 10253654   0.5102
 21    1          Y  4130763 10253654   0.4029

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.3 --out britomartis.hist

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --double-id \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.hist.prune.in \
    --make-bed --out britomartisnoALT.modern

92080 variants and 21 people pass filters and QC.

#Extract ALTA2 sample

--remove remove.txt


plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.hist.prune.in  \
    --make-bed --out britomartisMISS01.modern

132 variants and 32 people pass filters and QC.

92080 variants and 21 people pass filters and QC.

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.hist.prune.in  \
    --make-bed --remove remove.txt --out britomartisNOALT2.modern


#!/bin/bash

# Path to the input file
FILE_PATH="bam_paths_modern.txt"
# Read each line in the file
while read -r line; do
    # Extract the filename without path and extension
    filename=$(basename "$line" | cut -d'.' -f 1)

    # Extract the first four characters and the last four characters
    first_four=$(echo "$filename" | cut -c 1-4)
    last_four=$(echo "$filename" | rev | cut -c 1-4 | rev)

    # Get the middle part of the filename
    middle_part=$(echo "$filename" | cut -c 5- | rev | cut -c 5- | rev)

    # Output the formatted string
    echo "${first_four}_${middle_part}_${last_four}"
done < "$FILE_PATH"


plink --tped angsd_modern_plink.tped --tfam angsd_modern_plink.tfam --geno 0.1 --maf 0.05 --remove remove.txt  --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.3 --out britomartis.modernall

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.modernall.prune.in  \
    --make-bed --remove remove.txt --out britomartisNOALT2.modern

plink --tped angsd_modern_plink.tped --tfam angsd_modern_plink.tfam  --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.modernall.prune.in  \
    --make-bed --remove remove.txt --recode vcf --out britomartisNOALT2.modern

    plink --tped angsd_modern_plink.tped --tfam angsd_modern_plink.tfam  --geno 0.1 --double-id --pca \
        --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.modernall.prune.in  \
        --make-bed --remove remove.txt --out

# Calculate nucleotide diversity (Pi)
vcftools --vcf britomartisNOALT2.modern.vcf --site-pi --out piNOALT2_output
#vcftools --vcf britomartisNOALT2.modern.vcf --window-pi 50000 --out piNOALT2.modernworld
vcftools --vcf britomartisNOALT2.modern.vcf --TajimaD 1000 --out tajNOALT2.modernworld
vcftools --vcf britomartisNOALT2.modern.vcf --relatedness --out relateNOALT2.modernworld

awk '{if ($3 != 'nan') {sum += $3; count++}} END {if (count > 0) print sum/count; else print "No non-zero values"}' tajNOALT2.modernworld.Tajima.D

Tajimas D: 1.891


# Calculate the mean of the 5th column
mean_pi=$(awk '{sum += $3} END {print sum/NR}' piNOALT2_output.sites.pi)

# Output the result
echo "mean pi: $mean_pi"

mean pi: 0.237203







plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --recode vcf --double-id \
    --allow-extra-chr --set-missing-var-ids @:# \
    --make-bed --remove remove_sweden.txt --out britomartisSWEDEN.modern




plink --bfile britomartisSWEDEN.modern --double-id \
    --allow-extra-chr --set-missing-var-ids @:# --genome --out out


   FID1 IID1 FID2 IID2 RT    EZ      Z0      Z1 Z2  PI_HAT PHE       DST     PPC   RATIO
    9   1  13   1 UN    NA  0.9709  0.0291  0.0000  0.0146  -1  0.834913  0.8548  2.1568
     4   1  15   1 UN    NA  0.9737  0.0091  0.0171  0.0217  -1  0.839148  0.9966  2.4328
     3   1  15   1 UN    NA  0.9522  0.0478  0.0000  0.0239  -1  0.833183  0.9846  2.3382
     3   1   5   1 UN    NA  0.9374  0.0626  0.0000  0.0313  -1  0.837214  0.8075  2.1276
    10   1  13   1 UN    NA  0.9321  0.0603  0.0076  0.0378  -1  0.841313  1.0000  2.8529
     1   1   9   1 UN    NA  0.9592  0.0000  0.0408  0.0408  -1  0.841646  0.6812  2.0680
     4   1  12   1 UN    NA  0.9246  0.0533  0.0221  0.0487  -1  0.843179  0.9999  2.6494
     3   1  12   1 UN    NA  0.8884  0.1116  0.0000  0.0558  -1  0.836924  0.9148  2.2079
     9   1  20   1 UN    NA  0.8871  0.1129  0.0000  0.0565  -1  0.843400  0.9994  2.5373
     4   1   5   1 UN    NA  0.9206  0.0454  0.0340  0.0567  -1  0.844560  1.0000  2.7328
     1   1  10   1 UN    NA  0.8979  0.0770  0.0251  0.0636  -1  0.845403  0.9995  2.5405
    13   1  17   1 UN    NA  0.9239  0.0000  0.0761  0.0761  -1  0.844079  1.0000  2.6842
    14   1  20   1 UN    NA  0.9152  0.0000  0.0848  0.0848  -1  0.844971  0.9999  2.6613
    10   1  20   1 UN    NA  0.8330  0.1618  0.0052  0.0861  -1  0.848307  1.0000  2.7613
    17   1  20   1 UN    NA  0.8985  0.0166  0.0850  0.0932  -1  0.850837  1.0000  3.0909
     5   1  20   1 UN    NA  0.8947  0.0084  0.0968  0.1010  -1  0.852193  1.0000  3.0216
    15   1  20   1 UN    NA  0.8979  0.0000  0.1021  0.1021  -1  0.851406  1.0000  3.9734
     1   1  12   1 UN    NA  0.8888  0.0000  0.1112  0.1112  -1  0.853173  1.0000  3.6029
     1   1  15   1 UN    NA  0.8888  0.0000  0.1112  0.1112  -1  0.848272  1.0000  3.1053
     5   1  13   1 UN    NA  0.8869  0.0000  0.1131  0.1131  -1  0.847869  1.0000  4.1703
     1   1  13   1 UN    NA  0.8163  0.1343  0.0494  0.1165  -1  0.853563  1.0000  3.8698
     1   1   4   1 UN    NA  0.7657  0.2343  0.0000  0.1172  -1  0.850614  1.0000  3.2523
     1   1  17   1 UN    NA  0.8824  0.0000  0.1176  0.1176  -1  0.848207  0.9628  2.2744
    12   1  13   1 UN    NA  0.8675  0.0174  0.1151  0.1238  -1  0.855856  1.0000  4.8696
     4   1  13   1 UN    NA  0.7484  0.2516  0.0000  0.1258  -1  0.845718  1.0000  2.9270
    12   1  20   1 UN    NA  0.8218  0.1026  0.0756  0.1269  -1  0.855566  1.0000  3.7222
    13   1  15   1 UN    NA  0.8718  0.0000  0.1282  0.1282  -1  0.853549  1.0000  4.1823
     1   1   3   1 UN    NA  0.7431  0.2569  0.0000  0.1284  -1  0.846335  0.9993  2.5354
     1   1   5   1 UN    NA  0.8701  0.0000  0.1299  0.1299  -1  0.849784  1.0000  3.1839
     9   1  14   1 UN    NA  0.8625  0.0000  0.1375  0.1375  -1  0.853595  1.0000  2.8608
     9   1  10   1 UN    NA  0.7908  0.1377  0.0715  0.1404  -1  0.857452  1.0000  3.4369
     3   1  14   1 UN    NA  0.7695  0.1678  0.0627  0.1466  -1  0.858188  1.0000  3.5224
    13   1  20   1 UN    NA  0.7115  0.2828  0.0057  0.1471  -1  0.857202  1.0000  4.2191
     3   1  13   1 UN    NA  0.7053  0.2947  0.0000  0.1474  -1  0.839575  0.9993  2.5315
     4   1  16   1 UN    NA  0.8512  0.0000  0.1488  0.1488  -1  0.859283  1.0000  2.7114
     9   1  17   1 UN    NA  0.8454  0.0000  0.1546  0.1546  -1  0.859423  1.0000  2.6802
     4   1  14   1 UN    NA  0.7870  0.0917  0.1213  0.1671  -1  0.862278  1.0000  3.6131
     4   1  20   1 UN    NA  0.6619  0.3381  0.0000  0.1691  -1  0.852609  1.0000  3.5450
    10   1  17   1 UN    NA  0.8284  0.0033  0.1683  0.1700  -1  0.863567  1.0000  3.0617
    10   1  16   1 UN    NA  0.8265  0.0000  0.1735  0.1735  -1  0.855147  0.9997  2.5681
     4   1   9   1 UN    NA  0.6713  0.3038  0.0249  0.1768  -1  0.861883  1.0000  3.5300
     3   1  20   1 UN    NA  0.6363  0.3637  0.0000  0.1818  -1  0.848857  1.0000  2.7975
     3   1  17   1 UN    NA  0.7112  0.2044  0.0844  0.1866  -1  0.864419  1.0000  3.8245
     3   1  10   1 UN    NA  0.6180  0.3820  0.0000  0.1910  -1  0.861797  1.0000  3.5829
     4   1  10   1 UN    NA  0.6425  0.3087  0.0487  0.2031  -1  0.866156  1.0000  4.0055
    14   1  17   1 UN    NA  0.7954  0.0000  0.2046  0.2046  -1  0.860049  1.0000  2.7085
     1   1  20   1 UN    NA  0.6532  0.2196  0.1272  0.2370  -1  0.872562  1.0000  5.1176
     3   1  11   1 UN    NA  0.7456  0.0000  0.2544  0.2544  -1  0.872914  1.0000  3.8229
    10   1  14   1 UN    NA  0.7334  0.0000  0.2666  0.2666  -1  0.877077  1.0000  4.6890
     3   1   4   1 UN    NA  0.4571  0.5429  0.0000  0.2714  -1  0.865907  1.0000  4.6335
     9   1  16   1 UN    NA  0.6517  0.0814  0.2668  0.3075  -1  0.885446  1.0000  4.3295
     3   1  16   1 UN    NA  0.5426  0.2396  0.2178  0.3376  -1  0.888910  1.0000  5.5282
     3   1   9   1 UN    NA  0.3717  0.4889  0.1394  0.3838  -1  0.894175  1.0000  8.6458
     4   1  17   1 UN    NA  0.5368  0.1547  0.3084  0.3858  -1  0.897618  1.0000  6.8319
   FID1 IID1 FID2 IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO



   plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --double-id \
       --allow-extra-chr --set-missing-var-ids @:# --out out --recode vcf --out britomartisSWEDEN.


#Stats for modern population

plink --bfile britomartisSWEDEN.modern --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 10 5 0.7 --out britomartismo

plink --bfile britomartisSWEDEN.modern  --double-id \
--allow-extra-chr --set-missing-var-ids @:# --extract britomartismo.prune.in  --recode vcf --out britomartisSWEDEN.modern


# Calculate nucleotide diversity (Pi)
vcftools --vcf britomartisSWEDEN.modern.vcf --site-pi --out pi_output
vcftools --vcf britomartisSWEDEN.modern.vcf --window-pi 50000 --out pi.modernworld
vcftools --vcf britomartisSWEDEN.modern.vcf --TajimaD 50000 --out taj.modernworld
vcftools --vcf britomartisSWEDEN.modern.vcf --relatedness --out relate.modernworld

q

mean_pi=$(awk '{sum += $4} END {print sum/NR}' taj.modernworld.Tajima.D)
echo "Tajimas D: $mean_pi"

Tajimas D: -0.537471


# Calculate the mean of the 5th column
mean_pi=$(awk '{sum += $3} END {print sum/NR}' pi_output.sites.pi)

# Output the result
echo "mean pi: $mean_pi"

mean pi: 0.206222


#Stats for historical population

plink --bfile britomartisALLIND.hist  --double-id \
--allow-extra-chr --set-missing-var-ids @:#  --recode vcf --out britomartisALLIND

vcftools --vcf britomartisALLIND.vcf --site-pi --out pi_site.hist
vcftools --vcf britomartisALLIND.vcf --window-pi 50000 --out pi.hist
vcftools --vcf britomartisALLIND.vcf --TajimaD 1000000 --out taj.hist
vcftools --vcf britomartisALLIND.vcf --relatedness --out relate.hist


vcftools --vcf britomartisALLIND.vcf --TajimaD 1000 --out taj.hist
awk '{if ($3 != 'nan') {sum += $3; count++}} END {if (count > 0) print sum/count; else print "No non-zero values"}' taj.hist.Tajima.D

Tajimas D: 1.11864


# Calculate the mean of the 5th column
mean_pi=$(awk '{sum += $3} END {print sum/NR}' pi_site.hist.sites.pi)

# Output the result
echo "mean pi: $mean_pi"

mean pi: 0.430733

#Modern Sweden

plink --tped angsd_modern_plink.tped --tfam angsd_modern_plink.tfam --geno 0.1 --maf 0.05 --remove remove_world.txt  --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.3 --out britomartis.modernswed

plink --tped "$TPED_FILE" --tfam "$TFAM_FILE" --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.modernall.prune.in  \
    --make-bed --remove remove.txt --out britomartisNOALT2.modern

plink --tped angsd_modern_plink.tped --tfam angsd_modern_plink.tfam  --geno 0.1 --maf 0.05 --double-id --pca \
    --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.modernall.prune.in  \
    --make-bed --remove remove.txt --recode vcf --out britomartisNOALT2.modern


plink --tped angsd_modern_plink.tped --tfam angsd_modern_plink.tfam --geno 0.1 --maf 0.05 --recode vcf --remove remove_world.txt  --double-id --allow-extra-chr --set-missing-var-ids @:# --extract britomartis.modernall.prune.in --out britomartis.modernswed

vcftools --vcf britomartis.modernswed.vcf --site-pi --out piMS_site.hist
vcftools --vcf britomartis.modernswed.vcf --window-pi 50000 --out piMS.hist
vcftools --vcf britomartis.modernswed.vcf --TajimaD 1000 --out tajMS.hist
vcftools --vcf britomartis.modernswed.vcf --relatedness --out relateMS.hist


vcftools --vcf britomartisALLIND.vcf --TajimaD 1000 --out taj.hist
awk '{if ($3 != 'nan') {sum += $3; count++}} END {if (count > 0) print sum/count; else print "No non-zero values"}' tajMS.hist.Tajima.D

Tajimas D: 1.11864


# Calculate the mean of the 5th column
mean_pi=$(awk '{sum += $3} END {print sum/NR}' piMS_site.hist.sites.pi)

# Output the result
echo "mean pi: $mean_pi"

mean pi: 0.430733

#Relatedness analysis

plink --bfile britomartisALLIND.hist --double-id \
    --allow-extra-chr --set-missing-var-ids @:# --genome --outbritomartisALLIND.hist.relate

plink --bfile britomartisNOALT2.modern  --double-id --allow-extra-chr --set-missing-var-ids @:# --genome --out britomartisNOALT2.modern.relate
