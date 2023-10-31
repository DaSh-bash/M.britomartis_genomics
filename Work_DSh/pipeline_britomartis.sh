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
  WARNING: Skipping mount /var/apptainer/mnt/session/etc/resolv.conf [files]: /etc/resolv.conf doesn't exist in container
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


 bcftools filter -e FORMAT/PL[0] < 50'  reduced_samle_contemp_mcalling.Chr1.qfilter.vcf.gz


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
