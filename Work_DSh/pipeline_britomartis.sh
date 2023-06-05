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

module load GetOrganelle/1.7.7.0
module load GetOrganelleDB

get_organelle_from_reads.py -1 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_1.fastp.fastq.gz -2 /crex/proj/uppstore2017185/b2014034_nobackup/Dasha/M.britomartis_Conservation/00_Mapping_Calling_sarek/00_Benchmarking4x/SMAL_2_1967-L001_2.fastp.fastq.gz -t 4 -o SMAL_2_1967.mt -F animal_mt -R 10




###### Snakemake for organelle
### working in separate file

mamba env create --name mtgenome --file mapping.yaml
conda activate snakemake-tutorial
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
