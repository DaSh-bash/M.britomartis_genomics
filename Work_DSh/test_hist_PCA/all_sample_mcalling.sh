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
