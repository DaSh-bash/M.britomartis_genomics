# Some sh commands to start the different pipelines


############# SNP CALLING AND FILTERING USING SBATCH RUN FILES #################
# Note the code is cleaned from specific uppamx commands, so one have to add
# "module load XXX" depending on software, version etc

# Reference fasta:
fasta="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa"


# MAP & MARK DUPLICATE
for ind in $(cat all_individuals.txt)
do
 R1="fastq/"$id"_1.fastq.gz"
 R2="fastq/"$id"_2.fastq.gz"
 bampref="bam/$ind"
 map=$(sbatch -J bwa.$ind -t 1-00:00:00 -p node -n 20 run_bwa_mem.sh $fasta $R1 $R2 $ind $bampref |cut -f4 -d" ")
 sbatch -J mkdup.$ind -d afterok:$map -n 10 -t 1-00:00:00 -p core run_mkdupl.sh $ind $bampref
done

# BASE RECALIBRATE NEW SAMPLES
vcf1="SNPs.Kardos.vcf.gz"
vcf2="dogs.557.publicSamples.ann.chrAll.PASS.vcf.gz"
for ind in $(cat new_individuals.txt)
do
 bampref="bam/$ind.md"
 sbatch -J bqsr.$ind -t 1-00:00:00 -p node -n 20 run_bqsr.sh $fasta $vcf1 $vcf2 $bampref
done

# HAPLOTYPECALLER
for bam in $(cat all_bamfiles.txt)
do
 ind=`echo $bam |cut -f2 -d"/" |cut -f1 -d"."`
 vcfpref="gvcf/$ind"
 sbatch -J hc.$ind -p node -n20 -t 1-00:00:00 -p core run_haplotypeCaller.sh $fasta $bam $vcfpref.hc
done

# JOINT GENOTYPING
prefix="AllInd"
out="vcf/$prefix.vcf"
in="$prefix.gvcf.list"
joint=$(sbatch -J JointGenotyping -t 7-00:00:00 -C mem1TB -p node -n20 run_GenotypeGVCFs.sh $fasta $in $out |cut -f4 -d" ")
sbatch -J tabix -d afterok:$joint -t 3-00:00:00 -p core run_zip_and_tabix.sh $out



# EXTRACT ONLY SNPs and USE GATKs HARD FILTER
in="vcf/$prefix.vcf.gz"
outpref="vcf/$prefix"
filt=$(sbatch -J SNPHardFilt.$prefix -t 5-00:00:00 -p core -n10 run_GATK_hard_filter.sh $ref $in $outpref |cut -f4 -d" ")
sbatch -J tabixSNP -d afterok:$filt -t 1-00:00:00 -p core run_zip_and_tabix.sh $outpref.SNPs.vcf
sbatch -J tabixHF -d afterok:$filt -t 1-00:00:00 -p core run_zip_and_tabix.sh $outpref.SNPs.HF.vcf


# SAVE ONLY BIALLELIC
in="vcf/$prefix.SNPs.HF.vcf.gz"
outp="vcf/$prefix.SNPs.HF.Bi"
biall=$(sbatch -J biall -t 2-00:00:00 -p core run_extractBiAllelic.sh $in $outp |cut -f4 -d" ")
sbatch -J tabixBiAll -d afterok:$biall -t 1-00:00:00 -p core run_zip_and_tabix.sh $outp.recode.vcf


# REMOVE LOCI WITH ONLY HET OR HOM CALLS
in="vcf/$prefix.SNPs.HF.Bi.recode.vcf.gz"
outp="vcf/$prefix.SNPs.HF.Bi.HH"
remove="$prefix.OnlyHomOrHet.txt"
rmhomhet=$(sbatch -J rmHomHet -t 2-00:00:00 -p core run_removeOnlyHomOrHet.sh $in $remove $outp |cut -f4 -d" ")
sbatch -J tabixrmHomHet -d afterok:$rmhomhet -t 1-00:00:00 -p core run_zip_and_tabix.sh $outp.recode.vcf


# FILTER FOR COVERAGE
sbatch -J CheckCov -t 10:00:00 -p core $scrdir/run_coverageFromVCF.sh $prefix.vcf.gz $prefix.vcf.idepth
# Check the average mean depth:
awk '(NR>1){n++; sum+=$3}END{mean=sum/n; print mean}' $prefix.vcf.idepth
#24.7207
maxmeandepth=49.4414
minmeandepth=10
in="../vcf/$prefix.SNPs.HF.Bi.HH.recode.vcf.gz"
outp="../vcf/$prefix.SNPs.HF.Bi.HH.DP"
depth=$(sbatch -J depth -t 2-00:00:00 -p core run_coverageFilter.sh $in $minmeandepth $maxmeandepth $outp |cut -f4 -d" ")
sbatch -J tabixDepth -d afterok:$depth -t 1-00:00:00 -p core $scrdir/run_zip_and_tabix.sh $outp.recode.vcf


# PHREDSCORE AND MISSING DATA
GQ=30
maxmiss=0.95
in="vcf/$prefix.SNPs.HF.Bi.HH.DP.recode.vcf.gz"
outp="vcf/$prefix.SNPs.HF.Bi.HH.DP.GQ.miss"
phred=$(sbatch -J phredFilt -t 10:00:00 -p core run_filterPhredAndMissing.sh $in $GQ $maxmiss $outp |cut -f4 -d" ")
sbatch -J tabixPhred -d afterok:$phred -t 5:00:00 -p core run_zip_and_tabix.sh $outp.recode.vcf

# MINOR ALLELE COUNT
macval=2
in="vcf/$prefix.SNPs.HF.Bi.HH.DP.GQ.miss.recode.vcf.gz"
outp="vcf/$prefixSNPs.HF.Bi.HH.DP.GQ.miss.mac"
mac=$(sbatch -J MAC$macval -t 5:00:00 -p core run_filterMAC.sh $in $macval $outp |cut -f4 -d" ")
sbatch -J tabixMAC$macval -d afterok:$mac -t 5:00:00 -p core run_zip_and_tabix.sh $outp.recode.vcf

# HARDY (excess of heterozygous individuals)
hardy=0.001
in="vcf/$prefix.SNPs.HF.Bi.HH.DP.GQ.miss.mac.recode.vcf.gz"
out1="vcf/$prefixSNPs.HF.Bi.HH.DP.GQ.miss.mac"
out2="vcf/$prefix.final"
har=$(sbatch -J hardy -t 5:00:00 -p core $scrdir/run_filterHardy.sh $in $hardy $out1 $out2 |cut -f4 -d" ")
sbatch -J tabixHardy$hardy -d afterok:$har -t 5:00:00 -p core $scrdir/run_zip_and_tabix.sh $out2.recode.vcf





################################## SNAKEMAKE ###################################
smkdir="./snakemake/"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SNAKEMAKE PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Runs joint genotyping, SNP extraction, filtering and VEP
#Make sure the following programs are available on the system:
#snakemake/5.30.1
#GATK/3.8-0
#htslib/1.12
#vcftools/0.1.16
#vep/99
# Plot the DAG (directed acyclic graph)
snakemake --snakefile $smkdir/main.snakefile --dag | dot -Tsvg > dag/dag.$p1.svg
# Run full pipeline
snakemake --snakefile $smkdir/main.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -C {cluster.C} -t {cluster.time} -e {cluster.error} -o {cluster.output}" --cluster-config $smkdir/cluster_main.json


# ~~~~~~~~~~~~~~~~~~~~~~~~ OUTGROUP SNAKEMAKE PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~
# Pipeline that maps the outgroups onto the dog reference, makes allsites files
# and extracts the relevant sites (polymorphic in our wolves)
# Software needed:
#snakemake/5.30.1
#bwa/0.7.17
#samtools/1.10
#picard/2.23.4
#GATK/3.8-0
#BEDTools/2.29.2

# Plot DAG (directed acyclic graph)
snakemake --dag --snakefile $smkdir/outgroup.snakefile | dot -Tsvg > dag/dag.outgroups.svg
# Then start the actual pipeline
snakemake --snakefile $smkdir/outgroup.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -C {cluster.C} -t {cluster.time} -e {cluster.error} -o {cluster.output} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}" --cluster-config $smkdir/snake/cluster_outgroup.json
