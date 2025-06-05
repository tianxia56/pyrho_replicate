# pyrho_replicate

extract pop vcf first
then
## run process*.sh {pop} {chr}
## test run*demo.sh {pop} {chr}
conda activate smcpp_env
module load BCFtools/1.11-GCC-10.2.0
module load SAMtools/1.21-GCC-12.2.0

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/
20141020.strict_mask.whole_genome.bed

Alignment with Spence & Song (2019) for SMC++ and PyRhO Preparation:
The paper's workflow involves:
	1. Data QC & Phasing: (User mentioned this is done prior to these scripts):
 

1. 1000 Genomes Project Phase 3 Data: They started with this dataset.
2. Exclusion of Related Individuals: They removed individuals to ensure a set of unrelated samples within each population.
3. Variant Filtering (using VCFtools):
	○ Minimum allele count (MAC): --mac 1 (kept singletons initially, but later steps might effectively remove some).
	○ Minimum minor allele frequency (MAF): --maf 0.001 (minor allele frequency > 0.1%). This removes very rare variants that might be sequencing errors or provide little information for LD-based recombination inference.
	○ Missingness: --max-missing 0.99 (allowed up to 1% missing data per site).
	○ Minimum quality score: --minQ 30 (Phred-scaled quality score for the variant call).
	○ Biallelic SNPs: Restricted to biallelic SNPs.
4. Per-Population Filtering: The above VCFtools filtering was applied independently to each of the 26 populations.
5. Exclusion of Specific Genomic Regions (This is the key "masking" part):
	○ Centromeres: They excluded centromeric regions.
	○ Telomeres: They excluded telomeric regions.
	○ Human Leukocyte Antigen (HLA) region: This region on chromosome 6 is known for its extreme polymorphism, complex LD patterns, and strong selection, which can interfere with standard recombination map inference. They explicitly excluded the HLA region (chr6:28,477,797–33,448,354 in hg19).
Chromosome 8 Inversion Polymorphism: They excluded a large polymorphic inversion region on chromosome 8 (chr8:8,100,000–11,900,000 in hg19), as inversions suppress recombination in heterozygotes and create distinct LD patterns.



 
	2. Ancestral Allele Annotation: Your annotate_aa.py and Step 2 & 3 in process_chromosome_smc.sh handle this well.
	3. Masking: Step 5 in process_chromosome_smc.sh prepares the chromosome-specific mask.
	4. SMC++ for Demography:
		○ Run smc++ vcf2smc (Step 6).
		○ Run smc++ estimate with a uniform recombination rate to get an initial demographic model (model.final.json) (Step 7). This is precisely what your scripts do.
	5. SMC++ for Recombination Maps:
		○ Run smc++ cv using the demographic model from the previous step to infer a fine-scale, population-specific recombination map for each chromosome. This step is NOT currently in your scripts.
	6. PyRhO for Fine-scale Recombination Maps:
		○ Use the VCF (phased, AA-annotated), the demographic model from SMC++, and the recombination map (either from smc++ cv or another source like a sex-averaged map) as input to pyrho.

```
POP="PJL"
mkdir -p "./${POP}_slurm_logs" # Create if you want custom log dir name
sbatch --job-name="${POP}_smc_direct" \
       --output="./${POP}_slurm_logs/${POP}_smc_direct_%j.out" \
       --error="./${POP}_slurm_logs/${POP}_smc_direct_%j.error" \
       ./process_1kgp_smc.sh "${POP}"```

