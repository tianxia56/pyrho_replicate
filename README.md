# pyrho_replicate


## run process*.sh {pop} {chr}
## test run*demo.sh {pop} {chr}
conda activate smcpp_env
module load BCFtools/1.11-GCC-10.2.0
module load SAMtools/1.21-GCC-12.2.0

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/
20141020.strict_mask.whole_genome.bed

Alignment with Spence & Song (2019) for SMC++ and PyRhO Preparation:
The paper's workflow involves:
	1. Data QC & Phasing: (User mentioned this is done prior to these scripts)
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
