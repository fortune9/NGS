## Bash/Slurm scripts for calling variants/mutations from DNA-seq data

These scripts start with bam files containing mapped reads to a
genome. The order of running them should be:

	add_RG_and_sort.sbc -> markDup_realign_BQSR.sbc ->
	call_snps.gatk.ploidy.sbc or call_snps.samBcf.sbc

The last step is calling variants using either the
[GATK](https://software.broadinstitute.org/gatk/best-practices/)
method or the
[samtools/bcftools](http://samtools.sourceforge.net/mpileup.shtml)
method.

Run each script without any argument will simply print out a brief
usage. Note that these scripts were written to run as a SLURM job on a
computer cluster, but they can also be run as Linux Bash scripts. So,
for example, 
	
	sbatch add_RG_and_sort.sbc [args]
	bash add_RG_and_sort.sbc [args]

both work.

These scripts are not fully tested, please use with caution. You
feedback is always welcome.

## Author: Zhenguo Zhang, zhangz.sci@gmail.com

