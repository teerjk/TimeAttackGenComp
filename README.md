# TimeAttackGenComp #

## Fast comparison of sample genotypes for quality control ##

- - - -

### Purpose ###
Massively parallel sequencing involves many sample handling steps leading to the potential for inappropriate sample identity. This can result from sample swaps, sample mixing, or similar. This tool is designed to compare genotypes at defined positions (ideally, common SNPs). The level of genotype discordance can be used to infer sample relationships. Low discordance between a pair of samples indicates samples coming from the same individual; higher discordance indicates samples come from separate individuals. This measured genetic discordance can be compare to expected sample relationships to confirm sample pairs, or ensure all samples are from unique individuals.
A workflow is provided to start from FASTQ files, perform fast alignment with SNAP, call variants with BCFtools, and compare with a fast perl script. Several additional plotting features are added to assist with QC evaluation.

### Prerequisites ### 
1. Cromwell (tested with version 38): https://github.com/broadinstitute/cromwell
1. SNAP (tested with verion 1.0.3 and 1.0dev.104): https://github.com/amplab/snap
1. BCFTools (tested with version 1.9): https://github.com/samtools/bcftools
1. Java (to run cromwell)
1. Perl (to run comparison)
1. R (for plotting)

### Getting started ###
To run the end-to-end pipeline, a WDL script is provided. You will need to prepare a sample input file, fill in details in the input.json file, and ensure 'snap-aligner' and 'bcftools' are in your PATH. You may want to edit the WDL script to adjust resources defined in each "runtime" block.

1. Download repository to a known location.
1. Copy WDL script and input.json to working directory and change to working dir.
1. Create a sample input file: a tab-delimited text file, one sample per line with the following columns
    1. sample_name
    1. FASTQ_read1_full_path
    1. FASTQ_read2_full_path
1. Edit input.json file:
    1. output_name: define an output name for the summary files
    1. inputFastqFile: the full path to the input text file from step 3.
    1. output_dir: full path to directory for output files (will be created if it doesn't exist).
    1. target_bed: A BED file based on the appropriate genome reference defining positions to compare. A file based on GRCh37 coordinates within refGene regions that have 1000 Genomes allele frequency >15% has been provided for convenience.
    1. ref: define the folder for SNAP reference indexes
    1. ref_fastq: full path to reference multi FASTA. A .fai file should be in the same folder.
    1. compare: Adjust the path to point to your local copy of compare_simple.pl
    1. R_heatmap: Adjust the path to point to your local copy of heatmap.R
    1. plot_dist: Adjust the path to point to your local copy of plot_dist.R
1. Submit the WDL file to cromwell. An example bash wrapper is provided (adjust to point to your configuration and cromwell .jar file.
1. When satisfied, delete the temporary files in cromwell-executions/ and cromwell-workflow-logs/

### Interpreting Output ###
1. The main output is OUTPUT_NAME.snv.out.txt. This file is a tab-delim text showing the all-vs-all matrix of sample-sample genotype discordance. This file can be loaded into Excel for visualization.
2. A simple heatmap is plotted: OUTPUT_NAME.pdf. The output order is based on your inputFASTQFile. If samples expected to come from the same individual are listed together, the heatmap will show triangles of red low discordance along the diagonal. If no samples should match, you should only see red boxes on the diagonal (samples matching themselves).
3. VCF files with SNV calls are written to vcf/
4. Allele frequency files, histogram plots, and scatterplots colored by chromosome are included in af/

### Caveats ###
We have noticed that tumor samples with extensive Loss Of Heterozygoisty (LOH) can show high discordance with a matched normal. This is because we compare bi-allelic genotypes, and a heterozgous genotype in the normal would not match an LOH homozygote in the tumor. Viewing the af.chr.pdf allele frequency scatterplots is recommended to check for LOH.

We have also noticed that samples with low levels of sample mixing may not have high discrepancy when compared to another uncontaminated sample from the same individual. This is due to the calling of a genotype: low levels of contamination from another sample may not affect genotyping (by design to avoid error from contamination.) Viewing the allele frequency plots may be useful: look for 'pinched' allele frequency bands at ~5-10% and ~90-95%. If seen in all chromsomes, this suggests sample mixing, although other experiments should be run to confirm.

### Todo ###
1. Allow for running the workflow with BAM, VCF, or snv.smp input files to avoid alignment and genotype calling if not needed.
1. Address false mismatch due to LOH.
1. Prepare Docker containers to use within WDL.

### Citation ###
Preprint coming soon!
