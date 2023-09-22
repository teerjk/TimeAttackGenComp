# TimeAttackGenComp #

## Fast comparison of sample genotypes for quality control ##

- - - -

### Purpose ###
Massively parallel sequencing involves many sample handling steps leading to the potential for inappropriate sample identity. This can result from sample swaps, sample mixing, or similar. This tool is designed to compare genotypes at defined positions (ideally, common SNPs). The level of genotype discordance can be used to infer sample relationships. Low discordance between a pair of samples indicates samples coming from the same individual; higher discordance indicates samples come from separate individuals. This measured genetic discordance can be compare to expected sample relationships to confirm sample pairs, or ensure all samples are from unique individuals.
A WDL workflow is provided to start from FASTQ, BAM, or VCF files. When starting from FASTQ, fast alignment is performed with SNAP. Variants are then called with BCFtools, and compared with a fast Perl script. Several additional plotting features are included to assist with QC evaluation.

### Prerequisites ### 
1. Cromwell (tested with versions 38, 42): https://github.com/broadinstitute/cromwell
1. SNAP (tested with version 1.0.3, 1.0dev.104, and 2.0.1): https://github.com/amplab/snap
1. BCFTools (tested with version 1.9, 1.15.1): https://github.com/samtools/bcftools
1. Java (to run cromwell)
1. Perl (to run comparison)
1. R (for plotting)

SNAP genome references indexes need to be built before use. Check SNAP documentation, but this should work:
```
snap-aligner index ref.fa .
```


### Getting started ###
To run the end-to-end pipeline, a WDL script is provided. You will need to prepare a sample input file, fill in details in the input.json file, and ensure 'snap-aligner' and 'bcftools' are installed and in your PATH. (See below to use the provided docker file.) You may want to edit the WDL script to adjust resources defined in each "runtime" block and remove the "docker" blocks if you are not using the container.

1. Download repository to a known location.
1. Copy TimeAttackGenComp.wdl and TimeAttackGenComp_inputs.json to working directory and change to working dir.
1. Create an input file based on desired input: FASTQ, BAM, or VCF
    1. Create a FASTQ sample input file: a tab-delimited text file, one sample per line with the following columns
        1. sample_name
        1. FASTQ_read1_full_path
        1. FASTQ_read2_full_path
    1. OR: create a BAM sample input file: a tab-delimited text file, one sample per line with the following columns
        1. sample_name
        1. BAM_full_path
        1. BAI_full_path
    1. OR: create a VCF sample input file: a tab-delimited text file, one sample per line with the following columns
        1. sample_name
        1. BGzipped_VCF_full_path
        1. VCF.tbi_full_path
1. Edit TimeAttackGenComp_inputs.json file:
    1. Select the desired input block (FASTQ, BAM, or VCF) and delete the other two.
    1. output_name: define an output name for the summary files
    1. output_dir: full path to directory for output files (will be created if it doesn't exist).
    1. inputFile: the full path to the input text file from step 3.
    1. target_bed: An uncompressed BED file based on the appropriate genome reference defining positions to compare. Files based on GRCh37 and GRCh38 (lifted over) coordinates within refGene regions that have 1000 Genomes allele frequency >15% have been provided for convenience. Uncompress the included files if you wish to use them.
    1. ref: define the folder for SNAP reference indexes, which need to be built before running.
    1. ref_fastq: full path to reference multi FASTA. A .fai file should be in the same folder.
    1. compare: Adjust the path to point to your local copy of compare_simple.pl
    1. R_heatmap: Adjust the path to point to your local copy of heatmap.R
    1. plot_dist: Adjust the path to point to your local copy of plot_dist.R
1. Submit the WDL file to cromwell. An example bash wrapper is provided (adjust to point to your configuration and cromwell .jar file.
1. When satisfied, delete the temporary files in cromwell-executions/ and cromwell-workflow-logs/


### Using Docker ###
A Dockerfile is provided that includes SNAP, samtools, and downloads the latest version of this TimeAttackGenComp repository. We have successfully built with Docker and converted to Singularity .sif files. You can then define the docker/singularity image in the runtime blocks of TimeAttackGenComp.wdl. (This is currently defined as "timeattackgencomp_0.2", but you may need to change depending on how you define your images.) You will need to define your cromwell application configuration file to appropriately handle containers in your environment. See the [Cromwell Container documentation](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/) for details.

An example cromwell application configuration file "application.slurm.singularity.conf" is provided.
Substitute the following strings as appropriate:
  SINGULARITY_EXAMPLE  #substitute with the full path to your SIF file timeattackgencomp_0.2.sif
  EXAMPLE_DATA         #substitute with the full path to folder to be mounted. See Singularity documentation for details on the --bind command. Remove '/EXAMPLE_DATA:/EXAMPLE_DATA:ro,' if no other folder are needed to be mounted.

The Slurm endpoint was succesfully tested in our environment. The PBS endpoint has worked in the past, but should be considered experimental as we no longer have PBS. We recommend this configuration file is used as inspiration, as HPC configurations vary. Note that we use "pbs_walltime" to define walltime resource requests.

### Interpreting Output ###
1. The main output is OUTPUT_NAME.snv.out.txt. This file is a tab-delim text showing the all-vs-all matrix of sample-sample genotype discordance. This file can be loaded into Excel for visualization.
2. A simple heatmap is plotted: OUTPUT_NAME.pdf. The output order is based on your inputFASTQFile. If samples expected to come from the same individual are listed together in the input file, the heatmap will show triangles of red low discordance along the diagonal. If no samples should match, you should only see red boxes on the diagonal (samples matching themselves).
3. VCF files with SNV calls are written to vcf/
4. Allele frequency files, histogram plots, and scatterplots colored by chromosome are included in af/

### Caveats ###
We have noticed that tumor samples with extensive Loss Of Heterozygoisty (LOH) can show high discordance with a matched normal. This is because we compare bi-allelic genotypes, and a heterozgous genotype in the normal would not match an LOH homozygote in the tumor. Viewing the af.chr.pdf allele frequency scatterplots is recommended to check for LOH.

We have also noticed that samples with low levels of sample mixing may not have high discrepancy when compared to another uncontaminated sample from the same individual. This is due to the calling of a genotype: low levels of contamination from another sample may not affect genotyping (by design to avoid error from contamination.) Viewing the allele frequency plots may be useful: look for 'pinched' allele frequency bands at ~5-10% and ~90-95%. If seen in all chromsomes, this suggests sample mixing, although other experiments should be run to confirm.

### Todo ###
1. Allow for running the workflow with VCF or snv.smp input files to avoid genotype calling if not needed. DONE!
1. Address false mismatch due to LOH.
1. Prepare Docker containers to use within WDL. DONE!

### Citation ###
Hopefully coming soon!
