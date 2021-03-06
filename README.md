# Local Cromwell implementation of GATK4 germline variant calling pipeline
See the [GATK](https://gatk.broadinstitute.org/hc/en-us) website for more information on this toolset 
## Assumptions
- Using hg38 human reference genome build
- Running 'locally' i.e. not using HPC/SLURM scheduling, or containers. This repo was specifically tested on Pawsey Nimbus 16 CPU, 64GB RAM virtual machine, primarily running in the `/data` volume storage partition. 
- Starting from short-read Illumina paired-end fastq files as input

### Dependencies
The following versions have been tested and work, but GATK and Cromwell are regularly updated and so one must consider whether they would like to use newer versions of these tools. 
- BWA/0.7.15
- GATK v4.0.6.0
- SAMtools/1.5
- picard/2.9
- Python/2.7
- Cromwell v61

## Quick start guide
### Installing and preparing environment for GATK4 with Cromwell

1. Clone repository
```
git clone https://github.com/SarahBeecroft/cromwellGATK4.git
cd cromwellGATK4
chmod 777 *.sh
```

2. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you haven’t already. Create Conda environment using the supplied conda environment file

```
conda env create --file gatk4_pipeline.yml
```

3. Download the necessary .jar files
    - The Cromwell workfow orchestration engine can be downloaded from https://github.com/broadinstitute/cromwell/releases/ 
    - GATK can be downloaded from https://github.com/broadinstitute/gatk/releases. Unzip the file with `unzip` 
    - Picard can be downloaded from https://github.com/broadinstitute/picard/releases/


4. Upload the resource bundle file from IRDS using rclone or filezilla and unpack it with `tar xzvf resource.tar.gz`. Note that the `hg38_wgs_scattered_calling_intervals.txt` will need to be to generated using the following:

```
cd <your_resource_dir>
find `pwd` -name "scattered.interval_list" -print | sort > hg38_wgs_scattered_calling_intervals.txt
```

5. Set up the config files. Files that you need to edit with the correct paths to your data/jar files or other specific configurations are:
    - `Multisample_Fastq_to_Gvcf_GATK4_inputs_hg38.json`
    - `Multisample_jointgt_GATK4_inputs_hg38.json`
        - both json files will need the correct paths to your reference file locations, and the file specifying your inputs i.e. `samples.txt` or `gvcfs.txt`
    - `samples.txt`
    - `gvcfs.txt`
        - These are the sample input files (tab seperated)
        - The format for samples.txt is sampleID, sampleID_readgroup, path_to_fastq_R1_file, path_to_fastq_R2_file,
        - The format for gvcfs.txt is sample ID, gvcf, gvcf .tbi index file
        - Examples are included in this repo
        - NOTE: Having tabs, not spaces, is vital for parsing the file. Visual studio code tends to introduce spaces, so if you are having issues, check the file with another text editor such as sublime. 
    - `launch_cromwell.sh`
    - `launch_jointgt.sh`
        - These are the scripts which launch the pipeline. 
        - `launch_cromwell.sh` launches the fastq to gvcf stage
        - `launch_jointgt.sh` launched the gvcf joint genotyping to cohort vcf step. This is perfomed when you have run all samples through the fastq to gvcf stage.
        - Check the paths and parameters make sense for your machine
    - `local.conf`
        - the main tuneable parameters here are:
        	- `concurrent-job-limit = 5` this is the max number of concurrent jobs that can be spawned by cromwell. This depends on the computational resources available to you. 5 was determined to work reasonably well on a 16 CPU, 64GB RAM Nimbus VM (Pawsey). 
        	- `call-caching enabled = true`. Setting this parameter to `false` will disable call caching (i.e. being able to resume if the job fails before completion). By default, call caching is enabled. 
    - `cromwell.options`
        - `cromwell.options` requires editing to provide the directory where you would like the final workflow outputs to be written
    - `Multisample_Fastq_to_Gvcf_GATK4.wdl`
    - `ruddle_fastq_to_gvcf_single_sample_gatk4.wdl`
        - The paths to your jar files will need to be updated
        - The path to your conda `activate` binary will need to be updated (e.g. `/data/miniconda/bin/activate`)

6. Launch the job within a `screen` or `tmux` session, using `./launch_cromwell.sh`. When that has completed successfully, you can launch the second stage of the pipeline (joint calling) with `./launch_jointgt.sh`. Ensure you pipe the stdout and stderr to a log file using (for example) `./launch_cromwell.sh &> cromwell.log`

### Overview of the steps in `Multisample_Fastq_to_Gvcf_GATK4.wdl`
This part of the pipeline takes short-read, Illumina paired-end fastq files as the input. The outputs generated are sorted, duplicate marked bam files and their indices, duplicate metric information, and a GVCF file for each sample. The GVCF files are used as input for the second part of the pipeline (joint genotyping).

```
FastqToUbam
GetBwaVersion
SamToFastqAndBwaMem
MergeBamAlignment
SortAndFixTags
MarkDuplicates
CreateSequenceGroupingTSV
BaseRecalibrator
GatherBqsrReports
ApplyBQSR
GatherBamFiles
HaplotypeCaller
MergeGVCFs
```

### Overview of the steps in `Multisample_jointgt_GATK4.wdl`
This part of the pipeline takes GVCF files (one per sample), and performs joint genotyping across all of the provided samples. This means that old previously generated GVCFs can be joint-called with new GVCFs whenever you need to add new samples. The key output from this is a joint-genotyped, cohort-wide VCF file. This file can be used for a GEMINI database after normalisation with VT and annotation with a tool such as VEP or SNPEFF. 

The file `hg38.custom_100Mb.intervals` is required for this step of the pipeline to run. This is included in the git repo for convenience, but should be moved to your resource directory with all the other resource files. 

```
GetNumberOfSamples
ImportGVCFs
GenotypeGVCFs
HardFilterAndMakeSitesOnlyVcf
IndelsVariantRecalibrator
SNPsVariantRecalibratorCreateModel
SNPsVariantRecalibrator
GatherTranches
ApplyRecalibration
GatherVcfs
CollectVariantCallingMetrics
GatherMetrics
DynamicallyCombineIntervals
```
