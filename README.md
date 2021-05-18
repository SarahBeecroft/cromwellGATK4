# Local Cromwell implementation of GATK4 germline variant calling pipeline
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


4. Upload the resource bundle file from IRDS using rclone or filezilla and unpack it with `tar xzvf resource.tar.gz`

5. Set up the config files. Files that you need to edit with the correct paths to your data/jar files or other specific configurations are:
    - `Multisample_Fastq_to_Gvcf_GATK4_inputs_hg38.json`
    - `Multisample_jointgt_GATK4_inputs_hg38.json`
        - both json files will need the correct paths to you picard and GATK jar files
    - `samples.txt`
    - `gvcfs.txt`
        - These are the sample input files (tab seperated)
        - The format for samples.txt is sample ID, fastq R1 location and file name, fastq R2 location and file name
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

6. Launch the job within a `screen` or `tmux` session, using `./launch_cromwell.sh`. When that has completed successfully, you can launch the second stage of the pipeline (joint calling) with `./launch_jointgt.sh`

### Overview of the steps in the pipeline

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

### Dependencies
    -BWA/0.7.15
    -GATK >=4.0.6.0
    -SAMtools/1.5
    -picard/2.9
    -Python/2.7
    -Cromwell
