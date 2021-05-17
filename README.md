

## Quick start guide
### Installing and preparing environment for GATK4 with Cromwell

1. Clone repository
```
git clone https://github.com/SarahBeecroft/cromwellGATK4.git
cd cromwellGATK4
chmod 777 *.sh
chmod 777 *.pl
```

2. Install [https://docs.conda.io/en/latest/miniconda.html]Miniconda if you haven’t already. Create Conda environment using the supplied conda environment file

```
conda env create --file gatk4_pipeline.yml
```

3. Download the necessary .jar files
    + The Cromwell workfow orchestration engine can be downloaded from https://github.com/broadinstitute/cromwell/releases/ 
    + GATK can be downloaded from https://github.com/broadinstitute/gatk/releases. Unzip the file with `unzip` 
    + Picard can be downloaded from https://github.com/broadinstitute/picard/releases/


4. Upload the resource bundle file from IRDS using rclone or filezilla and unpack it with `tar xzvf resource.tar.gz`

5. Setup the `samples.txt` file, which is a tab-seperated list of file sample names, read group information, and fastq file locations. An example is included in the git repo. For fastq files with the naming convention of `D20_1111_HHWHWDSXY_TTCCTGTT-AGTATCTT_L002_R2.fastq.gz` all in the same directory, the `create_spreadsheet.pl` script can be used to automatically generate the `samples.txt` file. NOTE: Having tabs, not spaces, is vital for parsing the file. Visual studio code tends to introduce spaces, so if you are having issues, check the file with another text editor such as sublime. 

6. Edit the config/pipeline files to reflect the your environment:
    + Within `ruddle_fastq_to_gvcf_single_sample_gatk4.wdl` and `Multisample_jointgt_GATK4.wdl`, the paths to the jar files will need to be checked (default is /data/gatk4_multisample)
    + Edit `cromwell.options` to where you want the final output.
    + Edit `Multisample_Fastq_to_Gvcf_GATK4_inputs_hg38.json` and `Multisample_jointgt_GATK4_inputs_hg38.json` so that the “REFERENCE FILES” and “KNOWN SITES RESOURCES” point to where you copied and extracted the resource.tar.gz file.
    + Check `launch_cromwell.sh` is correct for your machine.
    + Setup the `gvcfs.txt` file, which is a tab-seperated list of file sample names, g.vcf.gz files, and g.vcf.gz.tbi files. 
  
  7. Launch the job within a `screen` or `tmux` session, using `./launch_cromwell.sh`. When that has completed successfully, you can launch the second stage of the pipeline (joint calling) with `./joint_gt.sh`
