#!/bin/bash

conda activate gatk4_pipeline

java -Dconfig.file=/data/cromwellGATK4/local.conf -jar /data/cromwellGATK4/cromwell-61.jar run /data/cromwellGATK4/Multisample_Fastq_to_Gvcf_GATK4.wdl \
   -i /data/cromwellGATK4/Multisample_Fastq_to_Gvcf_GATK4_inputs_hg38.json \
   -o /data/cromwellGATK4/cromwell.options
