#!/bin/bash

source activate gatk4_pipeline

java -Dconfig.file=/data/gatk4_multisample/local.conf -jar /data/gatk4_multisample/cromwell-61.jar run -i Multisample_jointgt_GATK4_inputs_hg38.json
