#!/bin/bash

conda activate gatk4_pipeline

java -Dconfig.file=/data/cromwellGATK4/local.conf -jar /data/cromwellGATK4/cromwell-61.jar run -i Multisample_jointgt_GATK4_inputs_hg38.json
