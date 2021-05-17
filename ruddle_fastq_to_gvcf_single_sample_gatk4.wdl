## Copyright Broad Institute, 2018
## 
## This WDL pipeline implements data pre-processing according to the GATK Best Practices 
## (June 2016).  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Software version requirements (see recommended dockers in inputs JSON)
## - GATK 4 or later
## - Picard (see gotc docker)
## - Samtools (see gotc docker)
## - Python 2.7
##
## Cromwell version support 
## - Successfully tested on v32
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.


## Adapted to Yale Ruddle HPC by Sander Pajusalu (sander.pajusalu@yale.edu)
##adapted memory specs to fit pawsey system (i.e. 4gb ram per core) and adjusted the JVM heap size to 80% of the available ram. Made any cores be used 100% to prevent sharing of cores. 

# WORKFLOW DEFINITION 
workflow Fastq_to_Gvcf_GATK4 {

  File sample_fastq_file 
  Array[Array[String]] sample_fastqs = read_tsv(sample_fastq_file)
  String sample_name


  String unmapped_bam_suffix
  String ref_name
 
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  File? ref_alt

  String bwa_commandline
  Int compression_level
  
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  File scattered_calling_intervals_list
  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  String base_file_name = sample_name + "." + ref_name


  # Get the version of BWA to include in the PG record in the header of the BAM produced 
  # by MergeBamAlignment. 
  call GetBwaVersion 


  # Align flowcell-level unmapped input bams in parallel
  scatter (fastq in sample_fastqs) {
    String readgroup = fastq[1]

    call FastqToUbam {
      input:
        sample_name = fastq[0],
        readgroup = readgroup,
        fastq1 = fastq[2],
        fastq2 = fastq[3],
        outpref = sample_name + "." + readgroup
    }

    # Map reads to reference
    call SamToFastqAndBwaMem {
      input:
        input_bam = FastqToUbam.unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = sample_name + "." + readgroup + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_alt = ref_alt,
	compression_level = compression_level
     }

    # Merge original uBAM and BWA-aligned BAM 
    call MergeBamAlignment {
      input:
        unmapped_bam = FastqToUbam.unmapped_bam,
        bwa_commandline = bwa_commandline,
        bwa_version = GetBwaVersion.version,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = sample_name + "." + readgroup + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        compression_level = compression_level
    }
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
     input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags {
    input:
      input_bam = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      compression_level = compression_level,
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
  }
  
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }  
  }  
  
  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      compression_level = compression_level  
 }

  scatter (interval_file in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = GatherBamFiles.output_bam,
        input_bam_index = GatherBamFiles.output_bam_index,
        interval_list = interval_file,
        output_filename = base_file_name + ".g.vcf.gz",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_filename = base_file_name + ".g.vcf.gz"
  }

  # Outputs that will be retained when execution is complete  
  output {
    Array[File] unmapped_bam = FastqToUbam.unmapped_bam
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File analysis_ready_bam = GatherBamFiles.output_bam
    File analysis_ready_bam_index = GatherBamFiles.output_bam_index
    File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
    File output_vcf = MergeGVCFs.output_vcf
    File output_vcf_index = MergeGVCFs.output_vcf_index
  } 
}

# TASK DEFINITIONS

# Generate ubam from fastq (per lane ubams)

task FastqToUbam {

  String fastq1
  String fastq2
  String sample_name
  String readgroup
  String outpref 

  command {
    source /data/miniconda/bin/activate gatk4_pipeline && \
    java -Xmx12G -Xms12G -jar /data/cromwellGATK4/picard.jar FastqToSam \
    FASTQ=${fastq1} \
    FASTQ2=${fastq2} \
    OUTPUT=${outpref}.unmapped.bam \
    READ_GROUP_NAME=${readgroup} \
    SAMPLE_NAME=${sample_name} \
    LIBRARY_NAME=${sample_name} \
    PLATFORM=ILLUMINA
  }
  output {
    File unmapped_bam = "${outpref}.unmapped.bam"
  }
}


# Get version of BWA
task GetBwaVersion {

  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    source /data/miniconda/bin/activate gatk4_pipeline && \
    bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  output {
    String version = read_string(stdout())
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File input_bam
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  Int compression_level

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    source /data/miniconda/bin/activate gatk4_pipeline && \
		java -Dsamjdk.compression_level=${compression_level} -Xms30000m -Xmx30000m -jar /data/gatk4_multisample/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=/dev/stdout \
			INTERLEAVE=true \
			NON_PF=true \
    | \
		${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
    | \
		samtools view -1 - > ${output_bam_basename}.bam

  >>>
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  String bwa_commandline
  String bwa_version
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  Int compression_level

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms12000m -Xmx12000m" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${aligned_bam} \
      --UNMAPPED_BAM ${unmapped_bam} \
      --OUTPUT ${output_bam_basename}.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER "unsorted" \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 2000000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  
  Int compression_level

  command {
    set -o pipefail
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms20000m -Xmx20000m" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "queryname" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms4000m -Xmx4000m" \
      SetNmAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  }

  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
#
#  Int compression_level
#
# # Task is assuming query-sorted input so that the Sesourcery and Supplementary reads get marked correctly.
  command {
  source /data/miniconda/bin/activate gatk4_pipeline && \
  /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Dsamjdk.compression_level=2 -Xms20000m" \
      MarkDuplicates \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --METRICS_FILE ${metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
  }

  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}
# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict  
  

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because b37 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    b37_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + b37_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + b37_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + b37_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  

  command { 
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Xms12000m -Dsamjdk.use_async_io_read_samtools=true" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename


  command {
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Xms8000m -Xmx8000m" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index


  command {  
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Xms12000m -Xmx12000m -Dsamjdk.use_async_io_read_samtools=true" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_bam_basename}.bam \
      -L ${sep=" -L " sequence_group_interval} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename

  Int compression_level

  command {
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms12000m -Xmx12000m" \
      GatherBamFiles \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
    }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  } 
}

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination

  String? java_options
  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])


  command <<<
  set -e
    source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/cromwelGATK4/gatk-4.2.0.0/gatk --java-options "-Xmx30000m -Xms30000m ${java_opt} -Dsamjdk.use_async_io_read_samtools=true --native-pair-hmm-threads=10 -pairHMM AVX_LOGLESS_CACHING" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -L ${interval_list} \
      -O ${output_filename} \
      -contamination ${default=0 contamination} \
      -ERC GVCF
  >>>

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_filename


  command <<<
  set -e
  source /data/miniconda/bin/activate gatk4_pipeline && \
    /data/gatk4_multisample/gatk-4.2.0.0/gatk --java-options "-Xmx60000m -Xms60000m"  \
      MergeVcfs \
      --INPUT ${sep=' --INPUT ' input_vcfs} \
      --OUTPUT ${output_filename}
  >>>

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}
