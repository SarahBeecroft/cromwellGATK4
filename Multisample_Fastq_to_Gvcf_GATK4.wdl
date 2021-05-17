# WORKFLOW DEFINITION 

import "ruddle_fastq_to_gvcf_single_sample_gatk4.wdl" as single_wf


workflow Multisample_Fastq_to_Gvcf_GATK4 {

  File inputSamplesFile
  Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)

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


  # Align flowcell-level unmapped input bams in parallel
  
  call make_uniq_samples_file {
    input:
      inputSamplesFile = inputSamplesFile,
   }

  scatter (sample in make_uniq_samples_file.uniq_samples) {

   call make_fastq_file {
    input:
      inputSamplesFile = inputSamplesFile,
      sample_name = sample
   }


   call single_wf.Fastq_to_Gvcf_GATK4 { 
      input: 

        sample_fastq_file = make_fastq_file.sample_fastq_file,
        sample_name = sample,

        
        unmapped_bam_suffix = unmapped_bam_suffix,

        ref_name = ref_name,
       
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_alt = ref_alt,

        bwa_commandline = bwa_commandline,
        compression_level = compression_level,
        
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,

        scattered_calling_intervals_list = scattered_calling_intervals_list
    }

  }
  # Outputs that will be retained when execution is complete  
  output {
    Array[Array[File]] unmapped_bam = Fastq_to_Gvcf_GATK4.unmapped_bam
    Array[File] duplication_metrics = Fastq_to_Gvcf_GATK4.duplication_metrics
    Array[File] bqsr_report = Fastq_to_Gvcf_GATK4.bqsr_report
    Array[File] analysis_ready_bam = Fastq_to_Gvcf_GATK4.analysis_ready_bam
    Array[File] analysis_ready_bam_index = Fastq_to_Gvcf_GATK4.analysis_ready_bam_index
    Array[File] analysis_ready_bam_md5 = Fastq_to_Gvcf_GATK4.analysis_ready_bam_md5
    Array[File] output_vcf = Fastq_to_Gvcf_GATK4.output_vcf
    Array[File] output_vcf_index = Fastq_to_Gvcf_GATK4.output_vcf_index
    File uniq_samples = make_uniq_samples_file.uniq_samples_file
  } 


}

# TASK DEFINITIONS

task make_uniq_samples_file {
  File inputSamplesFile

  command {
    cat ${inputSamplesFile} | cut -f1 | sort | uniq > uniq_samples.list
  }

  output {
    File uniq_samples_file = "uniq_samples.list"
    Array[String] uniq_samples = read_lines("uniq_samples.list")
  }
}


task make_fastq_file {
  File inputSamplesFile
  String sample_name

  command {
    cat ${inputSamplesFile} | grep "^${sample_name}\s" > ${sample_name}.fastqs.txt
  }

  output {
    File sample_fastq_file = "${sample_name}.fastqs.txt"
  }
}
