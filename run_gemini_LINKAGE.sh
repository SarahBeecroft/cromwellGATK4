#!/bin/bash

DATABASE=wes_hg38.database
PATIENT1=D18_0398 
PATIENT2=D18_0389
OUTPUT_FILE=test_linkage.txt

#Example use: ./run_gemini_linkage.sh
#This script takes a comma seperated list of chromosomal locations (i.e such as those generated from linkage analysis), and outputs variants in those regions, which are shared by user-specified patients.
for region in chr19:55439166-57646570, chr19:46734255-52115645
do
gemini query --region $region -q \
"select gene, qual, chrom, start, end, ref, alt, codon_change, aa_change, aa_length, impact, impact_severity, transcript, HGVSc, HGVSp, \
        exon,\
        intron, \
        type, \
        cDNA_position, \
        CADD_PHRED, \
        polyphen_pred, \
        sift_pred, \
        Loftool, \
        SpliceAI_pred_DP_AG, \
        SpliceAI_pred_DP_AL, \
        SpliceAI_pred_DP_DG, \
        SpliceAI_pred_DP_DL, \
        SpliceAI_pred_DS_AG, \
        SpliceAI_pred_DS_AL, \
        SpliceAI_pred_DS_DG, \
        SpliceAI_pred_DS_DL, \
        SpliceAI_pred_SYMBOL, \
        gnomADg, \
        gnomADg_AF, \
        gnomADg_AN, \
        gnomADg_AC, \
        gnomADg_AF_popmax, gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, \
        gnomADg_nhomalt, \
        ensembl_gene_id, \
        transcript, \
        is_exonic, \
        is_coding, \
        is_lof, \
        is_splicing, \
        is_canonical, \
        biotype, \
        ensp, \
        swissprot, \
        domains, \
        UNIPARC, \
        UNIPROT_ISOFORM, \
        GENE_PHENO, \
        PHENO, \
        CLIN_SIG, \
        Clinvar, \
        gt_depths.$PATIENT1, gt_ref_depths.$PATIENT1, gt_alt_depths.$PATIENT1, \
        gt_depths.$PATIENT2, gt_ref_depths.$PATIENT2, gt_alt_depths.$PATIENT2 \
        from variants where qual >=100 and (MAX_AF_POPS >= 0.01) and (gnomADg_AF_popmax >= 0.01)" \
       --show-samples \
       --sample-delim ";" \
       --gt-filter "(gt_types.$PATIENT1 != HOM_REF) and (gt_depths.$PATIENT1 >= 10) and \
        (gt_types.$PATIENT2 != HOM_REF) and (gt_depths.$PATIENT2 >= 10) and \
        (gt_types).(*).(==HOM_ALT).(count <=5) and (gt_types).(*).(==HET).(count <=5)" \
        --header \
        $DATABASE >> $OUTPUT_FILE
done
