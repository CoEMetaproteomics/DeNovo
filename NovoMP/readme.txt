code by Feng Xian (fengx90) 

## file descriptions


1) novoMP_Mlulti-tier_Filtering.R

## use output from BPS-Novor as input followed by multi-tier filtering: de novo score, charge state, peptide length, mass accuracy, retention time shift, ion mobility shifts

2) Convert_ook0_ccs.py, BPSNovor_output_to_IM2Deep.py, and Fragpipe_PSM_to_IM2Deep_for_Calibration_IncludingRT.py

## Convert_ook0_ccs.py is used to convert measured ion mobility (1/k0 values) to CCS values. 

## BPSNovor_output_to_IM2Deep.py: this file is used to convert BPS-Novor output format for CCS prediction in IM2Deep (e.g., format the peptide sequences with detected modifications). 

## Fragpipe_PSM_to_IM2Deep_for_Calibration_IncludingRT.py: this file is used to convert database search output (e.g. psm.tsv from FragPipe) for retention time calibration in DeepLC and CCS calibration in IM2Deep.

3) diamond_blastp.txt

## the file contains the codes for BlastP analysis using Diamond package.

4) Check_and_Format_tsv_output.py
## this file is used to check the integrity of BlastP analysis output.

5) Filter_BlastP_output_with_LCA_file.py
## this file is used to filter BlastP output based on the LCA output, in order to keep most relevant protein candidates that matched to de novo peptides.

6) Get_NCBI_sequences_fasta.py
## this file is used to automatically download protein sequences of mapped de novo peptides to form a novoMP-DB. 
