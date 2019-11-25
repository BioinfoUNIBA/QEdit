<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />  
  </head>
  <body>
<h1>get_DE_events.py</h1>
<h5>This scripts and its related files are part of the supplemental material for the paper<br>
  "Quantifying RNA editing in deep transcriptome datasets"</h5>
<p class-text="justify">
For control case studies by launching the get_DE_events.py script the user can filter REDItoolDnaRna.py outputs according to the following criteria:
<ul>
<li>RNAseq coverage per position (default <b>10 reads</b>)</li>
<li>Minimum editing frequency per position (default <b>10%</b>)</li>
For each editing candidate, the script applies the Mann–Whitney test to check the significance between the two groups, 
A and B. By default the test is carried out only if the number of editing events per position is equal to 50% of the samples per group. This treshold can be manually modified for both groups by playing with the -mtsA and -mtsB options respectively.
Optionally, p-values can be corrected using Benjamini–Hochberg or Bonferroni tests. 
</ul>
<p>Usage:</p> 
<pre>
usage: get_DE_events.py [-h] [-c MIN_COVERAGE] [-cpval PVALUE_CORRECTION]
                        [-input_file SAMPLES_INFORMATIONS_FILE]
                        [-gene_pos_file GENE_POS_FILE] [-f MIN_EDIT_FREQUENCY]
                        [-mtsA GROUPA_MIN_SAMPLE_TESTING]
                        [-mtsB GROUPB_MIN_SAMPLE_TESTING]
                        [-sig ONLY_SIGNIFICANT]
                        [-siglevel STATISTICAL_SIGNIFICANCE] [-linear]
                        [-graph] [-chr_col CHR_COLUMN] [-rsite RSITE]

optional arguments:
  -h, --help            show this help message and exit
  -c MIN_COVERAGE       Coverage-q30
  -cpval PVALUE_CORRECTION 1 --> Bonferroni correction / 2 --> Benjamini hochberg
  -input_file SAMPLES_INFORMATIONS_FILE Comma separated file e.g: Sample,Group,Type (e.g SRR1093527,GROUPA,BrainCerebellum..., SRR1088437,GROUPB,ArteryTibial... etc
  -gene_pos_file GENE_POS_FILE nonsynonymous_table_NONREP derived from Rediportal NOTE: THIS OPTION CAN BE USED ONLY IN COMBINATION with -graph
  -f MIN_EDIT_FREQUENCY
                        Editing Frequency
  -mtsA GROUPA_MIN_SAMPLE_TESTING
                        min percentage of groupA samples
  -mtsB GROUPB_MIN_SAMPLE_TESTING
                        min percentage of groupB samples
  -sig ONLY_SIGNIFICANT
                        Return only significant editing events
  -siglevel STATISTICAL_SIGNIFICANCE
                        cutoff level to reject H0 hypothesis default 0.05
  -linear               Enable linear statistical model
  -graph                R graph compatible table containing the following
                        columns: Site|Delta|Mannwhitney|pval|Benjamini Hochberg corrected pvalue|status
                        NOTE: THIS OPTION CAN BE USED ONLY IN
                        COMBINATION with -Gene_pos_file
  -chr_col CHR_COLUMN   If set to "yes" a chromosome_position column will be
                        aded to R graph table. NOTE: THIS OPTION IS SPECIFIC
                        FOR -graph & -Gene_pos_file COMBINATION
  -rsite RSITE          If set to "yes" all recoding sites will be shown in
                        the output table. NOTE: THIS OPTION ONLY WORKS IN
                        LINEAR AND DEFAULT MODE.
                                                                                      
<b>e.g.</b> python ../REDItools/accessory/get_DE_events.py -cpval 2 -input_file  sample_information.csv -sig yes
<p class-text="justify">The script will filter REDItoolDnaRna.py outputs for each sample contained in the 
SAMPLES_INFORMATIONS_FILE returning only significant editing events (pval <= 0.05)
in accordance with Benjamini hochberg correction.</p>

</pre>
</body>
</html> 
