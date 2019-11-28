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
This script compares REDItools output table arising from multiple samples and returns dysregulated RNA editing at recoding events by means of the Mann-Whitney U-test described in Silvestris et al. (2019) or the statistical pipeline proposed by Tran et al. (2019). REDItoools output table are pre-filtered according to these main following criteria.  
<ul>
<li>RNAseq coverage per position (default <b>10 reads</b>)</li>
<li>Minimum editing frequency per position (default <b>10%</b>)</li>
</ul>

<p class-text="justify">For each editing candidate, the script applies the MannWhitney test to check the significance between the two groups, A and B. <br> By default the test is carried out only if the number of editing events per position is equal to 50% of the samples per group. <br> This treshold can be manually modified (for both groups) by playing with the -mtsA and -mtsB options respectively. <br> Returned p-values can be corrected using Benjamini Hochberg or Bonferroni tests.</p> 
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
  <b>-h</b>, --help            show this help message and exit
  
  <b>-c</b> MIN_COVERAGE       Coverage-q30
  
  <b>-cpval</b> PVALUE_CORRECTION 1 --> Bonferroni correction / 2 --> Benjamini Hochberg
  
  <b>-input_file</b> SAMPLES_INFORMATIONS_FILE (.sif) Comma separated file e.g: Sample,Group,Type 
  (e.g SRR1093527,GROUPA,BrainCerebellum..., SRR1088437,GROUPB,ArteryTibial... etc)
  An example file is provided<a href="https://raw.githubusercontent.com/BioinfoUNIBA/QEdit/master/Example_files/csv_input_file"> here</a>
  
  <b>-gene_pos_file</b> GENE_POS_FILE nonsynonymous_table_NONREP derived from Rediportal 
  NOTE: A gene_pos file is required by -graph or -rsite.
  An example file can be found here <a href="https://raw.githubusercontent.com/BioinfoUNIBA/QEdit/master/Example_files/nonsynonymous_table_NONREP_2BS.txt">here</a>.
  
  <b>-f</b> MIN_EDIT_FREQUENCY Editing Frequency
  
  <b>-mtsA</b> GROUPA_MIN_SAMPLE_TESTING min percentage of groupA samples                      
  
  <b>-mtsB</b> GROUPB_MIN_SAMPLE_TESTING min percentage of groupB samples                      
  
  <b>-sig ONLY_SIGNIFICANT</b> Return only statistically significant editing events
  
  <b>-siglevel STATISTICAL_SIGNIFICANCE</b> cutoff level to reject H0 hypothesis default 0.05
  
  <b>-linear</b> Enable linear statistical model (Tran et al., 2019).
  
  <b>-graph</b> R graph compatible table containing the following columns: 
  Site|Delta|Mannwhitney|pval|Benjamini Hochberg corrected pvalue|status
  NOTE: THIS OPTION CAN BE USED ONLY IN COMBINATION with -Gene_pos_file
  
  <b>-chr_col CHR_COLUMN</b> If set to "yes" a chromosome_position column will be added to R graph table. 
  NOTE: THIS OPTION IS SPECIFIC FOR -graph & -Gene_pos_file COMBINATION
  
 <b> -rsite RSITE If set to "yes"</b> all recoding sites will be shown in the output table. 
 NOTE: THIS OPTION ONLY WORKS IN DEFAULT MODE.
                                                                                      
<b>e.g.</b> python ../REDItools/accessory/get_DE_events.py -cpval 2 -input_file  sample_information.csv -sig yes
<p class-text="justify">The script will filter REDItoolDnaRna.py outputs for each sample contained in the
SAMPLES_INFORMATIONS_FILE returning only significant editing events (pval <= 0.05)
in accordance with Benjamini Hochberg correction.</p>

</pre>
<h1>Accessory files</h1>
<ul>
  <li>sample_status_file_creator.py</li>
  This script generates a sample_information.csv file (.sif) compatible with get_DE_events.py. It requires:
  <ul>
    <li> A csv sample file containing the main informations about each sample to be used in the experiment. 
    An example of this file is included.</li>
    <li> Samples group1 (e.g. ArteryTibial) </li>
    <li> Samples group2 (e.g BrainCerebellum) </li>
    
    >python sample_status_file_creator.py csv_input_file, sample_group1, sample_group2
   </ul>
</ul>


</body>
</html> 

