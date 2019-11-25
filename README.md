<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />  
  </head>
  <body>
<div>
<h1>QEdit: RNA editing quantification in deep transcriptome data</h1>
</div> 
<p align="justify">
RNA editing is a co-/post-transcriptional mechanism involving the insertion/deletion or substitution of specific bases in precise RNA localizations. Substitutional RNA editing, mostly by the adenosine to inosine (A-to-I) deamination, is prominent in mammals.
It has profound functional consequences and its deregulation has been linked to a variety of human diseases including neurological and neurodegenerative disorders or cancer.
RNA editing can be profiled in deep transcriptome data but its detection is not trivial. Indeed the number of A-to-I candidates as well as the accuracy of predictions depends on the adopted computational strategy. Individual steps or specific software used can dramatically alter the quality of results.
Quantifying RNA editing is relevant to compare independent samples and study its potential role in different experimental conditions.
Here we provide simple scripts to calculate different metrics for RNA editing quantification.
Although a variety of programs to identify RNA editing candidates in RNAseq data has been released, we will profile inosinomes using <a href="https://github.com/BioinfoUNIBA/REDItools">REDItools</a> and <a href="http://srv00.recas.ba.infn.it/atlas/index.html">REDIportal</a>.
</p>
<p align="justify">The bioinformatics pipeline used to investigate RNA editing in RNAseq datasets is depicted in the following Figure. Below each step, we report main bioinformatics tools used, whose commad lines are included in this repository.
</p>  
<img src="fig1_qedit.jpg">
<h3>RNA editing detection</h3>
<h4>De novo approach</h4>
<p align="justify">
RNA editing candidates can be detected using REDItools. There are two current versions: 1) <a href="https://github.com/BioinfoUNIBA/REDItools">REDItools</a> 1.3 or 2) <a href="https://github.com/BioinfoUNIBA/REDItools2">REDItools 2.0</a>.
REDItools2 is a faster re-implementation of REDItools1 for HPC clusters. Its serial version is about ten times faster than REDItools1.
The complete workflow for detecting de novo RNA editing events with REDItools is described <a href="https://github.com/BioinfoUNIBA/REDItools#Nature%20Protocol%20scripts">here</a>.
</p>
</body>
</html>
