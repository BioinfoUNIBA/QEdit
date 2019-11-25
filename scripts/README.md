
Skip to content
Pull requests
Issues
Marketplace
Explore
@claudiologiudice
Learn Git and GitHub without any code!

Using the Hello World guide, you’ll start a branch, write comments, and open a pull request.

1
3

    3

BioinfoUNIBA/REDItools
Code
Issues 0
Pull requests 1
Actions
Projects 0
Wiki
Security
Insights
REDItools/accessory/

1

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"

2

  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

3

<html xmlns="http://www.w3.org/1999/xhtml">

4

  <head>

5

    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />  

6

  </head>

7

  <body>

8

<h1>get_DE_events.py</h1>

9

<h5>This scripts and its related files are part of the supplemental material for the paper<br>

10

  "Investigating RNA editing in deep transcriptome datasets with REDItools and REDIportal"</h5>

11

<p class-text="justify">

12

For control case studies by launching the get_DE_events.py script the user can filter REDItoolDnaRna.py outputs according to the following criteria:

13

<ul>

14

<li>RNAseq coverage per position (default <b>10 reads</b>)</li>

15

<li>Minimum editing frequency per position (default <b>10%</b>)</li>

16

For each editing candidate, the script applies the Mann–Whitney test to check the significance between the two conditions, 

17

control and HD. By default the test is carried out only if the number of editing events per position is equal to 50% of the samples per group. 

18

Optionally, p-values can be corrected using Benjamini–Hochberg or Bonferroni tests. 

19

</ul>

20

<p>Usage:</p> 

21

<pre>

22

get_DE_events.py [-h] [-c MIN_COVERAGE] [-cpval PVALUE_CORRECTION]

23

                        [-input_file SAMPLES_INFORMATIONS_FILE]

24

                        [-f MIN_EDIT_FREQUENCY] [-mts MIN_SAMPLE_TESTING]

25

                        [-sig ONLY_SIGNIFICANT] [-linear]

26

  

27

optional arguments:

28

  -h, --help                             show this help message and exit

29

  -c MIN_COVERAGE                        Coverage-q30

30

  -cpval PVALUE_CORRECTION 1 -->         Bonferroni correction / 2 --> Benjamini hochberg

31

  -input_file SAMPLES_INFORMATIONS_FILE  Comma separated file e.g: <b>Sample,Status</b>

32

  -f MIN_EDIT_FREQUENCY                  Editing Frequency

33

  -mts MIN_SAMPLE_TESTING                min percentage of each sample category

34

  -sig ONLY_SIGNIFICANT                  Return only significant editing events 

35

                                         (if -cpval flag is activated)

36

  -linear                                Calculate differential RNA editing according to Tran et al. (2019)

37

                                                                                        

38

<b>e.g.</b> python ../REDItools/accessory/get_DE_events.py -cpval 2 -input_file  sample_information.csv -sig yes

39

<p class-text="justify">The script will filter REDItoolDnaRna.py outputs for each sample contained in the 

40

SAMPLES_INFORMATIONS_FILE returning only significant editing events (pval <= 0.05)

41

in accordance with Benjamini hochberg correction.</p>

42

​

43

</pre>

44

</body>

45

</html> 

46

​

@claudiologiudice
Commit changes
Commit summary
Optional extended description
Commit directly to the master branch.
Create a new branch for this commit and start a pull request. Learn more about pull requests.

    © 2019 GitHub, Inc.
    Terms
    Privacy
    Security
    Status
    Help

    Contact GitHub
    Pricing
    API
    Training
    Blog
    About

