# coding=utf-8
# Copyright (c) 2019-2020 Claudio Lo Giudice <clalogiudice@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#################################### REDI OUT TABLE ########################################################
#Region		Position	Reference	Strand	Coverage-q30	MeanQ	BaseCount[A,C,G,T]	   #
#AllSubs	Frequency	gCoverage-q30	gMeanQ	gBaseCount[A,C,G,T]	gAllSubs	gFrequency #
############################################################################################################

###################################GET_DE_events_table#####################################################
#chromosome	position	editing_type	SRR3306830_CTRL		SRR3306831_CTRL	  SRR3306832_CTRL #
#SRR3306833_CTRL	SRR3306834_CTRL	SRR3306835_CTRL   SRR3306836_CTRL	SRR3306823_DIS            #	
#SRR3306824_DIS	    SRR3306825_DIS	SRR3306826_DIS	SRR3306827_DIS	SRR3306828_DIS	etc.              #
#[num_controls/num_disease]	delta_diff	pvalue (Mannwhitney)       				  #			  #
###########################################################################################################



import os, sys, argparse
from scipy import stats
from scipy.stats import wilcoxon, mannwhitneyu, fisher_exact
import numpy as np
import pandas as pd
import math
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument("-c", action = 'store', dest = 'min_coverage', 
		type = int, default=10,  help='Coverage-q30')
parser.add_argument("-cpval", action = 'store', dest = 'pvalue_correction',
                type = int, default = 0, help = '1 --> Bonferroni correction / 2 --> Benjamini hochberg')
parser.add_argument("-input_file", action = 'store', dest = 'samples_informations_file',
		type = str, default= 'empty', help = 'Comma separated file  e.g: Sample,Group,Type \
		SRR1093527,GROUPA,BrainCerebellum... SRR1088437,GROUPB,ArteryTibial... etc')
parser.add_argument("-gene_pos_file", action = 'store', dest = 'gene_pos_file',
		type = str, default= 'empty', help = 'nonsynonymous_table_NONREP derived from Rediportal\
					NOTE: A gene_pos_file is required by -graph or -rsite options. \
An example file can be found at "https://github.com/BioinfoUNIBA/QEdit/blob/master/Example_files/nonsynonymous_table_NONREP_2BS.txt"')
parser.add_argument("-f", action = 'store', dest = 'min_edit_frequency',
		type = float, default=0.1, help='Editing Frequency')
#parser.add_argument("-mts", action = 'store', dest = 'min_sample_testing',
#		type = float, default=50.0, help="min percentage of each sample category")
parser.add_argument("-mtsA", action = 'store', dest = 'groupA_min_sample_testing',
		type = float, default=50.0, help="min percentage of groupA samples")
parser.add_argument("-mtsB", action = 'store', dest = 'groupB_min_sample_testing',
		type = float, default=50.0, help="min percentage of groupB samples")
parser.add_argument("-sig", action = 'store', dest = 'only_significant',
		type = str, default = 'no', help = 'Return only significant editing events')
parser.add_argument("-siglevel", action = 'store', dest = 'statistical_significance',
        type = float, default = 0.05, help = 'cutoff level to reject H0 hypothesis default 0.05')
parser.add_argument("-linear", action = 'store_true', help = 'Enable linear statistical model')
parser.add_argument("-graph", action = 'store_true', help = 'R graph compatible table \
					containing the following columns: \
					Edited_Site | Delta_mean | log_padjstd	| color \
					NOTE: THIS OPTION CAN BE USED ONLY IN COMBINATION with -Gene_pos_file' )
parser.add_argument("-chr_col", action = 'store', dest = 'chr_column',
		type = str, default = 'no', help = 'If set to "yes" a chromosome_position column will be\
					aded to R graph table. \
					NOTE: THIS OPTION IS SPECIFIC FOR -graph & -Gene_pos_file COMBINATION' )
parser.add_argument("-rsite", action = 'store_true', help = 'If set to "yes" all recoding sites will be shown in\
					the output table. \
					NOTE: THIS OPTION ONLY WORKS IN DEFAULT MODE.')

args = parser.parse_args()
min_coverage = args.min_coverage
min_edit_frequency = args.min_edit_frequency
min_sample_testing_A = args.groupA_min_sample_testing
min_sample_testing_B = args.groupB_min_sample_testing
only_significants = args.only_significant
statistical_significance = args.statistical_significance
pvalue_correction = args.pvalue_correction
samples_informations_file = args.samples_informations_file
gene_pos_file = args.gene_pos_file
enable_linear_model = args.linear
enable_graph = args.graph
chr_column = args.chr_column
rsite = args.rsite


if args.graph == True:
	if (((args.graph == True) and (args.gene_pos_file == 'empty')) 
		or ((args.graph == False) and (not args.gene_pos_file == 'empty'))):
			parser.error('gene_pos_file from Rediportal is MISSING' + \
			' or -graph option has not been selected!' + "\n" + \
		'Please type "python get_DE_events.py -h" for more details on usage of this script.')

if args.rsite == True:
	if (((args.rsite == True) and (args.gene_pos_file == 'empty')) 
		or ((args.rsite == False) and (not args.gene_pos_file == 'empty'))):
			parser.error('gene_pos_file from Rediportal is MISSING' + \
			' or -rsite option has not been selected!' + "\n" + \
		'Please type "python get_DE_events.py -h" for more details on usage of this script.')

if args.samples_informations_file == 'empty':
	parser.error('sample_informations_file is MISSING!' + '\n' + \
	'Please type "python get_DE_events.py -h" for more details on usage of this script.')

def call_differential_editing_sites(config_file):
	stability_value = 0.03 #value below which you may use a lower coverage for adding more samples to increase power
	min_disease_people = 5 #min number people supporting higher coverage for whch you may base stability off measurements off of
	min_control_people = 5  #min number control poeple supporting higher coverage for which you may base stability off of
	min_disease_people_5_cov = 10 #min disease number of people of 5 coverage you must have if needing to use unstable 5x coverage
	min_control_people_5_cov = 10 #min control number of people of 5 coverage you must have if needing to use unstable 5x coverage
	editing_file= './temp.csv'
	output_file = './groupA_vs_groupB_linear_model.csv'
	editing_table = pd.read_csv(editing_file,sep='\t')
	config_table = pd.read_csv(config_file,sep=',',skiprows=1,header=None)
	all_people = config_table[0]
	#disease_people = config_table[0][config_table[1] == "DIS"].reset_index(drop = True) #TODO Change do disease!!!
	#control_people = config_table[0][config_table[1] == "CTRL"].reset_index(drop = True) #TODO Change to control!!!

	disease_people = config_table[0][config_table[1] == "GROUPB"].reset_index(drop = True) #TODO Change do disease!!!
	control_people = config_table[0][config_table[1] == "GROUPA"].reset_index(drop = True) #TODO Change to control!!!
	#now get just an editing table and coverage table
	edit_level_table = editing_table[all_people]
	#edit_level_table = editing_table[np.r_[all_people]]

	def get_editing_levels_for_cov_table(i):
	  info = i.astype(str).str.split(pat="\\^")
	  editing_levels = info.apply(lambda x: float('nan') if x[0] == "nan" else x[2])
	  return editing_levels
	cov_table = edit_level_table.apply(get_editing_levels_for_cov_table)
	cov_table = cov_table.apply(lambda x: pd.to_numeric(x)) #TODO check if as.numeric and pandas to_numeric do the same.

	def get_editing_levels(i):
	  info = i.astype(str).str.split(pat="\\^")
	  editing_levels = info.apply(lambda x: float('nan') if x[0] == "nan" else x[0])
	  return editing_levels
	edit_level_table = edit_level_table.apply(get_editing_levels)
	edit_level_table = edit_level_table.apply(lambda x: pd.to_numeric(x)) #TODO check precision on R and python

	#go down line by line and get the prevalence info and mean editing levels based off of stable coverages
	#WARNING I'm using float here, not integer allowing NaN values. Is ok?
	coverage_threshold_used = np.repeat(0.,edit_level_table.shape[0]) #will hold the coverage threshold required for this editing site
	stability_based_on = np.repeat(0.,edit_level_table.shape[0]) #will hold what coverage stability requirements were determined
	stable_mean_disease_editing_level = np.repeat(0.,edit_level_table.shape[0]) #mean autistic editing level using individuals passing coverage threshold
	stable_std_dev_disease_editing_level = np.repeat(0.,edit_level_table.shape[0]) #standard deviation of autistic editing level using individuals passing coverage threshold
	stable_mean_control_editing_level = np.repeat(0.,edit_level_table.shape[0]) #mean control editing level using individuals passing coverage threshold
	stable_std_dev_control_editing_level = np.repeat(0.,edit_level_table.shape[0]) #standard deviation of control editing level using individuals passing coverage threshold
	stable_number_disease_with_at_least_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #number of autistic individuals passing the coverage threshold
	stable_number_disease_nonzero_editing_and_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #number of autistic individuals without non zero editing level and passing coverage threshold
	stable_disease_prevalence = np.repeat(0.,edit_level_table.shape[0]) #proportion autistic individuals with nonzero editing
	stable_number_control_with_at_least_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #same as disease but for control subjects
	stable_number_control_nonzero_editing_and_min_coverage = np.repeat(0.,edit_level_table.shape[0])
	stable_control_prevalence = np.repeat(0.,edit_level_table.shape[0])
	stable_total_number_individuals_nonzero_editing_and_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #total number of disease and control subjects passing the coverage threshold and having nonzero editing level
	stable_mann_whitney_p_value = np.repeat(0.,edit_level_table.shape[0]) #wilcoxon rank sum test p value using individuals passing the coverage threshold
	stable_editing_level_effect_size = np.repeat(0.,edit_level_table.shape[0]) #difference between mean disease and mean control
	stable_frequency_fishers_p_value = np.repeat(0.,edit_level_table.shape[0]) #prevalence p value determined using two-tailed fisher's exact test
	stable_frequency_OR = np.repeat(0.,edit_level_table.shape[0]) #odds ratio of the fisher's exact teest
	stable_prevalence_effect_size = np.repeat(0.,edit_level_table.shape[0]) #difference in editing level prevalences between disease and control subjects
	#WARNING those are np arrays.

	for i in range(0,edit_level_table.shape[0]):
	  print i  #keep track of progress
	  disease_edit_row = edit_level_table.loc[i, disease_people]
	  control_edit_row = edit_level_table.loc[i, control_people]
	  disease_cov_row = cov_table.loc[i, disease_people]
	  control_cov_row = cov_table.loc[i, control_people]
	  #find what coverage we can base stability off of
	  number_disease_20_cov = disease_cov_row[disease_cov_row >= 20].count()
	  number_control_20_cov = control_cov_row[control_cov_row >=20].count()
	  number_disease_15_cov = disease_cov_row[disease_cov_row >= 15].count()
	  number_control_15_cov = control_cov_row[control_cov_row >= 15].count()
	  number_disease_10_cov = disease_cov_row[disease_cov_row >= 10].count()
	  number_control_10_cov = control_cov_row[control_cov_row >= 10].count()
	  number_disease_5_cov = disease_cov_row[disease_cov_row >= 5].count()
	  number_control_5_cov = control_cov_row[control_cov_row >= 5].count()
	  if number_disease_20_cov >= min_disease_people and number_control_20_cov >= min_control_people:
		stability_based_on[i] = 20
	  elif number_disease_15_cov >= min_disease_people and number_control_15_cov >= min_control_people:
		stability_based_on[i] = 15
	  elif number_disease_10_cov >= min_disease_people and number_control_10_cov >= min_control_people:
		stability_based_on[i] = 10
	  elif number_disease_5_cov >= min_disease_people_5_cov and number_control_5_cov >= min_control_people_5_cov:
		stability_based_on[i] = 5
	  else:
		#stability_based_on[i] = -99999 # there's no np.nan integer representation, only float. We use an invalid value.
		stability_based_on[i] = float('nan')

	  #need to deal with cases where there just are not enough disease individuals or control individuals to calculate mean
	  if np.isnan(stability_based_on[i]):

		coverage_threshold_used[i] = 5 #I warn users not to use editing sites that don't have any stability_based_on measurement. We include min coverage of 5 just to get statistical information anyways
		#stable_min_cov=5
		#otherwise we can now try to find the stable_min_cov that'll be used for calculation of all statistics'

	  else:
		current_stability_cov =  stability_based_on[i]
		stability_disease_mean = disease_edit_row[disease_cov_row >= current_stability_cov].mean()
		stability_control_mean = control_edit_row[control_cov_row >= current_stability_cov].mean()
		#print np.arange(5,stability_based_on[i]+1e-4,5)
		for j in np.arange(5,stability_based_on[i]+1e-4,5): #WARNING using 1e-4 allowing to include stop
		  disease_mean = disease_edit_row[disease_cov_row >= j].mean()
		  control_mean = control_edit_row[control_cov_row >= j].mean()
		  if np.absolute(disease_mean-stability_disease_mean) <=stability_value and np.absolute(control_mean-stability_control_mean) <=stability_value :
		    coverage_threshold_used[i] = j
		    break
	  #now let's calculate all our statics based on the stable coverage threshold
	  stable_min_cov = coverage_threshold_used[i]
	  disease_adju_edit_row = disease_edit_row[np.logical_and(np.logical_and((~np.isnan(disease_edit_row)), (~np.isnan(disease_cov_row))), (disease_cov_row >= stable_min_cov))]
	  disease_adju_cov_row = disease_cov_row[np.logical_and((~np.isnan(disease_cov_row)), (disease_cov_row >= stable_min_cov))]
	  control_adju_edit_row = control_edit_row[ np.logical_and(np.logical_and((~np.isnan(control_edit_row)), (~np.isnan(control_cov_row))), (control_cov_row >= stable_min_cov))]
	  control_adju_cov_row = control_cov_row[np.logical_and((~np.isnan(control_cov_row)), (control_cov_row >= stable_min_cov))]
	  stable_mean_disease_editing_level[i] = disease_adju_edit_row.mean()
	  stable_std_dev_disease_editing_level[i] = disease_adju_edit_row.std()
	  stable_mean_control_editing_level[i] = control_adju_edit_row.mean()
	  stable_std_dev_control_editing_level[i] = control_adju_edit_row.std()
	  stable_number_disease_with_at_least_min_coverage[i] = disease_adju_cov_row[disease_adju_cov_row >=stable_min_cov].count()
	  stable_number_disease_nonzero_editing_and_min_coverage[i] = disease_adju_cov_row[ (~np.isnan(disease_adju_cov_row)) & (disease_adju_cov_row >= stable_min_cov) & (disease_adju_edit_row > 0) ].count()
	  stable_disease_prevalence[i] = stable_number_disease_nonzero_editing_and_min_coverage[i]/stable_number_disease_with_at_least_min_coverage[i]
	  stable_number_control_with_at_least_min_coverage[i] = control_adju_cov_row[control_adju_cov_row >=stable_min_cov].count()
	  stable_number_control_nonzero_editing_and_min_coverage[i] = control_adju_cov_row[(~np.isnan(control_adju_cov_row)) & (control_adju_cov_row >= stable_min_cov) & (control_adju_edit_row > 0)].count()
	  stable_control_prevalence[i] = stable_number_control_nonzero_editing_and_min_coverage[i]/stable_number_control_with_at_least_min_coverage[i]
	  stable_total_number_individuals_nonzero_editing_and_min_coverage[i] = (stable_number_disease_nonzero_editing_and_min_coverage[i] + stable_number_control_nonzero_editing_and_min_coverage[i]).sum()
	  if (len(disease_adju_edit_row) >=1) & (len(control_adju_edit_row) >=1):
		if (np.all(disease_adju_edit_row.values == control_adju_edit_row.values)):
		  stable_mann_whitney_p_value[i] = float('nan')
		else:
		  temp, stable_mann_whitney_p_value[i] = mannwhitneyu(disease_adju_edit_row,control_adju_edit_row, alternative='two-sided')
	  else:
		stable_mann_whitney_p_value[i] = float('nan')
	  stable_editing_level_effect_size[i] =  np.absolute(stable_mean_disease_editing_level[i] - stable_mean_control_editing_level[i])
	  fisher_matrix = np.matrix([[stable_number_disease_nonzero_editing_and_min_coverage[i], stable_number_disease_with_at_least_min_coverage[i]-stable_number_disease_nonzero_editing_and_min_coverage[i]], [stable_number_control_nonzero_editing_and_min_coverage[i], stable_number_control_with_at_least_min_coverage[i]-stable_number_control_nonzero_editing_and_min_coverage[i]]])
	  stable_frequency_OR[i], stable_frequency_fishers_p_value[i] = fisher_exact(fisher_matrix)  
	  #print stable_frequency_OR[i]
	  #print stable_frequency_fishers_p_value[i]
	  stable_prevalence_effect_size[i] = np.absolute(stable_disease_prevalence[i] - stable_control_prevalence[i])

	#now put everything back together as a table
	header_info = editing_table[['chromosome','position','editing_type']]
	stats_table = pd.DataFrame(coverage_threshold_used)
	stats_table = stats_table.rename(columns={stats_table.columns[0]: 'coverage_threshold_used'})
	stats_table['stability_based_on'] = pd.DataFrame(stability_based_on)
	stats_table['stable_mean_disease_editing_level'] = pd.DataFrame(stable_mean_disease_editing_level)
	stats_table['stable_std_dev_disease_editing_level'] = pd.DataFrame(stable_std_dev_disease_editing_level)
	stats_table['stable_mean_control_editing_level'] = pd.DataFrame(stable_mean_control_editing_level)
	stats_table['stable_std_dev_control_editing_level'] = pd.DataFrame(stable_std_dev_control_editing_level)
	stats_table['stable_number_disease_with_at_least_min_coverage'] = pd.DataFrame(stable_number_disease_with_at_least_min_coverage)
	stats_table['stable_number_disease_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_number_disease_nonzero_editing_and_min_coverage)
	stats_table['stable_disease_prevalence'] = pd.DataFrame(stable_disease_prevalence)
	stats_table['stable_number_control_with_at_least_min_coverage'] = pd.DataFrame(stable_number_control_with_at_least_min_coverage)
	stats_table['stable_number_control_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_number_control_nonzero_editing_and_min_coverage)
	stats_table['stable_control_prevalence'] = pd.DataFrame(stable_control_prevalence)
	stats_table['stable_total_number_individuals_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_total_number_individuals_nonzero_editing_and_min_coverage)
	stats_table['stable_mann_whitney_p_value'] = pd.DataFrame(stable_mann_whitney_p_value)
	stats_table['stable_editing_level_effect_size'] = pd.DataFrame(stable_editing_level_effect_size)
	stats_table['stable_frequency_fishers_p_value'] = pd.DataFrame(stable_frequency_fishers_p_value)
	stats_table['stable_frequency_OR'] = pd.DataFrame(stable_frequency_OR)
	stats_table['stable_prevalence_effect_size'] = pd.DataFrame(stable_prevalence_effect_size)

	full_table = pd.concat([header_info, stats_table, editing_table[all_people]], axis=1)
	if os.path.exists('temp.csv'):
		os.remove('temp.csv')
	#write the full_table to output
	full_table.to_csv(output_file, sep='\t', index=False)

	print "job completed\n"

def Set_Chr_Nr(Chr):
    """ Sort by chromosome """
    if Chr: 
	New = Chr.lstrip('chr').split('_')[0]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New

def Sample_percentage(row):
	"""Percentage of samples from each type"""
	percentage = (len(filter(lambda x: x!= '-', row))/float(len(row)))*100
	return round(percentage)

	
def Sample_count(row):
	"""Number of samples from each type"""
	count = len(filter(lambda x: x!= '-', row))
	return count 

def get_bh(pvalue,siglevel): # rows list as input
	"""B-H correction """
	pvalue = sorted(pvalue, key = lambda x: x[-1])
	x=1
	for i in pvalue:
		nf=float(i[-1])*len(pvalue)
		fdr=round(nf/x,6)
		if fdr<=siglevel:
			i.append(fdr)
			i.append('yes')
		else:
			i.append(fdr) 
			i.append('no')
		x+=1
	return pvalue#,y,p


	
def get_b(pvalue,siglevel): 
	"""Bonferroni correction"""
	pvalue = sorted(pvalue, key = lambda x: x[-1])
	#y=0
	pp=1.0
	for i in pvalue:
		p=float(i[-1])*len(pvalue)
		if p<=siglevel:
			i.append(p)
			i.append('yes')
			if p<pp: pp=p
		else: 
			i.append(p)
			i.append('no')
	return pvalue#,y,pp

def get_ed_level(l,m):
	"""Get editing level (under_edited, overedited, nonsignificant)"""
	if float(m) < 1:
		return "nonsignificant"
	else:
		if np.sign(float(l)) == -1.0:
			return "under_edited"
		else:
			return "over-edited"

def call_sites(gene_pos_file):
	dic_aa_changes = {} # key(chrom,pos) values [gene, aa_change]
	with open(gene_pos_file, 'r') as g:
		for i in g:
			l = (i.strip()).split('\t')
			chrmsm, pstn, gene,  = l[0], int(l[1]), l[10]
			aa_change = ((l[12].split(',')[0]).split(':')[-1]).split('.')[1]
			dic_aa_changes.setdefault((chrmsm, pstn), [gene, aa_change])
	return dic_aa_changes	

def only_sig(row):
	"""Returns only significant events"""
	if row[-1] == 'yes':
		return row

def tuple_replace(i):
	if type(i) == tuple:
		return i[0]
	else:
		return i

def tuple_replace_bis(k):
        if type(k) == tuple:
                return k[1]
        else:
             	return k

def type_ed(x):
    if str(x).find('^AG') != -1 or str(x).find('^TC') != -1: 
		ed_type = x[1].split('^')[1]
		return ed_type

def remove_underscore(lis):
	lis = lis[:lis.index('_')]
	return lis

sample_informations = {}
with open(samples_informations_file, 'r') as f:
    for line in f:
        if line.startswith('SRR'):
            line = map(str.strip, line.split(','))
            sample_informations.setdefault(line[0], [line[1], line[2]]) #line[1]

cwd = filter(os.path.isdir, os.listdir(os.getcwd()))
all_available_sites = []
sample_edited_sites = {}
for directory in cwd:
    if directory.startswith('SRR'):
		path = list(os.walk(directory + '/editing/'))
		table = path[1][0] + '/' + path[1][-1][-1] 
		with open(table,'r') as a:
			for line in a:
				if line.startswith('chr'):
					s = map(str.strip, line.split("\t"))
					if s[7].split(' ')[0] in ['AG','TC']:
						site, freq, coverage, ed_type = s[0] + "_" + s[1], s[8], s[4], s[7]
						if site not in all_available_sites: 
							all_available_sites.append(site)
						if (int(coverage) >= min_coverage) and (float(freq) >= min_edit_frequency):
							if ed_type == 'AG': 
								gnum_cov = eval(s[6])[2]
								freq_ed_type_gnum_cov = '%s^%s^%s' %(freq, ed_type, gnum_cov)								
								sample_edited_sites.setdefault((directory, site), []).append((freq, freq_ed_type_gnum_cov))
							else:
								cnum_cov = eval(s[6])[1]
								freq_ed_type_cnum_cov = '%s^%s^%s' %(freq, ed_type, cnum_cov)	
								sample_edited_sites.setdefault((directory, site), []).append((freq, freq_ed_type_cnum_cov))							


table_columns = map(lambda x: x + '_' + '_'.join(sample_informations[x]), sorted(sample_informations.keys()))

disease = [i.replace('_GROUPB','') for i in table_columns if i.upper().find('GROUPB') != -1]
controls = [i.replace('_GROUPA','') for i in table_columns if i.upper().find('GROUPA') != -1]

if enable_linear_model:
	outtable=''
	header = ['chromosome', 'position', 'editing_type'] + map(remove_underscore, controls) + map(remove_underscore, disease)
	outtable += '\t'.join(header) + '\n'
	for chrom in sorted(all_available_sites, key = lambda x: Set_Chr_Nr(x)):
		row = [chrom]
		for col in header[2:]:
			row.append(sample_edited_sites.get((col.split('_')[0],chrom), ['-'])[0])
		ctrls = zip(*(zip(controls,row[1:])))[1]
		dss = zip(*(zip(disease,row[len(ctrls)+1:])))[1]
		ctrls_freq = map(tuple_replace, ctrls)
		dss_freq = map(tuple_replace, dss)

		row_b = map(tuple_replace_bis, row)
		row_b = row_b[0].split('_') + row_b[2:]
		if any(filter(None,map(type_ed,row))):
			ed = filter(None,map(type_ed,row))[0]
			row_b.insert(2,ed)
		else:
			row_b.insert(2,'-')

		outtable += '\t'.join(map(str,row_b)).replace('-','NA') + '\n'


	with open('temp.csv','w') as t:
		t.write(outtable)
		t.close()

	call_differential_editing_sites(samples_informations_file) 
	

else:
	header = ['chromosome', 'position', 'editing_type'] + controls + disease + ['[groupA_samples/groupB_samples]'] \
	+ ['delta_diff'] + ['pvalue (Mannwhitney)'] 

	if pvalue_correction == 1:
		header += ['pvalue Bonferroni corrected', 'significant']
	if pvalue_correction == 2:
		header += ['pvalue BH corrected', 'significant']	

	processed_rows = []
	unprocessed_rows = []
	for chrom in sorted(all_available_sites, key = lambda x: Set_Chr_Nr(x)):
		row = [chrom]
		for col in header[3:header.index('[groupA_samples/groupB_samples]')]:
			row.append(sample_edited_sites.get((col.split('_')[0],chrom), ['-'])[0])
		ctrls = zip(*(zip(controls,row[1:])))[1]
		dss = zip(*(zip(disease,row[len(ctrls)+1:])))[1] 	
		ctrls_freq = map(tuple_replace, ctrls)
		dss_freq = map(tuple_replace, dss)
		row.append(str([Sample_count(ctrls), Sample_count(dss)]))
		if (Sample_percentage(ctrls) >= min_sample_testing_A) and (Sample_percentage(dss) >= min_sample_testing_B):
			#row.append([Sample_count(ctrls), Sample_percentage(ctrls), min_sample_testing_A, Sample_count(dss), Sample_percentage(dss), min_sample_testing_B])
			ctrls_mean = sum(map(float, filter(lambda x: x!= '-', ctrls_freq)))/len(filter(lambda x: x!= '-', ctrls_freq))
			dss_mean = sum(map(float, filter(lambda x: x!= '-', dss_freq)))/len(filter(lambda x : x!= '-', dss_freq))
			delta_diff = ctrls_mean - dss_mean  #abs(ctrls_mean - dss_mean)
			pvalue=stats.mannwhitneyu(ctrls_freq, dss_freq, alternative='two-sided')
			row.append(str(round(delta_diff, 3))) #delta_diff
			row.append(str(round(pvalue[1], 6))) #MW pvalue
			processed_rows.append(row)
		else:
			if pvalue_correction == 0:
				row += ['-', '-']
			else:
				row += ['-', '-', '-', '-']
			unprocessed_rows.append(row)

	if pvalue_correction == 1:
		processed_rows = get_b(processed_rows, statistical_significance) #default 0.05
	if pvalue_correction == 2:
		processed_rows = get_bh(processed_rows, statistical_significance) #default 0.05

	total_rows = processed_rows + unprocessed_rows	

	if enable_graph:
		dic_aa_changes = call_sites(gene_pos_file)
		header = ['Site', 'Delta', 'Mannwhitney pval','Benjamini Hochberg corrected pvalue', 'status']
		if chr_column == 'yes' : header.insert(0,'Chrom_pos')
		
		out = filter(None, map(only_sig, get_bh(processed_rows, statistical_significance)))
		#out = get_bh(processed_rows, statistical_significance)
		print '\t'.join(header)

		for l in sorted(out, key = lambda x: Set_Chr_Nr(x[0])):
			l = map(tuple_replace_bis, l)
			chr_pos = (l[0].split('_')[0], int(l[0].split('_')[1])) 
			site = dic_aa_changes.get(chr_pos, '-')
			if site != '-': site = '_'.join(site)
			delta = l[-4]
			pmw = l[-3]
			padjstd = l[-2]
			adjstd_plog10 = -(round(math.log(float(l[-2])), 3)) 	
			ed_level = get_ed_level(delta, adjstd_plog10)
			if header[0] == 'Chrom_pos':
				chr_pos = '_'.join(map(str,list(chr_pos)))
				l =  [chr_pos, site, delta, pmw, padjstd, ed_level]
			else:
				l = [site, delta, adjstd_plog10, ed_level] 
			print '\t'.join(map(str,l))
	
	else:
		if (pvalue_correction != 0 and only_significants == 'no'):
			out = total_rows
		elif (pvalue_correction != 0 and only_significants == 'yes'):
			out = filter(None, map(only_sig, processed_rows))
		else:
			out = total_rows #without correction

		if rsite:
			dic_aa_changes = call_sites(gene_pos_file)
			header.insert(2,'Recoding_site')

			print '\t'.join(header)
			for l in sorted(out, key = lambda x: Set_Chr_Nr(x[0])):
				chr_pos = (l[0].split('_')[0], int(l[0].split('_')[1])) 
				site = dic_aa_changes.get(chr_pos, '-')
				l = l[0].split('_') + l[1:] 
				l.insert(2, '_'.join(site))
				if any(filter(None,map(type_ed,l))):
					ed = filter(None,map(type_ed,l))[0]
					l.insert(2,ed)
				else:
					l.insert(2,'-')
				l = map(tuple_replace_bis, l)
				print '\t'.join(map(str,l))
		else:
			print '\t'.join(header)
			for l in sorted(out, key = lambda x: Set_Chr_Nr(x[0])):
				chr_pos = (l[0].split('_')[0], int(l[0].split('_')[1])) 
				l = l[0].split('_') + l[1:] 
				if any(filter(None,map(type_ed,l))):
					ed = filter(None,map(type_ed,l))[0]
					l.insert(2,ed)
				else:
					l.insert(2,'-')
				l = map(tuple_replace_bis, l)
				print '\t'.join(map(str,l))
