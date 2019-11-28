import os, random, string, sys

try:
        sample_status_file = sys.argv[1]
except:
        sys.exit('<sample_status_file>')
 
def randomStringDigits(stringLenght=8):
        return ''.join(random.choice(string.digits) for i in range(stringLenght))

with open(sample_status_file, 'r') as e:
        e.readline()
        for line in e:
                line = line.split(',')
                srr_dir = line[0].strip()
                dna_rna_folder = srr_dir + '/editing/DnaRna_' + randomStringDigits()
                if not os.path.exists(dna_rna_folder):
                         os.makedirs(dna_rna_folder)
                os.system('cp ../tables/%s %s' %(srr_dir, dna_rna_folder))  #<------- LINE TO MODIFY
#If REDItools output tables (e.g SRR...) are stored in a different folder, please modify the previous line accordingly.
