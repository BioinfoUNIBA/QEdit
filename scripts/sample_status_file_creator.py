import sys

try:
        input_file = sys.argv[1]
        sample_group1 = sys.argv[2]
        sample_group2 = sys.argv[3]
except:
        sys.exit('<csv_input_file, sample_group1, sample_group2>')

output_file = sample_group1 + '_vs_' + sample_group2 + '.sif'
with open(output_file,'w') as e:
        samples = []
        with open(input_file,'r') as f:
            header = 'Sample,Group,Type\n'
            e.write(header)
            for line in f:            
                line = line.split(',')
                body_site = line[0]
                srr = line[3]                     
                if (body_site == sample_group1):
                        samples.append((srr, 'GROUPB', body_site))
                elif(body_site == sample_group2):
                        samples.append((srr, 'GROUPA', body_site))
        for sample in sorted(samples, key = lambda x:x[1]):                     
            e.write(','.join((sample[0],sample[1],sample[2])) + '\n')
