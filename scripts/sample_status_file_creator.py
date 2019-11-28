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
