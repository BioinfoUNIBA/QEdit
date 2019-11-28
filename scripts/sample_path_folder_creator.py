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
