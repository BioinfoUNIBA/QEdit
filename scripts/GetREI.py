
import sys, os, time, math, random, getopt, operator, string, errno
from collections import Counter
try: import pysam
except: sys.exit('Pysam module not found.')


pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python GetREI.py [options]
Options:
-i		Table file [mandatory]
-r		RECODING positions [mandatory]
-p		Min coverage of recoding sites [0]
-o		Output file base name (outfile_%s)
-h		Print this help

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:r:p:h")
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

cRec=0 #Min coverage of recoding sites
tabfile=''
fastafile=''
recodingfile=''
basename='outfile_%s'%(pid)

for o, a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-i": tabfile=a
	elif o == "-o": basename=a
	elif o == "-r": recodingfile=a
	elif o == "-p": cRec=int(a)
	else:
		assert False, "Unhandled Option"

ds={'0':'-','1':'+','2':'+','+':'1','-':'0'}




def eRes(res,lf):
	x=0
	for i in res:
		for j in i:
			if j!='ND':x+=1
	if x==lf: return 1
	return 0 

################################################

if not os.path.exists(tabfile):
	sys.exit('Table file %s not found.' %(tabfile))

table=pysam.Tabixfile(tabfile)
allchrs=table.contigs

if not os.path.exists(recodingfile):
	usage()
	sys.exit('Recoding file %s not found.' %(recodingfile))

outfile2=basename+'.rec'

sys.stderr.write('Reading Recoding...\n')
AG=[0,0]
o=open(outfile2,'w')
head='\t'.join(['Chr','Start/End','Cov','MisMatch','EdFreq','Gene','AAchange','AAchangePos'])
o.write('#'+head+'\n')
f=open(recodingfile)
for i in f:
	l=(i.strip()).split('\t')
	chr=l[0]
	cc=(int(l[1])-1,int(l[1]))
	gene=l[10]
	la=((l[12].split(',')[0]).split(':')[-1]).split('.')[1]
	aa=la[0]+la[-1]
	if chr not in allchrs: continue
	pos=[row.split('\t') for row in table.fetch(chr, cc[0], cc[1])]
	if len(pos)==0:
		res = [l[0],l[1],'-','-','%.5f' %(0.0),gene,aa,la]
		o.write('\t'.join(res)+'\n')
		continue
	if pos[0][7]=='-':
		res = [l[0],l[1],pos[0][4],pos[0][7],'%.5f' %(0.0),gene,aa,la]
		o.write('\t'.join(res)+'\n')
		continue
	elif int(pos[0][4]) < cRec:
		res = [l[0],l[1],pos[0][4],pos[0][7],pos[0][8],gene,aa,la]
		o.write('\t'.join(res)+'\n')
		continue
	mm=pos[0][7].split()[0]
	if mm not in ['AG','TC']: continue
	c=eval(pos[0][6])
	if mm=='AG':
		AG[0]+=c[0]
		AG[1]+=c[2]
		try: ed=float(c[2])/(c[2]+c[0])
		except: ed=0.0
	elif mm=='TC':
		AG[0]+=c[3]
		AG[1]+=c[1]
		try: ed=float(c[1])/(c[1]+c[3])
		except: ed=0.0
	res = [l[0],l[1],pos[0][4],'AG','%.5f' %(ed),gene,aa,la]
	o.write(' '.join(res)+'\n')
f.close()

try: REI=float(AG[1])/(sum(AG))
except: REI=0.0
o.write('>REI %.10f %.10f\n' %(REI,REI*100,))
table.close()
o.close()
sys.stderr.write('ALL done.\n')





