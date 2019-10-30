
import sys, os
import pysam

try:
	edtablefile=sys.argv[1] #"/lustrehome/epicardi/home/sra/SRR607337/SRR607337_dna.txt.gz"
	redifile=sys.argv[2] #"/lustrehome/epicardi/home/FINALediting/TABLE1_complete_corrected2.txt" #sys.argv[2]
except:
	sys.exit('<REDItool table> <REDIportal file>')

def BaseCount(seq,ref,mfr,VNUC):
	b={'A':0,'C':0,'G':0,'T':0}
	subs=[]
	subv=[]
	for i in seq.upper():
		if b.has_key(i): b[i]+=1
	for i in b:
		if not b.has_key(ref): continue
		if b[i]!=0 and i!=ref:
			vv=float(b[i])/(b[i]+b[ref])
			subv.append((b[i],vv,ref+i))
	subv.sort()
	subv.reverse()
	for i in subv:
		if i[0]>=VNUC and i[1]>=mfr: subs.append(i[2])
	freq=0.0
	if len(subs)==0: subs.append('-')
	else: freq=subv[0][1]	
	return sum(b.values()),[b['A'],b['C'],b['G'],b['T']],' '.join(subs),'%.2f'%(freq)

tab=pysam.Tabixfile(edtablefile)

nGall=0
nAGall=0
nPos=0
tPos=0
f=open(redifile)
for i in f:
	if i.startswith('Region'): continue
	if i.startswith('chromosome'): continue
	l=(i.strip()).split('\t')
	tPos+=1
	chr=l[0]
	pos=int(l[1])
	#print pos
	try: res=[kk.split('\t') for kk in tab.fetch(reference=chr,start=pos-1,end=pos)]
	except: res=[]
	if len(res)==0:
		nGall+=0
		nAGall+=0
		nPos+=1
	else:
		ref=res[0][2]
		dis=eval(res[0][6])
		if ref=='A':
			ng=dis[2]
			na=dis[0]
			nGall+=ng
			nAGall+=(na+ng)
			nPos+=1
		elif ref=='T':
			ng=dis[1]
			na=dis[3]
			nGall+=ng
			nAGall+=(na+ng)
			nPos+=1
		else:
			continue
f.close()

try: oved=float(nGall)/nAGall
except: oved=0.0

sys.stderr.write('%s %i %i %f\n' %(edtablefile,tPos,nPos,oved))


