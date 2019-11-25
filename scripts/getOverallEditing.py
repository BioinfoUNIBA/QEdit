
import sys, os
import pysam

try:
	edtablefile=sys.argv[1]
	redifile=sys.argv[2]
except:
	sys.exit('<REDItool table indexed by Tabix> <REDIportal file>')

try: tab=pysam.Tabixfile(edtablefile)
except: sys.exit('REDItools table not indexed.')
	
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
sys.stderr.write('TableName #AllPositions #UsedPositions OverallEditing OverallEditingx100\n')
sys.stderr.write('%s %i %i %.3f %.3f\n' %(edtablefile,tPos,nPos,oved,oved*100))


