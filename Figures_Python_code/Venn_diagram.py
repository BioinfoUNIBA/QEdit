import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

try:
	groupA = sys.argv[1]
	groupB = sys.argv[2]
except:
	sys.exit('<groupA, groupB>')

circleA = []
circleB = []
with open(groupA, 'r') as e:
	circleA = map(str.strip, e.readlines())
with open(groupB, 'r') as f:
	circleB = map(str.strip, f.readlines())
	
labelA = ' '.join(e.name.split('.')[0].split('_')[1:])
labelB = ' '.join(f.name.split('.')[0].split('_')[1:])

venn2([set(circleA), set(circleB)], \
set_labels= (labelA, labelB), set_colors=('orange','blue'))

if e.name.find('trimmed') != -1:
	title = labelA.split(' ')[0] + ' ' + labelA.split(' ')[2].upper() + " AG+TC"
else:
	title = labelA.split(' ')[0] + ' ' + 'UNTRIMMED ' + "AG+TC"
plt.suptitle(title, fontweight="bold")
plt.savefig('%sVS%s' %(labelA, labelB))
