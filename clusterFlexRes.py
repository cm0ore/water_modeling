from __future__ import division
from __future__ import print_function
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import sys
import os
import math
import matplotlib.ticker as plticker

def flexRes(csv):    #returns a list wtih all of the most flexible residues from the whole ensemble with (chain.id, resnum) 
  flexRes = []
  for line in csv:
    line = line.split(',')
    resID = (line[2].strip(),int(line[1])) #format is ('A', '45')
    if resID not in flexRes:
      flexRes.append(resID)

  return flexRes


def rmsfArray(pdb_list):
  bigRmsfArray = []
  pdbs = []
  for line in pdb_list:
    pdb = line.strip()
    rmsfArray = []
    pdbs.append(pdb)
    rmsfFile = open('/wynton/home/rotation/cmoore/best_ensembles/rmsfs/%s_rmsf.csv' % pdb)
    for line in rmsfFile:
      line = line.split(',')
      if line[0] == 'dataset':
        continue
      rmsf = float(line[3]) 
      rmsfArray.append(rmsf)
    rmsfArray = np.array(rmsfArray)
    bigRmsfArray.append(rmsfArray)
  bigRmsfArray = np.vstack(bigRmsfArray)
  return bigRmsfArray, pdbs     


pdb_list = open(sys.argv[1])
flexResFile = open(sys.argv[2])
flexResList = flexRes(flexResFile)
bigRmsfArray, pdbs = rmsfArray(pdb_list)

stdDevArray = np.std(bigRmsfArray, axis=0)
result = np.where(stdDevArray > 0.25) #the residues with the most standard deviation in RMSF
add3 = len(result[0])
add3 = [3]*add3
result = result[0] + add3
print(result) #prints the most variable residue numbers

stdDevArray = np.std(bigRmsfArray, axis=0) 

np.save('bigRmsfArray.npy', bigRmsfArray)
fig, ax = plt.subplots(figsize=(30, 5))
im = ax.imshow(bigRmsfArray, vmin=0, vmax=7, aspect="auto")
ax.set_title('RMSF per residue for each ensemble')
ax.set_yticks(np.arange(bigRmsfArray.shape[0]))
ax.set_yticklabels(pdbs)
loc = plticker.MultipleLocator(base=5.0) # this locator puts ticks at regular intervals
ax.xaxis.set_major_locator(loc)

ax.set_xlabel("Residue")
cbar = ax.figure.colorbar(im, orientation='vertical', pad=0.05, fraction=0.005)
cbar.ax.set_ylabel('RMSF', rotation=-90, va="bottom")
plt.savefig('seriesRMSFS.pdf', dpi=300)
plt.show()
flexResFile.close()  
pdb_list.close()

    
