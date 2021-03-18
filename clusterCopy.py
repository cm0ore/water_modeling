#use Kmeans clustering to cluster all of the water coordinates for each ensemble in k clusters. Across the ligand series, update the center xyz coordinates of each cluster depending on the positions of waters in that cluster. 
#2/13/21
#Camille Moore

from __future__ import division
from __future__ import print_function
import iotbx.pdb
from scipy import stats
import numpy as np
import sys
import os
import math
import iotbx.pdb.amino_acid_codes
from scipy.cluster.vq import vq, kmeans, whiten
from scipy.spatial import KDTree


def cluster_maker(pdb, initialArray):
  global apo_waterList
  pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdb)
  water_list = []
  waterArrays = {}
  water_dic = {}
  waterModelCodebook = []  #should ultimately be the length of all the waters and have elements=model.id of each water
  for model in pdb_obj.hierarchy.models():
    water_dic[model.id], waterArrays[model.id] = 0,[]
    for chain in model.chains():
      for rg in chain.residue_groups():
        if rg.atom_groups()[0].resname in aa_resnames: continue
        for ag in rg.atom_groups():
          if ag.resname == 'HOH':
            for atom in ag.atoms():
              water_dic[model.id] += 1
              waterModelCodebook.append(model.id)
              waterArrays[model.id].append(atom.xyz)
              water_list.append(atom.xyz)
  
  if file_base != '7KQO':  #the label of your apo pdb
    water_list = np.array(water_list)
    water_list = np.vstack((water_list, apo_waterList))
  if file_base == '7KQO':
    max_model = max(water_dic, key=lambda key: water_dic[key])
    max_waters = water_dic[max_model]
    print('maximum waters is %s, making %s clusters' % (max_waters, max_waters))
    water_list = np.array(water_list)
    apo_waterList = water_list
    initialArray = np.array(waterArrays[max_model])
  codebook, distortion = kmeans(water_list, initialArray)
  print('num clusters is:', codebook.shape)
  return codebook, water_list, waterModelCodebook, initialArray, apo_waterList

def clusterVariance(codebook, water_list, waterModelCodebook):
  varCodebookIndex, varMatrix = vq(water_list, codebook) 
  num_rows = codebook.shape[0] #gives number of clusters
  watersPerClus = [0]*num_rows 
  
  bigVarList = [0]*num_rows
  for num in range(varCodebookIndex.shape[0]):    #iteratre through all the waters in water_list
    clusterIndex = varCodebookIndex[num] #gives proper cluster index to add the individual water variances to 
    variance = varMatrix[num]
    bigVarList[clusterIndex] += variance  #ultimately each bigVarList[clusterIndex] entry will have the summed positional variance of each cluster from all waters in that cluster
    watersPerClus[clusterIndex] += 1
  
  bigVarList = np.array(bigVarList)/np.array(watersPerClus)     #now each cluster has the positional variance 
  with open('water_clusters.pml', "a") as file:
    for num in range(num_rows):
      print('pseudoatom cluster_%s, chain=%s, state=-1, b=%.3f, pos=[%s,%s,%s]' % (num, file_base, watersPerClus[num], codebook[num][0], codebook[num][1],codebook[num][2]), file=file)
  file.close()
   
  return watersPerClus, bigVarList  #watersPerClus a vector of dimension cluster_number (420) that has the occupancy (number of waters) of each cluster. 


def protein_coords(pdb):
  pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdb)
  protein_list = []
  protein_dic = {}
  for model in pdb_obj.hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        if rg.atom_groups()[0].resname not in aa_resnames: continue
        for ag in rg.atom_groups():
          if chain.id == 'A' and int(rg.resseq.strip()) < 3: continue    #because only the apo and 7kqp have the first two residues for some reason
          if ag.resname[0] != 'H':         #heavy atoms only
            for atom in ag.atoms():
              protein_list.append(atom.xyz)
              #atom_name = (ag.resname, rg.resseq, chain.id, atom.name.strip())   #format is ('ALA', '21', 'A', 'CA')
              atom_name = (ag.resname, rg.resseq.strip(), chain.id)
              protein_dic[atom.xyz] = atom_name 
  protein_list = np.array(protein_list)
  print('number of models:', model.id)
  return protein_list, protein_dic, model.id


def colorCoder(interactionScoreDic, rmsfArray, waterArray):    #make arrays of residue interactions for each cluster coded by the stability and consistency of the residues 
  waterStdDev = np.std(waterArray, axis=0)
  rmsfStdDev = np.std(rmsfArray, axis=0)
  rmsfAvg = np.average(rmsfArray, axis=0)
  waterAvg = np.average(waterArray, axis=0)
  variableArray = np.empty([3,int(waterArray.shape[1])])
  stableArray = np.empty([3,int(waterArray.shape[1])])
  unstableArray = np.empty([3,int(waterArray.shape[1])])

  for clusterNum in interactionScoreDic:
    if waterStdDev[clusterNum] > 0.3: #i only want to focus on the varying clusters
      array = variableArray
    elif waterStdDev[clusterNum] < 0.21 and waterAvg[clusterNum] > 0.72: #stable cluster
    #if waterAvg[clusterNum] > 0.72: 
      array = stableArray
    elif waterStdDev[clusterNum] < 0.21 and waterAvg[clusterNum] < 0.72: #consistently unoccupied cluster
    #elif waterAvg[clusterNum] < 0.72:
      array = unstableArray
    else:
      continue 
    for interaction in interactionScoreDic[clusterNum]:
       
      proteinID = interaction[0]  #('ALA', '21', 'A')
      numInteractions = interaction[1]
      if proteinID[2] == 'A':
        chainNum = int(proteinID[1]) -1  
      elif proteinID[2] == 'B':
        chainNum = int(proteinID[1]) + 166

      if rmsfStdDev[chainNum] > 0.04:
        array[0][clusterNum] += numInteractions
      elif rmsfStdDev[chainNum] < 0.015:
      #if 10 > 9:	 
        if rmsfAvg[chainNum] < 3.05:
          array[1][clusterNum] += numInteractions
        else:  
          array[2][clusterNum] += numInteractions

      
  print(variableArray, variableArray.shape)
  return variableArray, stableArray, unstableArray 
  #return stableArray, unstableArray 


def interactionScore(tree, waterCluster, waterOccupancy, models, InteractionFrequency):
 
  interactions = {}     #clustered water number  --> [(('ALA', '21', 'A'), (('GLY', '20', 'A')] aka closest residue to those coordinates 
  numb = 0
  count = -1
  for water_coordinates in waterCluster:
    count += 1
    if waterOccupancy[count] == 0:
      print('no waters in cluster', count, ', skipping')
      continue
    interactions[numb] = []
    nearestResidueDistance, nearestResidueIndex = tree.query(water_coordinates, int(models))  #we want the n nearest protein residues, where n=number of structures in the ensemble
    for num in range(int(models)):
      distance = nearestResidueDistance[num]
      if distance <= 3:
        res_coords = tree.data[nearestResidueIndex[num]]
        res_coords = tuple(res_coords)
        interactions[numb].append(resids[res_coords])
    numb += 1

#record frequency of interaction for each residue on a per coordinate basis 
  for coordinate in interactions:
    frequency = {}
    for residue in interactions[coordinate]:   #interactions has all the interactions of the current ensemble
      if residue not in frequency:
        frequency[residue] = interactions[coordinate].count(residue)  #frequency of each residue interaction for this cluster for this ensemble
        frequency[residue] /= int(models)  #norm the frequency by the number of models
      added = 0
      for resTuple in InteractionFrequency[coordinate]:  #the big list of interactions for this cluster
        if residue == resTuple[0]:
          resTuple[1] += frequency[residue]
          added = 1
       
      if added != 1:
        InteractionFrequency[coordinate].append([residue, frequency[residue]])
      
  return InteractionFrequency


################################
aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
hydrophobic_rank ={'CYS':1, 'ILE':2, 'PHE':2, 'VAL':3, 'LEU':4, 'MET':4, 'TRP':4, 'HIS':5, 'TYR':6, 'ALA':7, 'GLY':8, 'THR':9, 'SER':10, 'PRO':11, 'ARG':11, 'ASN':12, 'GLN':13, 'ASP':13, 'GLU':13, 'LYS':14}
initialArray = 0
apo_waterList = 0
apoBfactors = np.array([0]*420)
waterArray = []
waterOccupancy = []
InteractionFrequency = {} #water cluster number --> [[('ALA', '21', 'A'), number_of_interactions], (('GLY', '20', 'A'), number_of_interactions)]
for cluster in range(420):
  InteractionFrequency[cluster] = []

#rmsfArray = np.load('/wynton/home/rotation/cmoore/best_ensembles/bigRmsfArray.npy')
pdb_list = open(sys.argv[1])

big_x = []
big_y = []
big_z = []

for line in pdb_list:
  line = line.strip()
  print(line)
  file_name = line
  file_base = os.path.basename(file_name).split('_')[0]  
  rmsf_file = open('/wynton/home/rotation/cmoore/best_ensembles/rmsfs/%s_rmsf.csv' %file_name)
  aligned_file = ('%s_aligned.pdb' %file_base)
  
  cluster_coords, water_list, waterModelCodebook, initialArray, apo_waterList = cluster_maker(aligned_file, initialArray)
  big_x.append(cluster_coords[:,0].reshape(420))
  big_y.append(cluster_coords[:,1].reshape(420))
  big_z.append(cluster_coords[:,2].reshape(420))
  
  waterBfactors, waterVariance = clusterVariance(cluster_coords, water_list, waterModelCodebook)
  if file_base == '7KQO':
    apoBfactors = np.array(waterBfactors)
    waterBfactors = np.array(waterBfactors)
  protein_xyzs, resids, models = protein_coords(aligned_file)
  waterArray.append(waterVariance)
  if file_base != '7KQO':
    waterBfactors = np.array(waterBfactors)
    waterBfactors = waterBfactors - apoBfactors
  models = int(models)
  waterBfactors = np.true_divide(waterBfactors, models)
  #print(waterBfactors)  #just normed occupancy for now 
  waterOccupancy.append(waterBfactors)

  #tree = KDTree(protein_xyzs)
  #InteractionFrequency = interactionScore(tree, cluster_coords, waterBfactors, models, InteractionFrequency)  
  
pdb_list.close()
st = ['big_x', 'big_y', 'big_z']
num = 0
for array in [big_x, big_y, big_z]:
  array = np.row_stack(array)
  np.save('%s.npy' % st[num], array)
  num += 1

waterOccupancy = np.row_stack(waterOccupancy)
variable, stable, unstable = colorCoder(InteractionFrequency, rmsfArray, waterOccupancy)
np.save('varInteractions.npy', variable)
np.save('stableInteractions.npy', stable)
np.save('unstableInteractions.npy', unstable)




