# water_modeling
Scripts for modeling water after ensemble refinement. These were optimized to be used for a ligand series (aka comparing water positions/occupancy of many ensembles).

## 1. clusterCopy.py
This file uses the K-means algorithm to cluster waters for each ensemble in the ligand series. It requires a list of ALIGNED ensemble pdb files. These pdbs must be aligned in x,y,z space before running this script. I used pymol align to do this.
Also note that the number of clusters made by the algorithm is determined by the max number of waters in a single structure from the apo ensemble. The apo ensemble should be the first pdb in the list.

  ### Input: 
    * .txt file with list of aligned .pdbs. You need to run this script from the directory where these pdbs are. Example txt file is best_ensemble_filenames.txt
    
  ### Outputs: 
    * water_clusters.pml which makes a pseudoatom at each cluster center for all ensemble in the series (total pseudoatoms = numberOfPdbs * numberOfClusters).
      form is pseudoatom cluster_%s, chain=%s, state=-1, b=%.3f, pos=[%s,%s,%s] % (number, pdb, positionalVarianceofThatCluster, x, y, z)
    * big_x.npy, big_y.npy, and big_z.npy which are arrays with dimensions (numberOfPdbs, numberOfClusters) and each array has all the x or y or z coordinates of each cluster center for the ligand series
    * waterOccupancy.npy which has the occupancy of each cluster across the ligand series. Same dimensions as big_x or big_y or big_z. Occupancy is (numberOfWaters/cluster)/numberOfStructuresInEnsemble
 
 ## 2. ensem_protein_rmsf.py 
 This file calculates protein residue RMSFs for ensemble pdbs. It writes a csv or pml file so that residues will have b_factors=rmsf in pymol
   
   ### Inputs:
    * ensemble .pdb
    * "csv" or "pml" to determine output 
    * max_rmsf (optional)
   
   ### Outputs:
    * csv or pml file
   ```
   Example Usage: libtbx.python ensem_protein_rmsf.py example.pdb csv > example.csv
   ```
 ## 3. clusterFlexRes.py 
 This file makes a big RMSF array for the ligand series. Note that you have to have run ensem_protein_rmsf.py first in order to generate csv files for each 
 ensemble. This script can be run in jupyter (see Plots.ipynb) to additionally display a heatmap. 
  
   ### Inputs:
    * .txt file with list of aligned .pdbs. 
    * update the filepath to your directory to where your exampleRMSF.csv files are
    
   ### Outputs:
    * bigRmsfArray.npy which has the rmsf of each residue for each pdb in the ligand series. Dimensions are (numberOfPdbs, numberOfResidues).
    * seriesRMSFS.pdf which has the heatmap of the bigRMSFArray.
    
 ## 4. Plots.ipynb 
 This is a notebook that shows cool heatmap visuals and leiden clusters the big_x,big_y,big_z arrays. It also makes a file called 
 spatialColoredWater.pml that colors the pseudoatoms from water_clusters.pml according to the positional variance of that cluster across the ligand series. 
 It also scales the size of the pseudoatoms according to their occupancy. Run this notebook after you've run the above scripts.
 

    
