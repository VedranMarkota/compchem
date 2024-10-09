#!/usr/bin/env python3

#Adapted from https://www.macinchem.org/reviews/clustering/clustering.php 

#Here, we are generating Morgan fingerprints and clustering them. Then, the fingerprint leaders(centroids) are extracted by their indices from the original sdf file

from rdkit import Chem
from rdkit.Chem import AllChem

#Define clustering setup
def ClusterFps(fps,cutoff=0.7):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

#Import molecules; Renamed_Filtered.sdf is a molecular library file 
ms = [x for x in Chem.SDMolSupplier('Renamed_Filtered.sdf',removeHs=False)]

#loading sdf file and generating Morgan fingerprints
fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,2048) for x in Chem.ForwardSDMolSupplier('Renamed_Filtered.sdf',removeHs=False) if x is not None]

#clusters by fingerprints, cutoff is similarity score between 0 and 1, default is 0.7
clusters=ClusterFps(fps,cutoff=0.7)

print(clusters)

#prints the fist sub element of every element in the file
for i in clusters:
    print(i[0])

#creates an array from the first elements
a=[]
for i in clusters:    
    print(i[0])
    a.append(i[0])

#pick molecules from the original file Renamed_Filtered.sdf by indices and puts them into the new file RenamedFiltered_Centroids.sdf
picks = [ms[x] for x in a]

writer = Chem.SDWriter('RenamedFiltered_Centroids.sdf')
for mol in picks:
    writer.write(mol)
