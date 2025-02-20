#!/usr/bin/env python
# coding: utf-8

import sys
import math
import numpy as np
from matplotlib import pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

#manually setup topology, trajectory file and atoms involved in dihedrals
u = mda.Universe('dry.pdb', 'dry.dcd')
#print(len(u.trajectory))
c1 = 'resname LIG and name C10'
c2 = 'resname LIG and name N1'
c3 = 'resname LIG and name C4'
c4 = 'resname LIG and name O2'

def calculate_dihedrals_manually(atoms_list):
    ag1 = u.select_atoms(atoms_list[0])
    ag2 = u.select_atoms(atoms_list[1])
    ag3 = u.select_atoms(atoms_list[2])
    ag4 = u.select_atoms(atoms_list[3])
    angle = calc_dihedrals(ag1.positions, ag2.positions,
                               ag3.positions, ag4.positions,
                               box=ag1.dimensions)
    angle = np.rad2deg(angle)
    return angle

v0_list = [c1, c2, c3, c4]
v0_arr = [] 
for ts in u.trajectory:
    v0_angle = float(calculate_dihedrals_manually(v0_list))
    v0_arr.append(v0_angle)
    
with open('glycosidic-bond.txt','w') as f:
    for i, v0 in enumerate(v0_arr):
        output = str(i) + ',' + str(v0) + '\n'
        f.write(output)

frames = np.genfromtxt('glycosidic-bond.txt', usecols = 0, delimiter = ',')
v0 = np.genfromtxt('glycosidic-bond.txt', usecols = 1, delimiter = ',')
time = np.arange(1, len(frames)+1)*2
fig, ax = plt.subplots()
plt.axhspan(-90, 90, label='sin', alpha=0.2, facecolor ='white')
plt.axhspan(-180, -90, label='anti', alpha=0.2, facecolor ='chartreuse')
plt.axhspan(90, 180, alpha=0.2, facecolor ='chartreuse')
sc = ax.scatter(time, v0, marker='o', color='k')
plt.ylabel(r'$\chi$ (°)')
plt.xlabel('Time (ns)')
plt.xlim(0, 500)
plt.ylim(-180, 180)
plt.yticks([-180, -90, 0, 90, 180]) 
plt.legend()
#plt.tight_layout()
plt.title('MRS2905 base conformation', fontweight='bold')
plt.savefig('glycosidic-bond.png', dpi=300)
