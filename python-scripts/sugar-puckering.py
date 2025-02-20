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
c1 = 'resname LIG and name C4'
c2 = 'resname LIG and name C9'
c3 = 'resname LIG and name C8'
c4 = 'resname LIG and name C5'
o = 'resname LIG and name O2'

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

def pseudorotational_angle(v1, v2, v3, v4, v0):
    p = ((float(v4)-float(v0))-(float(v3)-float(v1)))/(2*float(v2)*(math.sin(math.radians(36))+math.sin(math.radians(72))))
    if float(v2)<0:
            P = math.degrees(math.atan(p))+180
    elif p<0:
            P = math.degrees(math.atan(p))+360
    else:
            P = math.degrees(math.atan(p))
    return round(P, 2)

def maximum_puckering(P, v2):
	v_max = abs(round((math.degrees(math.radians(float(v2))/math.cos(math.radians(P)))),2))
	return v_max

v1_list = [o, c1, c2, c3]
v2_list = [c1, c2, c3, c4]
v3_list = [c2, c3, c4, o]
v4_list = [c3, c4, o, c1]
v0_list = [c4, o, c1, c2]
v1_arr = [] 
v2_arr = []
v3_arr = [] 
v4_arr = []
v0_arr = [] 
for ts in u.trajectory:
    v1_angle = float(calculate_dihedrals_manually(v1_list))
    v1_arr.append(v1_angle)
    v2_angle = float(calculate_dihedrals_manually(v2_list))
    v2_arr.append(v2_angle)
    v3_angle = float(calculate_dihedrals_manually(v3_list))
    v3_arr.append(v3_angle)
    v4_angle = float(calculate_dihedrals_manually(v4_list))
    v4_arr.append(v4_angle)
    v0_angle = float(calculate_dihedrals_manually(v0_list))
    v0_arr.append(v0_angle)

for v1, v2, v3, v4, v0 in zip(v1_arr, v2_arr, v3_arr, v4_arr, v0_arr):
    P = pseudorotational_angle(v1, v2, v3, v4, v0)
    v_max = maximum_puckering(P, v2)
    with open('sugar-puckering.txt','a') as f:
        output = str(P) + ',' + str(v_max) + '\n'
        f.write(output)

P = np.genfromtxt('sugar-puckering.txt', usecols = 0, delimiter = ',')
v_max = np.genfromtxt('sugar-puckering.txt', usecols = 1, delimiter = ',')
time = np.arange(1, len(P)+1)*2
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.set_theta_zero_location("N")
ax.set_theta_direction('clockwise')
ax.set_rlabel_position(90)
#colors = time / 0.2
sc = ax.scatter(np.deg2rad(P), v_max,marker='o', c=time, cmap='plasma')
cb = fig.colorbar(sc, pad=0.1)
cb.set_label('Time (ns)', rotation=270, labelpad = 15)
plt.title('MRS2905 ribose puckering', fontweight='bold')
plt.savefig('polarplot.png', dpi=300)




