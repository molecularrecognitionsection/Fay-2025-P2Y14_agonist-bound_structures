from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolAlign
from rdkit.Chem import AllChem
from rdkit import RDConfig
from rdkit.Chem import PandasTools
import os
import pandas as pd
from rdkit.Chem import rdMolDescriptors
import matplotlib.pyplot as plt
import sys

def import_database():
    sdfFile =sys.argv[1]

    print("Database import...")

    df = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule',includeFingerprints=False,removeHs=False)

    df_2 = PandasTools.LoadSDF(sdfFile,molColName='Molecule_opt',includeFingerprints=False,removeHs=False)

    df_merge = pd.concat([df, df_2], axis=1)

    print("Done!")
    df_merge=df_merge.iloc[:5]
    return df_merge


def optimize(mol):
    value=AllChem.MMFFOptimizeMolecule(mol,maxIters=1000,mmffVariant="MMFF94s")
    return value

def rmsd(x, y):
    l=[]
    l_l=[]
    for i in range(x.GetNumAtoms()):
        l_for_tuple=[]
        l_for_tuple.append(i)
        l_for_tuple.append(i)
        t=tuple(l_for_tuple)
        l.append(t)
        l_l.append(l)
    try:
        rmsd_value=rdMolAlign.CalcRMS(x,y,map=l_l)
    except:
        rmsd_value=99
    return rmsd_value

def plot(df_merge):
    ax = df_merge['rmsd'].hist(rwidth=0.8,bins=50)  
    ax.set_xlabel("RMSD")
    ax.set_ylabel("Frequency")
    ax.grid(False)
    fig = ax.get_figure()
    fig.savefig('plot_rmsd_test.png',dpi=600)


def main():
    df_merge=import_database()
    print("Pose minimization, this will take some time...")
    df_merge["Optimization"] = df_merge['Molecule_opt'].apply(optimize)
    print("Minimization completed")
    df_merge['rmsd'] = df_merge.apply(lambda x: rmsd(x['Molecule'], x['Molecule_opt']), axis=1)
    plot(df_merge)
    PandasTools.WriteSDF(df_merge, 'out.sdf', molColName='Molecule_opt', idName=None, properties=list(df_merge.columns), allNumeric=False)

if __name__ == "__main__":
    main()
