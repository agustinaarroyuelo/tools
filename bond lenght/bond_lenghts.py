import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import sys, pymol
from pymol import cmd
stdout = sys.stdout
pymol.finish_launching()
sys.stdout = stdout 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# Compute bond length for C $\alpha$ - H $\alpha$, C $\alpha$ - N y C $\alpha$ - C $\beta$

def get_bond_lenghts(atom_name, pdb_file):
    try:
        cmd.load(pdb_file)
        n_states = cmd.count_states('all')
        n_residues = cmd.count_atoms('name CA')
        distance = np.zeros((n_states, n_residues))

        for s in range(1,n_states+1):    
            conf = []
            for i in range(1, n_residues+1):
                try:
                    dst = cmd.get_distance('resi {} and name CA'.format(i), 
                                           'resi {} and name {}'.format(i, atom_name), state=s)
                    conf.append(dst)
                except:
                    conf.append(np.nan)
            distance[s-1] = conf
        return distance
    except:
        print('Cannot open file {}. Please try a diferent PDB file.'.format(pdb_file))

def plot_lenghts(protein_name, atom_name, distance):
    
    if atom_name == "C":
        distance = distance[:,:-1 ]
    elif atom_name == 'N':
        distance = distance[:,1:distance.shape[1]]
    else:
        distance = distance[~np.isnan(distance)].reshape(distance.shape[0], -1)
        
    kde = distance.flatten()
    mu = kde.mean()
    sigma = kde.std()
    
    plt.figure(figsize=(18, 4))
    
    plt.subplot(121)
    plt.title('{} - Bond lenghts'.format(protein_name))
    sns.tsplot(data = distance, interpolate=False, err_style="ci_bars", ci=50)
    plt.ylabel("C$\\alpha$ - {} Distance ($\\AA$)".format(atom_name)) 
    plt.xlabel("Residue")

    plt.subplot(122)
    plt.title('{} - Bond lenght Density'.format(protein_name))
    sns.kdeplot(kde)
    plt.xlim(mu - sigma*4, mu + sigma*4)
    plt.plot(0, label='$\\mu$ = {:.3f}\n$\\sigma$ = {:.3f}'.format(mu, sigma), alpha=0)
    plt.xlabel("C$\\alpha$ - {} Distance ($\\AA$)".format(atom_name))
    plt.legend(loc='best', fontsize=14)
    plt.show()

dist = get_bond_lenghts('HA','2k39.pdb')

plot_lenghts('Ubiquitin', 'HA', dist)

