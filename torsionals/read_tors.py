import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import sys, pymol
from pymol import cmd, stored
stdout = sys.stdout
pymol.finish_launching()
sys.stdout = stdout 
import matplotlib.pyplot as plt
plt.style.use(['seaborn-darkgrid', 'seaborn-colorblind'])
import seaborn as sns
import numpy as np
import os

def get_phi(obj, resi, state):
    """
    Get the dihedral angle phi (OR-C1-O'x-C'x) 
    """
    atom1 = '%s and resi %s and name C' % (obj, resi-1)
    atom2 = '%s and resi %s and name N' % (obj, resi)
    atom3 = '%s and resi %s and name CA' % (obj, resi)
    atom4 = '%s and resi %s and name C' % (obj, resi)
    value = cmd.get_dihedral(atom1, atom2, atom3, atom4, state=state)
    return value


def get_psi(obj, resi, state): 
    """
    Get the dihedral angle psi (C1-O'x-C'x-C'x-1)
    """
    atom1 = '%s and resi %s and name N' % (obj, resi)
    atom2 = '%s and resi %s and name CA'% (obj, resi)
    atom3 = '%s and resi %s and name C' % (obj, resi)
    atom4 = '%s and resi %s and name N' % (obj, resi+1)
    value = cmd.get_dihedral(atom1, atom2, atom3, atom4, state=state)
    return value 

def get_chi(obj, resi, state): 
    """
    Get the dihedral angle chi (C1-O'x-C'x-C'x-1)
    """
    atom1 = '%s and resi %s and name N' % (obj, resi)
    atom2 = '%s and resi %s and name CA'% (obj, resi)
    atom3 = '%s and resi %s and name CB' % (obj, resi)
    atom4 = '%s and resi %s and name CG' % (obj, resi)
    value = cmd.get_dihedral(atom1, atom2, atom3, atom4, state=state)
    return value 

def read_torsionals(resi, pdb_file_name):
    '''
    Read torsional angles phi and psi for a given residue in a protein ensamble.
    ----------
    Parameters:
    ----------
    resi: int. Residue number in a protein sequence
    pdb_file_name: string. PDB file 
    '''
    
    cmd.load(pdb_file_name)
    pose = os.path.splitext(pdb_file_name)[0]
    n_states = cmd.count_states('all')
    tors = np.zeros((n_states, 2))
    
    for s in range(1, n_states+1):
        phis = get_phi(pose, resi, s)
        psis = get_psi(pose, resi, s)
        #chis = get_chi(pose, resi, s)
        tors[s-1] = [phis, psis]        
    
    plt.figure(figsize=(18, 4))
    plt.subplot(121)
    mu = tors[:,0].mean()
    sigma = tors[:,0].std()
    sns.kdeplot(tors[:,0], gridsize=50000)
    plt.plot(0, label='$\\mu$ = {:.3f}\n$\\sigma$ = {:.3f}'.format(mu, sigma), alpha=0)
    plt.xlim(-180,180)
    plt.xlabel('$\\phi$')
    plt.legend(loc='best', fontsize=14)
    plt.subplot(122)
    mu = tors[:,1].mean()
    sigma = tors[:,1].std()
    sns.kdeplot(tors[:,1], gridsize=50000)
    plt.plot(0, label='$\\mu$ = {:.3f}\n$\\sigma$ = {:.3f}'.format(mu, sigma), alpha=0)
    plt.xlim(-180, 180)
    plt.xlabel("$\\psi$")
    plt.legend(loc='best', fontsize=14)
    plt.savefig("kde_{}_torsionals_residue_{}.png".format(pose, resi), dpi = 300)

    return tors, pose

tors, pose = read_torsionals(25, '2k39.pdb')
