import math
import matplotlib.pyplot as plt

def SASA_angstrom(dico):
    """
    Calculate the total solvent-accessible surface area (SASA) in angstroms.

    Parameters:
    - dico (dict): A dictionary containing the SASA values for each residue.

    Returns:
    - sasa (float): The total SASA in angstroms.
    """
    sasa = 0
    for res in dico.keys():
        sasa += dico[res]['SASA']
    return sasa

def occurence(dico):
    """
    Count the occurrence of each amino acid in a dictionary.

    Args:
        dico (dict): A dictionary containing information about amino acids.

    Returns:
        dict: A dictionary containing the count of each amino acid.

    """
    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    amino_counts={aa:0 for aa in amino_acids}
    for res in dico.keys():
       amino_counts[dico[res]['resname']]+=1
    return amino_counts

def surface(dico_m1, dico_c):
    """
    Calculate the surface and interface residues based on the given dictionaries.

    Args:
        dico_m1 (dict): Dictionary containing information about chains and residues.
        dico_c (dict): Dictionary containing information about chains and residues.

    Returns:
        tuple: A tuple containing two dictionaries. The first dictionary represents the surface residues,
               and the second dictionary represents the interface residues.
    """
    surface = {}
    interface = {}
    for c in dico_m1['chains']:
        for res in dico_m1[c]['reslist']:
            if dico_m1[c][res]['rSASA'] > 0.25:
                surface[res] = dico_m1[c][res]
                Delta = dico_c[c][res]['rSASA'] - dico_m1[c][res]['rSASA']
                if Delta < 0:
                    interface[res] = dico_m1[c][res]
    return surface, interface


def deltaRSASA(dico_m1, dico_c):
    """
    Calculate the difference in relative solvent accessible surface area (rSASA) between two dictionaries.

    Parameters:
    - dico_m1 (dict): The first dictionary containing the rSASA values for each residue in each chain.
    - dico_c (dict): The second dictionary containing the rSASA values for each residue in each chain.

    Returns:
    - delta (dict): A dictionary containing the difference in rSASA values for each residue that has a rSASA value greater than 0.25.

    """
    delta = {}
    for c in dico_m1['chains']:
        for res in dico_m1[c]['reslist']:
            if dico_m1[c][res]['rSASA'] > 0.25:
                Delta = dico_c[c][res]['rSASA'] - dico_m1[c][res]['rSASA']
                delta[res] = Delta
    return delta



def distance_mono_complex(dico_m1, dico_m2, dico_c):
    """
    Calculate the difference in euclidean distance between the monomer and complex structures.
    """
    c1,c2=dico_c['chains']
    distances_chain1 = {}
    distances_chain2 = {}

    # Calculate distances for the first chain
    i=0
    for res in dico_m1[c1]['reslist']:
        if res in dico_c[c1]['reslist']:
            for atom in dico_m1[c1][res]['atomlist']:
                if atom in dico_c[c1][res]['atomlist']:
                    x1 = dico_m1[c1][res][atom]['x']
                    y1 = dico_m1[c1][res][atom]['y']
                    z1 = dico_m1[c1][res][atom]['z']
                    x2 = dico_c[c1][res][atom]['x']
                    y2 = dico_c[c1][res][atom]['y']
                    z2 = dico_c[c1][res][atom]['z']
                    distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                    distances_chain1[i] = distance
                    i+=1

    # Calculate distances for the second chain
    j=0
    for res in dico_m2[c2]['reslist']:
        if res in dico_c[c2]['reslist']:
            for atom in dico_m2[c2][res]['atomlist']:
                if atom in dico_c[c2][res]['atomlist']:
                    x1 = dico_m2[c2][res][atom]['x']
                    y1 = dico_m2[c2][res][atom]['y']
                    z1 = dico_m2[c2][res][atom]['z']
                    x2 = dico_c[c2][res][atom]['x']
                    y2 = dico_c[c2][res][atom]['y']
                    z2 = dico_c[c2][res][atom]['z']
                    distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                    distances_chain2[j] = distance
                    j+=1

    return distances_chain1, distances_chain2
            

