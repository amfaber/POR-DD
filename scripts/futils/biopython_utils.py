import Bio.PDB as bp
import numpy as np

def remove_hetatoms(prot):
    for chain in prot:
        hetres = [res for res in chain if res.id[0] != " "]
        [chain.detach_child(res.id) for res in hetres]

def load_prot(protpath):
    return bp.PDBParser().get_structure("", protpath)[0]

def atom_coords(prot, removeH = False):
    if removeH:
        return np.stack([atom.coord for atom in prot.get_atoms() if atom.element != "H"])
    else:
        return np.stack([atom.coord for atom in prot.get_atoms()])