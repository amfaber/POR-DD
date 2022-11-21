from rdkit.Chem import Draw
import numpy as np 
def mols_to_tiff(save_location, mols, square_side, clearconfs = False):
    if clearconfs:
        print("Removing conformation of every mol in provided list IN PLACE.")
        [mol.RemoveConformer(0) for mol in mols]
    pngs = []
    for i in range(int(np.ceil(len(mols)/square_side**2))):
        png = Draw.MolsToGridImage(mols[i*square_side**2: (i+1)*square_side**2], molsPerRow = square_side, returnPNG = False,
        legends = list(map(str, range(i*square_side**2, (i+1)*square_side**2))))
        pngs.append(png)
    pngs[0].save(save_location, save_all = True, append_images = pngs[1:])