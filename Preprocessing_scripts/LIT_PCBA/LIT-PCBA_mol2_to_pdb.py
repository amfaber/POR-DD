# %%
import os
from openbabel import pybel
from pathlib import Path

# %%

def convert(LIT_PCBA_path):
    walker = os.walk(LIT_PCBA_path)
    next(walker)
    for curdir, dirs, files in walker:
        prots = [Path(file) for file in files if file.endswith("protein.mol2")]
        for prot in prots:
            protpath = curdir / prot
            mol = next(pybel.readfile("mol2", str(protpath)))
            outname = str(protpath).replace(".mol2", "_tmp.pdb")
            # print(outname)
            # return
            mol.write("pdb", outname, overwrite = True)
            os.system(f"reduce {outname} > {outname.replace('_tmp.pdb', '.pdb')}")
            os.system(f"rm {outname}")

#%%

if __name__ == "__main__":
    convert("../data/raw_data/LIT-PCBA/AVE_unbiased")
    # convert("../data/raw_data/LIT-PCBA/full_data")


# %%
