#%%
import os
from os.path import join
from futils import PROT, ROOT

path = ROOT / "data/equibind_processed/4_mil_mol/chopped"
outdir = ROOT / "data/gnina_processed/4_mil_mol/chopped"
dir, _, files = next(os.walk(path))
for file in files:
    os.system(f"gnina -l {join(dir, file)} -r {PROT} -o {join(outdir, file)} --minimize > {join(outdir, file.split('.')[0])}.txt")


# %%
