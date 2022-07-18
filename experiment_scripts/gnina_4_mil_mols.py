#%%
import os
from os.path import join
from futils import PROT

path = "../data/equibind_processed/4_mil_mol/chopped"
outdir = "../data/gnina_processed/4_mil_mol/chopped"
dir, _, files = next(os.walk(path))
for file in files:
    os.system(f"gnina -l {join(dir, file)} -r {PROT} -o {join(outdir, file)} --minimize > {join(outdir, file.split('.')[0])}.txt")


# %%
