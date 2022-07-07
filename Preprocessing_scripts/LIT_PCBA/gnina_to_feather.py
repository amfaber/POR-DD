#!/usr/bin/env python
#%%
import rdkit.Chem as Chem
import os
from os.path import join
import pandas as pd
# %%
inputdir = "../data/gnina_processed/LIT-PCBA"
walker = os.walk(inputdir)
for curdir, folders, files in walker:
    if not curdir.split(os.path.sep)[-1] == "OPRK1":
        continue
    if not "gnina_out.sdf" in files:
        continue
    datapath = join(curdir, "gnina_out.sdf")
    supp = Chem.SDMolSupplier(datapath, sanitize = False)
    all_data = []
    for i, mol in enumerate(supp):
        all_data.append({})
        if mol is None:
            name = supp.GetItemText(i).split("\n")[0]
            all_data[-1]["_Name"] = name
            all_data[-1]["_MolFileComments"] = f"{i}"
            continue
        for prop in mol.GetPropNames(includePrivate = True):
            all_data[-1][prop] = mol.GetProp(prop)
    df = pd.DataFrame(all_data)
    df.to_feather(join(curdir, "gnina_properties.feather"))

# %%
