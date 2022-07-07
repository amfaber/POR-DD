#!/usr/bin/env python
import os
from openbabel import pybel
import Bio.PDB as bp
import rdkit.Chem as Chem
from os.path import join
import pandas as pd

inputdir = "../../data/gnina_processed/LIT-PCBA"
rawdir = "../../data/raw_data/LIT-PCBA/AVE_unbiased/"
walker = os.walk(inputdir)
full_df = pd.DataFrame()

for curdir, folders, files in walker:
    if "gnina_properties.feather" in files:
        df = pd.read_feather(join(curdir, "gnina_properties.feather"))
        dir, ligs, prot = curdir.split(os.path.sep)[-3:]
        df["directory"] = dir
        df["ligands"] = ligs
        df["protein"] = prot
        full_df = pd.concat((full_df, df))
full_df.reset_index(inplace=True)
full_df = full_df.astype({"minimizedAffinity": float, "minimizedRMSD": float, "CNNscore": float, "CNNaffinity": float, "CNN_VS": float, "CNNaffinity_variance": float, "index": int})
full_df.to_feather(join(inputdir, "full.feather"))

full_df_filtered = full_df.sort_values("CNN_VS", ascending = False).drop_duplicates(["_Name", "directory", "ligands"])
full_df_filtered["active_groundtruth"] = full_df_filtered["ligands"].str.contains("_active_")
full_df_filtered.reset_index().rename(columns = {"level_0": "unfiltered_index"}).to_feather(join(inputdir, "filtered.feather"))


n_df = {}
for curdir, folders, files in os.walk("../../data/raw_data/LIT-PCBA/AVE_unbiased"):
    smiles = [file for file in files if file.endswith(".smi")]
    if not smiles:
        continue
    actives = 0
    inactives = 0
    for smile in smiles:
        if "_active_" in smile:
            actives += len(Chem.SmilesMolSupplier(join(curdir, smile)))
        else:
            inactives += len(Chem.SmilesMolSupplier(join(curdir, smile)))
    total = actives + inactives
    n_df[curdir.split(os.path.sep)[-1]] = (actives, total)

n_df = pd.DataFrame(n_df).T.rename(columns = {0: "actives", 1: "total"})
n_df.to_csv(join("../../data/raw_data/LIT-PCBA/AVE_unbiased", "n.csv"))