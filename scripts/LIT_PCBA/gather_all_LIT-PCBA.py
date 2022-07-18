#!/usr/bin/env python
import os
import rdkit.Chem as Chem
from os.path import join
import pandas as pd
from futils import ROOT
from pathlib import Path

DEFAULT_INPUT_DIR = ROOT / "data/gnina_processed/LIT-PCBA"
# DEFAULT_INPUT_DIR = join(os.path.dirname(__file__), DEFAULT_INPUT_DIR)
DEFAULT_RAW_DIR = ROOT / "data/raw_data/LIT-PCBA/AVE_unbiased/"
# DEFAULT_RAW_DIR = join(os.path.dirname(__file__), DEFAULT_RAW_DIR)
def main(input_dir = DEFAULT_INPUT_DIR, raw_dir = DEFAULT_RAW_DIR):
    walker = os.walk(input_dir)
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
    full_df.to_feather(join(input_dir, "full.feather"))

    full_df_filtered = full_df.sort_values("CNN_VS", ascending = False).drop_duplicates(["_Name", "directory", "ligands"])
    full_df_filtered["active_groundtruth"] = full_df_filtered["ligands"].str.contains("_active_")
    full_df_filtered.reset_index().rename(columns = {"level_0": "unfiltered_index"}).to_feather(join(input_dir, "filtered.feather"))

    n_df = {}
    for curdir, folders, files in os.walk(raw_dir):
        smiles = [file for file in files if file.endswith(".smi")]
        if not smiles:
            continue
        actives = 0
        inactives = 0
        for smile in smiles:
            if "_active_" in smile:
                actives += len(Chem.SmilesMolSupplier(join(curdir, smile), titleLine = False))
            else:
                inactives += len(Chem.SmilesMolSupplier(join(curdir, smile), titleLine = False))
        total = actives + inactives
        n_df[curdir.split(os.path.sep)[-1]] = (actives, total)

    n_df = pd.DataFrame(n_df).T.rename(columns = {0: "actives", 1: "total"})
    n_df.to_csv(join(raw_dir, "n.csv"))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--input_dir", default = DEFAULT_INPUT_DIR, help="input directory")
    p.add_argument("--raw_dir", default = DEFAULT_RAW_DIR, help="raw directory")
    args = p.parse_args()
    main(args.input_dir, args.raw_dir)
# %%
