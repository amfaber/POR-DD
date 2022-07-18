#!/usr/bin/env python
import rdkit.Chem as Chem
import os
from os.path import join
import pandas as pd
from futils import ROOT

DEFAULT_INPUT_DIR = ROOT / "data/gnina_processed/LIT-PCBA"

def main(inputdir = DEFAULT_INPUT_DIR, targets = None):
    walker = os.walk(inputdir)
    for curdir, folders, files in walker:
        if targets is not None and not curdir.split(os.path.sep)[-3] in targets:
            continue
        if not "gnina_out.sdf" in files:
            continue
        print(f"Starting {os.path.abspath(curdir)}")
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

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--input_dir", default = DEFAULT_INPUT_DIR, help="input directory")
    p.add_argument("--targets", default = None, nargs = "*", help = "Only process some of the targets ")
    args = p.parse_args()
    main(args.input_dir, args.targets)