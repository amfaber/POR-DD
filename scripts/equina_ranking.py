#!/usr/bin/env python
from equina import *
import sys
from rdkit import Chem
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type = int, default = None)
    parser.add_argument("--fraction", type = float, default = None)
    parser.add_argument("--sort_by", default = "CNN_VS")
    parser.add_argument("--just_gnina", action = "store_true")
    parser.add_argument("--split_gnina", type = int, default = None)

    known_args, rest = parser.parse_known_args()
    if not ((known_args.fraction > 0) and (known_args.fraction < 1)):
        print("fraction has to between 0 and 1")
        sys.exit(1)
    args, cmdlineargs = multiligand_inference.parse_arguments(rest)
    args = multiligand_inference.get_default_args(args, cmdlineargs)
    
    lig_extension = args.ligands_sdf.split(".")[1]
    if lig_extension == "sdf":
        supp = Chem.SDMolSupplier(args.ligands_sdf)
    elif lig_extension == "smi":
        supp = Chem.SmilesMolSupplier(args.ligands_sdf, titleLine = False)
    else:
        print(f"unrecognized file format: {lig_extension}")
        exit(1)

    equina(args, make_parquet = True, just_gnina = known_args.just_gnina, split_gnina = known_args.split_gnina)
    pathoutput = Path(args.output_directory)
    # results = pd.DataFrame([(float(mol.GetProp("CNNaffinity")), i) for i, mol in enumerate(supp)])
    
    results = pd.read_parquet(str(pathoutput / "props.parquet"))
    results.sort_values(known_args.sort_by, inplace = True, ascending = False)
    top_n = int(known_args.fraction * len(results))
    best_indexes = np.array(results.index[:top_n])
    with open(pathoutput / "success.txt") as file:
        orig_idx = np.array([int(line.split(" ")[0]) for line in file.readlines()])
    new_output = str(pathoutput / "refinement")
    
    new_input = str(pathoutput / "refinement" / f"input.{lig_extension}")
    with open(new_input, "w") as file:
        for idx in orig_idx[best_indexes]:
            file.write(supp.GetItemText(int(idx)))
    args.output_directory = new_output
    args.ligands_sdf = new_input
    equina(args, repeats = known_args.repeats, make_parquet = True, just_gnina = known_args.just_gnina, split_gnina = known_args.split_gnina)
    
    