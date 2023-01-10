#!/usr/bin/env python
from rdkit.Chem import Draw
import numpy as np
def mols_to_tiff(save_location, mols, square_side, clearconfs = False, legends = None, highlight_mols = None):
    if clearconfs:
        if __name__ != "__main__":
            print("Removing conformation of every mol in provided list IN PLACE.")
        [mol.RemoveConformer(0) for mol in mols]
    pngs = []
    for i in range(int(np.ceil(len(mols)/square_side**2))):
        cur_mols = mols[i*square_side**2: (i+1)*square_side**2]
        cur_idx = range(i*square_side**2, min((i+1)*square_side**2, len(mols)))
        if legends is None:
            legend_getter = str
        else:
            legend_getter = lambda i: str(legends[i])
        if highlight_mols is None:
            highlighters = {}
        else:
            highlighters = dict(
                highlightAtomLists = [list(range(mol.GetNumAtoms())) if idx in highlight_mols else [] for idx, mol in zip(cur_idx, cur_mols)],
                highlightBondLists = [list(range(mol.GetNumBonds())) if idx in highlight_mols else [] for idx, mol in zip(cur_idx, cur_mols)],
            )
        png = Draw.MolsToGridImage(cur_mols,
        subImgSize = (400, 400),
        molsPerRow = square_side,
        returnPNG = False,
        legends = list(map(legend_getter, cur_idx)),
        **highlighters
        )
        pngs.append(png)
    pngs[0].save(save_location, resolution = 200, save_all = True, append_images = pngs[1:], quality = 100)

if __name__ == "__main__":
    from rdkit.Chem import SDMolSupplier
    import argparse
    from pathlib import Path
    import pandas as pd
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--sort_by", default = "CNNaffinity")
    p.add_argument("-n", type = int, default = 100)
    p.add_argument("--make_parquet", action = "store_true")
    p.add_argument("--parquet_path", default = None)
    p.add_argument("--name_col", default = None)
    args = p.parse_args()
    input = Path(args.input)
    if input.suffix != ".sdf":
        input = input / "gnina.sdf"
    if args.make_parquet:
        import sdf_to_parquet
        sdf_to_parquet.extract_single_file(str(input), args.parquet_path)
    
    if args.parquet_path is not None:
        parquet_path = Path(args.parquet_path)
    else:
        parquet_path = input.parent / "props.parquet"
    
    supp = SDMolSupplier(str(input))
    if parquet_path.exists():
        df = pd.read_parquet(parquet_path)
        df.sort_values(args.sort_by, inplace = True, ascending = False)
        # with open(input.parent / "success.txt") as file:
        #     orig_idx = np.array([int(line.split(" ")[0]) for line in file.readlines()])
        # best_idx = orig_idx[df.index[:args.n]]
        mols = [supp[int(idx)] for idx in df.index[:args.n]]
    else:
        print("Couldn't find parquet. Using SDF order")
        mols = [mol for _, mol in zip(range(args.n), supp)]
    
    if args.sort_by is not None:
        if args.name_col is not None:
            legends = list(zip(df[args.sort_by], df[args.name_col]))
        else:
            legends = list(df[args.sort_by])

    else:
        legends = None
    mols_to_tiff(str(input.parent / "mol_images.pdf"), mols, 5, clearconfs = True, legends = legends)

