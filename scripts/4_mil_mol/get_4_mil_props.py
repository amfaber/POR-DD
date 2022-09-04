#%%
import rdkit.Chem as Chem
import pandas as pd
import numpy as np
from pathlib import Path
import os
curdir, folders, files = next(os.walk("/home/qzj517/POR-DD/data/gnina_processed/4_mil_mol/chopped"))
sdfs = [file for file in files if file.endswith(".sdf")]
sorted_sdfs = [sdf.lstrip("0") for sdf in sorted([sdf.zfill(15) for sdf in sdfs])]
dfs = []
for sdf in sorted_sdfs:
    print(f"starting {sdf}")
    molprops = [mol.GetPropsAsDict(includePrivate = True) for mol in Chem.SDMolSupplier(str(Path(curdir) / sdf)) if mol]
    dfs.append(pd.DataFrame(molprops))
full_df = pd.concat(dfs)
full_df.to_parquet("/home/qzj517/POR-DD/data/gnina_processed/4_mil_mol/molprops.parquet")