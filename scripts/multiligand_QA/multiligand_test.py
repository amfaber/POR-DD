#%%
from futils import PROT, LIGS, ROOT
import multiligand_inference

multiligand_inference.main([
    "-r", PROT,
    "-l", LIGS,
    "-o", ROOT / "data/testing_clean_regularly"
])
# %%
# from futils import PROT, LIGS
# import rdkit.Chem as Chem
# supp = Chem.MultithreadedSDMolSupplier(LIGS)
# %%
