import pandas as pd
import os
import re
from futils import ROOT
from rdkit import Chem

def load_index(filename = "INDEX_general_PL_data.2020", path = ROOT / "comparison/PDBBind_index"):
    os.chdir(path)
    with open(filename) as file:
        data = {"name": [], "resolution": [], "release_year": [], "aff": [],
                "affinity_type": [], "Kd/Ki": [], "reference": [], "ligand_name": []}
        all_lines = file.readlines()
        delimiters = ["  ", "  ", "  ", "  ", "<=|>=|=|>|<|~", "// ", " ", "\n"]
        for line in all_lines:
            if line.startswith("#"):
                continue
            remaining_line = line

            for delim, ind in zip(delimiters, data):
                output, remaining_line = re.split(delim, remaining_line, maxsplit = 1)
                output = output.strip(" ()")
                data[ind].append(output)
    prefixes = {"m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
    converter = lambda s: float(s[:-2])*prefixes[s[-2]]
    df = pd.DataFrame(data)
    df["Kd/Ki"] = df["Kd/Ki"].apply(converter)
    df = df.astype({"aff": float})
    df["aff"] = -df["aff"]
    return df

def load_bulk2(fp):
    out = {}
    with open(fp) as file:
        s = file.read()
    s = s.split("---BEGIN DATA---\n")[1]
    s = s.split("## ")[1:]
    for mol in s:
        mol = re.split(" |\n", mol, 2)
        idx, name = mol[:2]
        idx = int(idx)
        out[idx] = {"name": name}
        for match in re.findall("[a-zA-Z]+: -?[0-9]+\.?[0-9]*", mol[2]):
            key, val = match.split(": ")
            val = float(val)
            out[idx][key] = val
    return pd.DataFrame(out).T

def load_mols_across_seeds(dirpath):
    from rdkit.Chem import SDMolSupplier
    join = os.path.join
    subfolders = next(os.walk(dirpath))[1]
    output = {}
    for folder in subfolders:
        supp = SDMolSupplier(join(dirpath, folder, "output.sdf"))
        mols = [mol for mol in supp]
        output[folder] = mols
    return output

def load_gnina_across_seeds(dirpath):
    join = os.path.join
    files = next(os.walk(dirpath))[2]
    output = {}
    for file in files:
        output[file] = load_bulk2(join(dirpath, file))
    return output

def mols_to_df(mollist):
    dicts = [mol.GetPropsAsDict(includePrivate = True) for mol in mollist]
    return pd.DataFrame(dicts)