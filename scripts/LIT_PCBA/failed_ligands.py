import os
from os.path import join
import futils
import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem import AllChem

def get_failed_entries(path):
    entries = []
    for curdir, folders, files in os.walk(path):
        if not "failed.txt" in files:
            continue
        keys = os.path.relpath(curdir, path).split(os.path.sep)
        with open(join(curdir, "failed.txt")) as f:
            
            data = tuple(tuple(line.strip().split(" ")) for line in f.readlines())
            # print(data)
            entries.append((*keys, data))
    df = pd.DataFrame(entries)
    df.columns = ["target", "ligand", "protein", "failed"]
    return df

def check_if_all_proteins_fail_on_identical_ligands(df):
    amount_of_failed_per = df.groupby(["target", "ligand"])["failed"].aggregate(lambda seq: len(set(seq)))
    return (amount_of_failed_per == 1).all()

def get_first_protein(df):
    return df.groupby(["target", "ligand"])["failed"].aggregate(lambda series: series.iloc[0])

def put_all_failed_into_one_file(df, input_dir, output_dir):
    all_mols = []
    for (target, ligand), index_name_pairs in df.iteritems():
        # print(f"target = {target}, ligand = {ligand}")
        path = join(input_dir, target, ligand) + ".smi"
        supp = Chem.SmilesMolSupplier(path, sanitize = False, titleLine = False)
        # print(f"supp len = {len(supp)}")
        for idx, name, in index_name_pairs:
            idx = int(idx)
            # print(f"idx = {idx}")
            if idx == len(supp) and name == "None":
                continue
            mol = supp[idx]
            if not mol.GetProp("_Name") == name:
                print(f"Target = {target}, ligand = {ligand}, idx = {idx}, name = {name}")
                raise ValueError(f"{mol.GetProp('_Name')} != {name}")
            mol.SetProp("origin", join(target, ligand))
            all_mols.append(mol)
    writer = Chem.SmilesWriter(join(output_dir, "all_failed.smi"))
    writer.SetProps(list(all_mols[0].GetPropNames()))
    for mol in all_mols:
        writer.write(mol)
    return all_mols

import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def generate_for_all_with_time_limit(mols, output_dir, timeout):
    mols = [AllChem.AddHs(mol) for mol in mols]
    ps = AllChem.ETKDGv2()
    for i, mol in enumerate(mols):
        with open(join(output_dir, "embed_attempts.txt"), "a") as f:
            try:
                with time_limit(timeout):
                    f.write(f"{i} {mol.GetProp('_Name')} {AllChem.EmbedMolecule(mol, ps)}\n")
            except TimeoutException:
                f.write(f"{i} {mol.GetProp('_Name')} Timeout\n")