# %%
import rdkit.Chem as Chem
import os
from os.path import join
import shutil
import yaml
import tempfile
import numpy as np
import pandas as pd
import sys
from futils import ROOT


# %%
ligfn = ROOT / "data/raw_data/cyp_screen/test_3D_opt_1216.sdf"
supp = Chem.SDMolSupplier(ligfn)
first_n = 10
protfn = ROOT / "data/raw_data/por_structures/3ES9_bioassem_1.pdb"

top_dir = ROOT / "data/raw_data/multiligand_QA"

multi_in_dir = join(top_dir, f"multi_input")
multi_in = join(multi_in_dir, f"first_{first_n}.sdf")
multi_outs = [join(top_dir, "multi_output_0"), join(top_dir, "multi_output_1")]
multi_path = ROOT / "EquiBind/multiligand_inference.py"

single_location = ROOT / "EquiBind/" #THIS WILL BE cd'ed INTO BEFORE RUNNING THE "SINGLE" SCRIPT AND OUT OF AFTERWARDS.
path_to_top_dir_from_single_location = ROOT / "data/raw_data/multiligand_QA"
single_outs = [join(path_to_top_dir_from_single_location, "single_output_0"), join(path_to_top_dir_from_single_location, "single_output_1")]
single_in = join(path_to_top_dir_from_single_location, "single_input")
single_path = "inference.py"

save_location = join(top_dir, "comparison.csv")

cleanup = False

suppress_runs = False

pd.options.display.float_format = "{:.2g}".format

batch_size = 5


# %%
os.makedirs(multi_in_dir, exist_ok = True)
[os.makedirs(multi_out, exist_ok = True) for multi_out in multi_outs]
os.makedirs(single_in, exist_ok = True)
[os.makedirs(single_out, exist_ok = True) for single_out in single_outs]

with open(multi_in, "w+") as multi_file:
    for i in range(first_n):
        lig = supp[i]
        name = lig.GetProp("_Name").replace(" ", "_")
        single_dir = join(single_in, name)
        os.makedirs(single_dir, exist_ok = True)
        shutil.copy2(protfn, join(single_dir, "protein.pdb"))
        # with Chem.SDWriter(join(single_dir, "ligand.sdf")) as w:
        #     w.write(lig)
        with open(join(single_dir, "ligand.sdf"), "w+") as w:
            w.write(supp.GetItemText(i))
        multi_file.write(supp.GetItemText(i))

# %%
# with HiddenPrints(suppress_runs):
suppressor = " > /dev/null" if suppress_runs else ""
for i, multi_out in enumerate(multi_outs):
    multi_cmd = f"python {multi_path} -o {multi_out} -r {protfn} -l {multi_in}"\
    f" --no_skip --batch_size 3 --n_workers_data_load 0 --seed {i} --batch_size {batch_size}" + suppressor
    os.system(multi_cmd)

for i, single_out in enumerate(single_outs):
    with tempfile.NamedTemporaryFile("w+") as tmp:
        yamldict = {'run_dirs': ['flexible_self_docking'],
                    'inference_path': single_in,
                    'test_names': 'timesplit_test',
                    'output_directory': single_out,
                    'run_corrections': True,
                    'use_rdkit_coords': False,
                    'save_trajectories': False,
                    'seed': i,
                    'num_confs': 1}
        yaml.dump(yamldict, tmp)
        tmp.seek(0)
        remember_path = os.getcwd()
        os.chdir(single_location)
        single_cmd = f"python {single_path} --config={tmp.name}" + suppressor
        os.system(single_cmd)
        os.chdir(remember_path)

# %%
single_ligs = []
for single_out in single_outs:
    walker = os.walk(single_out)
    next(walker)
    single_ligs.append([])
    for path, folders, files in walker:
        lig = Chem.SDMolSupplier(join(path, files[0]), removeHs = False, sanitize = False)[0]
        single_ligs[-1].append(lig)



# %%
multi_ligs = [[lig for lig in Chem.SDMolSupplier(join(m_outs, "output.sdf"), removeHs = False, sanitize = False)] for m_outs in multi_outs]

# %%
df = pd.DataFrame({
"multi_ligs_0": {lig.GetProp("_Name"): lig for lig in multi_ligs[0]},
"multi_ligs_1": {lig.GetProp("_Name"): lig for lig in multi_ligs[1]},
"single_ligs_0": {lig.GetProp("_Name"): lig for lig in single_ligs[0]},
"single_ligs_1": {lig.GetProp("_Name"): lig for lig in single_ligs[1]},
})

# %%
if cleanup:
    shutil.rmtree(multi_in_dir, ignore_errors = True)
    [shutil.rmtree(multi_out, ignore_errors = True) for multi_out in multi_outs]
    shutil.rmtree(single_in, ignore_errors = True)
    [shutil.rmtree(single_out, ignore_errors = True) for single_out in single_outs]

# %%
coords = df.applymap(lambda lig: lig.GetConformer().GetPositions())
n_atoms = df.iloc[:, 0].apply(lambda lig: lig.GetNumAtoms())

# %%
func = lambda comp1, comp2: [np.sqrt(((coords[comp1] - coords[comp2])**2).apply(np.sum))/n_atoms]
series = func("single_ligs_0", "multi_ligs_0") + func("multi_ligs_1", "multi_ligs_0") + func("single_ligs_0", "single_ligs_1")
# series = [np.sqrt(((coords["single_ligs_0"] - coords["multi_ligs_0"])**2).apply(np.sum))/n_atoms.iloc[:, 0]]
# series += [np.sqrt(((coords["multi_ligs_1"] - coords["multi_ligs_0"])**2).apply(np.sum))/n_atoms.iloc[:, 0]]
colnames = ["single_vs_multi_0", "multi_1_vs_multi_0", "single_0_vs_single_1"]

# %%
output = pd.DataFrame({col: s for col, s in zip(colnames, series)})
output.to_csv(save_location)
print(output)

# %%



