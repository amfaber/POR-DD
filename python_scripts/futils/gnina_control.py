import os
def gnina_dir(input = "/home/amfaber/equibind/data/results/PDBBind/",
              protein_input = "/home/amfaber/equibind/data/PDBBind_processed/",
              output = "/home/amfaber/gnina_own/results/equibind_dock_real",
              proteinpost = "_protein_processed.pdb",
              mode = "minimize"):
    names = os.listdir(input)
    names = [name for name in names if not name.endswith("failed")]
    dones = os.listdir(output)
    names = [name for name in names if not name in dones]
    for name in names:
        outputdir = os.path.join(output, name)
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        ligpath = os.path.join(input, name, os.listdir(os.path.join(input, name))[0])
        recpath = os.path.join(protein_input, f"{name}/{name}{proteinpost}")
        cmd = f"gnina -l {ligpath} "
        cmd += f"-r {recpath} "
        cmd += f"--{mode} "
        cmd += r"| grep 'Affinity\|RMSD\|CNN' "
        cmd += f"> {outputdir}/{name}.txt"
        #return cmd
        os.system(cmd)