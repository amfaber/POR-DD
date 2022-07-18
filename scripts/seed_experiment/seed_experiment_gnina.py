import os
from re import sub
def gnina_seed_dir(dirpath):
    rec_path = "/home/qzj517/POR-DD/data/raw_data/cyp_screen/protein.pdb"
    output = "/home/qzj517/POR-DD/data/gnina_processed/seed_experiment"
    
    join = os.path.join
    subfolders = next(os.walk(dirpath))[1]
    for folder in subfolders:
        lig_path = join(dirpath, folder, "output.sdf")
        output_path = join(output, folder + ".txt")
        cmd = f"gnina -r {rec_path} -l {lig_path} --minimize > {output_path}"
        os.system(cmd)

if __name__ == "__main__":
    gnina_seed_dir("/home/qzj517/POR-DD/data/equibind_processed/seed_experiment")