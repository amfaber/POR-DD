#%%
import sys
import os
from os.path import join
import multiligand_inference
import rdkit.Chem as Chem
sys.path.append("../EquiBind/datasets")
from multiple_ligands import Ligands

def experiment(inputdir, equibind_output, gnina_output):
    walker = os.walk(inputdir)
    next(walker)
    
    args, cmdline = multiligand_inference.parse_arguments([
        "--batch_size", "32",
        # "--no_skip"
    ])
    args = multiligand_inference.get_default_args(args, cmdline)
    model = multiligand_inference.load_model(args)

    for curdir, folders, files in walker:
        prots = [join(curdir, file) for file in files if file.endswith(".pdb")]
        ligs = [join(curdir, file) for file in files if file.endswith(".smi")]
        for lig in ligs:

            dataset = Ligands(lig, None, args)
            n_dataloaders = 12 if len(dataset) > 1000 else 0
            args.n_workers_data_load = n_dataloaders

            for prot in prots:
                ligname = os.path.splitext(os.path.basename(lig))[0]
                protname = os.path.splitext(os.path.basename(prot))[0]
                outputdir = join(
                    "LIT-PCBA",
                    os.path.basename(curdir),
                    ligname,
                    protname
                    )
                equibind_lig = join(equibind_output, outputdir)

                args.output_directory = equibind_lig
                args.ligands_sdf = lig
                args.rec_pdb = prot

                dataset.skips = multiligand_inference.find_previous_work(args)
                
                dataset.rec_graph = multiligand_inference.load_rec(args)


                work_was_done = multiligand_inference.main(args = args, lig_dataset = dataset, model = model)
                if work_was_done:
                    equibind_lig = join(equibind_lig, "output.sdf")
                    os.makedirs(join(gnina_output, outputdir), exist_ok = True)
                    gninacmd = f"gnina -r {prot} -l {equibind_lig} "\
                        f"-o {join(gnina_output, outputdir, 'gnina_out.sdf')} "\
                        "--minimize "\
                        f"> {join(gnina_output, outputdir, 'gnina_out.txt')}"
                    os.system(gninacmd)
                else:
                    print(f"No work was done for LIG = {ligname}, PROT = {protname}\nSkipping.")
                # break
            # break
        # break

#%%
if __name__ == "__main__":
    x = experiment("../data/raw_data/LIT-PCBA/AVE_unbiased", "../data/equibind_processed", "../data/gnina_processed")
    print(x)

# %%
