#%%
import sys
import os
from os.path import join
import multiligand_inference
import rdkit.Chem as Chem
sys.path.append("../EquiBind/datasets")
from multiple_ligands import Ligands
import argparse

#%%
def experiment(inputdir, equibind_output, gnina_output, targets_to_run = None):
    walker = os.walk(inputdir)
    next(walker)
    
    args, cmdline = multiligand_inference.parse_arguments([
        "--batch_size", "32",
        # "--no_skip"
    ])
    args = multiligand_inference.get_default_args(args, cmdline)
    model = multiligand_inference.load_model(args)

    for curdir, folders, files in walker:
        if targets_to_run is not None and curdir.rsplit("/", 1)[1] not in targets_to_run:
            continue
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
    from futils import ROOT
    parser = argparse.ArgumentParser(description="Run LIT-PCBA experiment.")
    parser.add_argument("-i", "--inputdir", default = ROOT / "data/raw_data/LIT-PCBA/AVE_unbiased", help="Directory containing input files.")
    parser.add_argument("-eo", "--equibind_output", default = ROOT / "data/equibind_processed", help="Directory to write equibind output to.")
    parser.add_argument("-go", "--gnina_output", default = ROOT / "data/gnina_processed", help="Directory to write gnina output to.")
    parser.add_argument("--targets", default = None, nargs="*", help="Targets to run.")
    args = parser.parse_args()
    # print(args)
    experiment(args.inputdir, args.equibind_output, args.gnina_output, args.targets)

# %%
