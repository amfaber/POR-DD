#!/usr/bin/env python
import os
import multiligand_inference
from pathlib import Path
import argparse
import sdf_to_parquet
from rdkit import Chem
FLOAT_COLS = ["minimizedAffinity", "minimizedRMSD", "CNNscore", "CNNaffinity", "CNN_VS", "CNNaffinity_variance"]

# def make_supplier(input):
#     if input.split(".")[1] == "smi":
#         return Chem.SmilesMolSupplier(input, titleLine = False)
#     elif input.split(".")[1] == "sdf":
#         return Chem.SDMolSupplier(input)

def do_gnina(r, l, outdir, split = None, make_parquet = False):
    if split is None:
        gnina_cmd = f"gnina -r {r} -l {l} -o {os.path.join(outdir, 'gnina.sdf')} --minimize"
        os.system(gnina_cmd)
        if make_parquet:
            sdf_to_parquet.extract_single_file(str(outdir / 'gnina.sdf'), str(outdir / "props.parquet"),
                float_cols = FLOAT_COLS)

    else:
        # with open("/home/qzj517/POR-DD/Internal_assays/Enamine_3mil_screen/equina_output/problem_gninas.txt") as file:
        #     problems = [int(line) for line in file.readlines()]
        supp = Chem.SDMolSupplier(l)
        os.makedirs(Path(outdir) / "gnina_inps", exist_ok = True)
        os.makedirs(Path(outdir) / "gnina_outs", exist_ok = True)
        os.makedirs(Path(outdir) / "props", exist_ok = True)
        for i in range((len(supp) + split - 1 ) // split):
            end_idx = min((i + 1)*split - 1, len(supp))
            # if end_idx not in problems:
            #     continue
            cur_in = Path(outdir) / "gnina_inps" / f"{end_idx}.sdf"
            cur_out = Path(outdir) / "gnina_outs" / f"{end_idx}.sdf"
            parquet_path = Path(outdir) / "props" / f"{end_idx}.parquet"
            with open(cur_in, "w") as file:
                for j in range(min(split, len(supp) - i*split)):
                    file.write(supp.GetItemText(i*split + j))
            gnina_cmd = f"gnina -r {r} -l {cur_in} -o {cur_out} --minimize"
            os.system(gnina_cmd)
            if make_parquet:
                sdf_to_parquet.extract_single_file(str(cur_out), str(parquet_path),
                    float_cols = FLOAT_COLS)


def equina(multiligand_inference_args, repeats = None, keeps = None, make_parquet = False, just_gnina = False, split_gnina = None):
    args = multiligand_inference_args
    passed_outdir = args.output_directory
    if repeats is None:
        out_dir = Path(passed_outdir)
        if not just_gnina:
            work_was_done = multiligand_inference.main(args = args, keeps = keeps)
        if just_gnina or work_was_done:
            do_gnina(args.rec_pdb, str(out_dir / 'output.sdf'), out_dir, split = split_gnina, make_parquet = make_parquet)
            # gnina_cmd = f"gnina -r {args.rec_pdb} -l {out_dir / 'output.sdf'} -o {out_dir / 'gnina.sdf'} --minimize"
            # os.system(gnina_cmd)
        # if make_parquet:
        #     sdf_to_parquet.extract_single_file(str(out_dir / 'gnina.sdf'), str(out_dir / "props.parquet"), float_cols = FLOAT_COLS)
    else:
        for i in range(repeats):
            out_dir = Path(passed_outdir) / f"rep{i}"
            args.output_directory = out_dir
            args.seed = i
            if not just_gnina:
                work_was_done = multiligand_inference.main(args = args, keeps = keeps)
            if just_gnina or work_was_done:
            # out_dir = Path(out_dir)
                do_gnina(args.rec_pdb, str(out_dir / 'output.sdf'), out_dir, split = split_gnina, make_parquet = make_parquet)
            #     gnina_cmd = f"gnina -r {args.rec_pdb} -l {out_dir / 'output.sdf'} -o {out_dir / 'gnina.sdf'} --minimize --seed {i}"
            #     os.system(gnina_cmd)
            # if make_parquet:
                # sdf_to_parquet.extract_single_file(str(out_dir / 'gnina.sdf'), str(out_dir / "props.parquet"),
                # float_cols = FLOAT_COLS)

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type = int, default = None)
    parser.add_argument("--make_parquet", action="store_true")
    parser.add_argument("--just_gnina", action = "store_true")
    parser.add_argument("--split_gnina", type = int, default = None)

    known_args, rest = parser.parse_known_args()
    args, cmdlineargs = multiligand_inference.parse_arguments(rest)
    args = multiligand_inference.get_default_args(args, cmdlineargs)

    equina(args, repeats = known_args.repeats, make_parquet = known_args.make_parquet, split_gnina = known_args.split_gnina)