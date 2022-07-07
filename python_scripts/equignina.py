#!/usr/bin/env python
import argparse
import os
from models import equibind
import multiligand_inference
from futils import ROOT
from pathlib import Path

args, cmdlineargs = multiligand_inference.parse_arguments()
args = multiligand_inference.get_default_args(args, cmdlineargs)
out_dir = args.output_directory
equibind_out_dir = ROOT / Path("data/equibind_processed") / out_dir
args.output_directory = str(equibind_out_dir)
multiligand_inference.main(args = args)
gnina_out_dir = ROOT / Path("data/gnina_processed") / out_dir
os.makedirs(gnina_out_dir, exist_ok=True)
gnina_cmd = f"gnina -r {args.rec_pdb} -l {equibind_out_dir / 'output.sdf'} -o {gnina_out_dir / 'gnina.sdf'} --minimize > {gnina_out_dir / 'gnina.txt'}"
os.system(gnina_cmd)