#!/usr/bin/env python
#%%

import os
from os.path import join

def gnina_dir(inputdir, outputdir, prot):
    dir, _, files = next(os.walk(inputdir))
    files = [file for file in files if file.endswith(".sdf")]
    for file in files:
        os.system(f"gnina -l {join(dir, file)} -r {prot} -o {join(outputdir, file)} --minimize > {join(outputdir, file.split('.')[0])}.txt")

def chop(inputfile, outputdir, size = 50*1048576):
    with open(inputfile, "rb") as file:
        index = -1
        while True:
            chunk = file.read(size)
            if chunk == b'':
                break
            body, leftover = chunk.rsplit(b"$$$$\n", 1)
            file.seek(-len(leftover), 1)
            n = body.count(b"$$$$\n") + 1
            index += n
            name = join(outputdir, f"{index}.sdf")
            with open(name, "wb") as f:
                f.write(body)

def gather(inputdir, outputfile):
    with open(outputfile, "wb") as file:
        for file in os.listdir(inputdir):
            with open(join(inputdir, file), "rb") as f:
                file.write(f.read())
                # file.write(b"$$$$\n")

def gather_shell(inputdir, outputfile, remove = True):
    if remove:
        os.system(f"rm {outputfile}")
    os.system(f"ls {join(inputdir, '*.sdf')} -v | xargs cat >> {outputfile}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Chop SDF file into smaller SDF files")
    subparsers = parser.add_subparsers(dest="command")
    chop_parser = subparsers.add_parser("chop", help="Chop SDF file into smaller SDF files")
    chop_parser = subparsers.add_parser("chop")
    chop_parser.add_argument("inputfile", help="SDF file to chop")
    chop_parser.add_argument("outputdir", help="Directory to save chopped SDF files")
    chop_parser.add_argument("--size", help="Size of each chunk", type=int, default=50*1048576)
    
    gnina_dir_parser = subparsers.add_parser("gnina_dir")
    gnina_dir_parser.add_argument("inputdir", help="Directory of SDF files")
    gnina_dir_parser.add_argument("outputdir", help="Directory to write gnina output to")
    gnina_dir_parser.add_argument("prot", help="Protein path")

    gather_parser = subparsers.add_parser("gather")
    gather_parser.add_argument("inputdir", help="Directory of SDF files")
    gather_parser.add_argument("outputfile", help="File to write SDF to")

    gather_shell_parser = subparsers.add_parser("gather_shell")
    gather_shell_parser.add_argument("inputdir", help="Directory of SDF files")
    gather_shell_parser.add_argument("outputfile", help="File to write SDF to")
    gather_shell_parser.add_argument("--remove", help="Remove output file before writing", action="store_true")

    args = parser.parse_args()

    locals()[args.command](**vars(args))

    # if args.command == "chop":
    #     chop(args.inputfile, args.outputdir, args.size)
    # elif args.command == "gnina_dir":
    #     gnina_dir(args.inputdir, args.outputdir, args.prot)
    # elif args.command == "gather":
    #     gather(args.inputdir, args.outputfile)
    # elif args.command == "gather_shell":
    #     gather_shell(args.inputdir, args.outputfile)