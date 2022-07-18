#%%
import rdkit.Chem as Chem
from os.path import join

def chop(inputfile, outputdir, size = 50*1048576):
    with open(inputfile, "rb") as file:
        index = -1
        while True:
            chunk = file.read(size)
            if chunk == b'':
                break
            body, leftover = chunk.rsplit(b"$$$$\n", 1)
            file.seek(-len(leftover), 1)
            n = body.count(b"$$$$\n")
            index += n
            name = join(outputdir, f"{index}.sdf")
            with open(name, "wb") as f:
                f.write(body)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Chop SDF file into smaller SDF files")
    parser.add_argument("inputfile", help="SDF file to chop")
    parser.add_argument("outputdir", help="Directory to save chopped SDF files")
    parser.add_argument("--size", help="Size of each chunk", type=int, default=50*1048576)
    args = parser.parse_args()
    chop(args.inputfile, args.outputdir, args.size)