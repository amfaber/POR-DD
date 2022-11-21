import multiligand_inference
from pathlib import Path
import os
import sys
def equina(ligands, receptor, output, additional_arglist = [], **kwargs):
    additional_arglist = list(map(str, additional_arglist))
    ligands = str(ligands)
    receptor = str(receptor)
    output = str(output)
    args, cmd = multiligand_inference.parse_arguments(
        ["-r", receptor, "-l", ligands, "-o", output] + additional_arglist
        )
    args = multiligand_inference.get_default_args(args, cmd)
    for key in kwargs:
        if key not in args:
            raise KeyError("Unknown argument")
    args.__dict__.update(kwargs)
    output = Path(output)
    multiligand_inference.main(args = args)
    gnina_cmd = f"gnina -l {output / 'output.sdf'} -r {receptor} -o {output / 'gnina.sdf'} \
--minimize"
    os.system(gnina_cmd)

    