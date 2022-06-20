from pathlib import Path
from os.path import relpath
_script_dir = Path(__file__).parent
_cwd = Path.cwd()
PROT_ABS = "../../data/raw_data/por_structures/3ES9_1_reduced.pdb"
LIGS_ABS = "../../data/raw_data/cyp_screen/test_3D_opt_1216.sdf"
LIGS_ABS = str((_script_dir / LIGS_ABS).resolve())
PROT_ABS = str((_script_dir / PROT_ABS).resolve())
LIGS = relpath(LIGS_ABS, _cwd)
PROT = relpath(PROT_ABS, _cwd)
