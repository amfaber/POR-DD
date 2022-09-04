#!/usr/bin/env python
#SBATCH --job-name=LIT-PCBA-benchmark
#SBATCH --cpus-per-task=12
#SBATCH --time=1-00:00:00
#SBATCH --exclude=a00862,a00860
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qzj517@alumni.ku.dk
#SBATCH --output=LIT_PCBA_benchmark.out
#SBATCH --error=LIT_PCBA_benchmark.err
import sys
import os
print(os.getcwd())
print(sys.path)
sys.path.append(".")
import failed_ligands
import rdkit.Chem as Chem
mols = Chem.SmilesMolSupplier("failed_ligands/all_failed.smi")
failed_ligands.generate_for_all_with_time_limit(mols, "failed_ligands", 10)