equibind_inference.py -o ../data/equibind_processed/chainA -r ../data/raw_data/por_structures/3ES9_bioassem_1.pdb -l ../data/raw_data/cyp_screen/test_3D_opt_1216.sdf
gnina -o ../data/gnina_processed/chainA/chainA.sdf -r ../data/raw_data/por_structures/3ES9_bioassem_1.pdb -l ../data/equibind_processed/chainA/output.sdf --minimize > ../data/gnina_processed/chainA/chainA.txt

gnina -o ../data/gnina_processed/chainA_full_gnina/chainA_full_gnina.sdf -r ../data/raw_data/por_structures/3ES9_bioassem_1.pdb\
 -l ../data/equibind_processed/chainA/output.sdf --autobox_ligand ../data/equibind_processed/chainA/output.sdf\
  > ../data/gnina_processed/chainA_full_gnina/chainA_full_gnina.txt