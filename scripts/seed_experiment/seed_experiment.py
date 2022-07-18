import multiligand_inference 
if __name__ == "__main__": 
    def make_args(seed):
        return ["-l", "/home/qzj517/POR-DD/data/raw_data/cyp_screen/test_3D_opt_1216.sdf",
    "-r", "../data/raw_data/por_structures/3ES9_bioassem_1.pdb", "-o",
    f"/home/qzj517/POR-DD/data/equibind_processed/seed_experiment/seed_{seed}", 
    "--n_workers_data_load", "4", "--batch_size", "32", "--no_skip", "--seed", f"{seed}"]

    for i in range(20):
        multiligand_inference.main(make_args(i))
