#%%
from rdkit import Chem
import matplotlib.pyplot as plt
import sys
import futils
from futils import ROOT
sys.path.append(ROOT / "gflownet/mols")
import gflownet
import gzip
import pickle
import torch
import os
from pathlib import Path
import numpy as np
from futils import mols_to_tiff
from futils import equina
from futils import loaders


#%%
def inference(input, output, sample_random, number, device = None):
    if device is None:
        device = "cuda:0" if torch.cuda.is_available() else "cpu"
    input = Path(input)
    os.makedirs(output, exist_ok = True)

    with gzip.open(input / "_0/info.pkl.gz", "rb") as file:
        args = pickle.load(file)["args"]


    dataset = gflownet.Dataset(args, ROOT.ds / "gflownet/mols/data/blocks_PDB_105.json", device)

    model = gflownet.make_model(args, dataset.mdp)
    model.to(device)

    if sample_random:
        dataset.set_sampling_model(model, None)
        sampled_mols = [dataset._get_sample_model()[-1][-2].mol for i in range(number)]
        with Chem.SmilesWriter(str(output / "sampled_from_random.smi"), includeHeader = False) as w:
            [w.write(mol) for mol in sampled_mols]


    with gzip.open(str(input / "_0/params.pkl.gz")) as file:
        params = pickle.load(file)

    currparams = model.parameters()
    for old, new in zip(currparams, params):
        old.data = torch.tensor(new, dtype = dataset.mdp.floatX, device = device)

    dataset.set_sampling_model(model, None)

    sampled_mols = [dataset._get_sample_model()[-1][-2].mol for i in range(number)]
    with Chem.SmilesWriter(str(output / "sampled_from_trained.smi"), includeHeader = False) as w:
        [w.write(mol) for mol in sampled_mols]


def generate_pictures(sample_path, figure_path, get_random = True):
    if get_random:
        random = [mol for mol in Chem.SmilesMolSupplier(f"{sample_path}/sampled_from_random.smi")]
    trained = [mol for mol in Chem.SmilesMolSupplier(f"{sample_path}/sampled_from_trained.smi")]
    square = 4
    mols_to_tiff.mols_to_tiff(f"{figure_path}/random.tiff", random, square)
    mols_to_tiff.mols_to_tiff(f"{figure_path}/trained.tiff", trained, square)

def score_samples(sample_input, por_struture, get_random = True):
    sample_input = Path(sample_input)
    if get_random:
        equina.equina(sample_input / "sampled_from_random.smi", por_struture,
                            sample_input / "random_affinity", batch_size = 32)
    equina.equina(sample_input / "sampled_from_trained.smi", por_struture,
                        sample_input / "trained_affinity", batch_size = 32)

#%%
def plot_sample_distribution(sample_path, figure_path, measure = "CNN_VS"):
    sample_path = Path(sample_path)
    plt.rc("axes", labelsize = 20)
    plt.rc("xtick", labelsize = 13)
    plt.rc("ytick", labelsize = 13)
    # four_mil_data = pd.read_parquet(ROOT / "data/gnina_processed/4_mil_mol/molprops.parquet")
    trained_mols = list(Chem.SDMolSupplier(str(sample_path / "trained_affinity" / "gnina.sdf")))
    random_mols = list(Chem.SDMolSupplier(str(sample_path / "random_affinity" / "gnina.sdf")))
    
    trained_df = loaders.mols_to_df(trained_mols)
    random_df = loaders.mols_to_df(random_mols)

    plt.hist(trained_df[measure], histtype = "step", density = True, bins = 20)
    plt.hist(random_df[measure], histtype = "step", density = True, bins = 20)
    if measure == "CNN_VS":
        xlabel = "CNN VS from EquiBind + gnina"
    elif measure == "CNNaffinity":
        xlabel = "CNN affinity from EquiBind + gnina"
    else:
        xlabel = None
    plt.xlabel(xlabel)
    plt.ylabel("Probability density")
    # plt.title("3.3 million Zinc samples vs 1000 gflownet samples")
    plt.title("Trained vs untrained gflownet")
    plt.legend(["Gflownet samples trained", "Gflownet samples random"])
    plt.savefig(figure_path / f"trained_vs_random_{measure}.png", bbox_inches = "tight")
    plt.close()

def plot_training(input_path, figure_path):
    input_path = Path(input_path)
    figure_path = Path(figure_path)
    with gzip.open(input_path / "_0" / "info.pkl.gz") as file:
        train_losses = pickle.load(file)["train_losses"]
    total, terminal, flow, _, _ = list(zip(*train_losses))
    fig, axes = plt.subplots(1, 3, figsize = np.array([6*3, 4.8]))
    mv_avg = 20
    axes[0].plot(np.convolve(total, np.ones(mv_avg), "valid")/mv_avg)
    axes[0].set_title("Total loss")
    axes[1].plot(np.convolve(terminal, np.ones(mv_avg), "valid")/mv_avg)
    axes[1].set_title("Terminal loss")
    axes[2].plot(np.convolve(flow, np.ones(mv_avg), "valid")/mv_avg)
    axes[2].set_title("Flow loss")
    plt.savefig(figure_path / "training_losses.png", bbox_inches = "tight")
    plt.close()

#%%
#%%

def characterize(input, results_path = None, device = None, por_structure = None):
    
    if results_path is None:
        results_path = ROOT / "gflownet/mols/results"
    else:
        results_path = Path(results_path)
    
    if por_structure is None:
        por_structure = ROOT / "data/raw_data/por_structures/3QE2_1_reduced.pdb"
    
    characterization_path = results_path / input / "characterization"
    sample_path = characterization_path / "samples"
    figure_path = characterization_path / "figures"
    os.makedirs(sample_path, exist_ok = True)
    os.makedirs(figure_path, exist_ok = True)

    plot_training(results_path / input, figure_path)
    
    inference(results_path / input, sample_path,
        sample_random = True, number = 1000, device = device
        )
    
    generate_pictures(sample_path, figure_path)
    score_samples(sample_path, por_structure)
    plot_sample_distribution(sample_path, figure_path, "CNN_VS")
    plot_sample_distribution(sample_path, figure_path, "CNNaffinity")

