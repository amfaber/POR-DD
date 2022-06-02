import pandas as pd
import os
import re

def load_index(filename = "INDEX_general_PL_data.2020", path = "/home/amfaber/comparison/PDBBind_index"):
    os.chdir(path)
    with open(filename) as file:
        data = {"name": [], "resolution": [], "release_year": [], "aff": [],
                "affinity_type": [], "Kd/Ki": [], "reference": [], "ligand_name": []}
        all_lines = file.readlines()
        delimiters = ["  ", "  ", "  ", "  ", "<=|>=|=|>|<|~", "// ", " ", "\n"]
        for line in all_lines:
            if line.startswith("#"):
                continue
            remaining_line = line

            for delim, ind in zip(delimiters, data):
                output, remaining_line = re.split(delim, remaining_line, maxsplit = 1)
                output = output.strip(" ()")
                data[ind].append(output)
    prefixes = {"m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
    converter = lambda s: float(s[:-2])*prefixes[s[-2]]
    df = pd.DataFrame(data)
    df["Kd/Ki"] = df["Kd/Ki"].apply(converter)
    df = df.astype({"aff": float})
    df["aff"] = -df["aff"]
    return df

def load_equibind(path = "/home/amfaber/gnina_own/results/equibind_dock_real/", skip = []):
    fields = ["name", "affinity1", "affinity2", "RMSD", "CNNscore", "CNNaffinity", "CNNvariance"]
    out = {field : [] for field in fields}
    for name in os.listdir(path):
        if name in skip:
            continue
        with open(os.path.join(path, name, f"{name}.txt")) as file:
            contents = file.read()
            values = re.findall(r"-?\d+\.?\d*", contents)
            if not len(values) == 6:
                continue
            out["name"].append(name)
            for field, val in zip(fields[1:], values):
                out[field].append(val)
    return pd.DataFrame(out).astype({key: val for key, val in zip(out, [str] + [float]*6)})

def load_scoreonly(path = "/home/amfaber/gnina_own/results/equibind_score_only", skip = []):
    fields = ["name", "Affinity", "CNNscore", "CNNaffinity", "CNNvariance"]
    out = {field : [] for field in fields}
    for name in os.listdir(path):
        if name in skip:
            continue
        with open(os.path.join(path, name, f"{name}.txt")) as file:
            contents = file.read()
            values = re.findall(r"-?\d+\.?\d*", contents)
            if not len(values) == 4:
                continue
            out["name"].append(name)
            for field, val in zip(fields[1:], values):
                out[field].append(val)
    return pd.DataFrame(out).astype({key: val for key, val in zip(out, [str] + [float]*6)})

def load_gnina(path = "/home/amfaber/gnina_own/results/gnina_dock/", skip = []):
    fields = ["mode", "affinity", "CNNpose", "CNNaffinity"]
    out = {"name": []}
    out.update({field : [] for field in fields})
    for name in os.listdir(path):
        if name in skip:
            continue
        with open(os.path.join(path, name, f"{name}.txt")) as file:
            contents = file.read()
            end = re.search(r"-+\+-+\+-+\+-+\n", contents).span()[1]
            values = contents[end:].split()
            
            out["name"] += [name]*9
            for i, val in enumerate(values):
                out[fields[i%4]].append(val)
    
    return pd.DataFrame(out).astype({key: typ for key, typ in zip(out, [str, int] + [float]*3)})

def load_bulk(fp):
    out = {}
    with open(fp) as file:
        s = file.read()
    for match in re.findall("[a-zA-Z]+: -?[0-9]+\.?[0-9]*", s):
        key, value = match.split(": ")
        try:
            out[key].append(float(value))
        except KeyError:
            out[key] = [float(value)]
    return pd.DataFrame(out)

def load_bulk2(fp):
    out = {}
    with open(fp) as file:
        s = file.read()
    s = s.split("---BEGIN DATA---\n")[1]
    s = s.split("## ")[1:]
    for mol in s:
        mol = re.split(" |\n", mol, 2)
        idx, name = mol[:2]
        idx = int(idx)
        out[idx] = {"name": name}
        for match in re.findall("[a-zA-Z]+: -?[0-9]+\.?[0-9]*", mol[2]):
            key, val = match.split(": ")
            val = float(val)
            out[idx][key] = val
    return pd.DataFrame(out).T