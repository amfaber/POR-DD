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