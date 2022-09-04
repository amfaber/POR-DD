import pandas as pd

def get_subset(df, dir):
    return df[df["directory"] == dir]
def get_data(df, n_df):
    dirs = df["directory"].unique()
    output_data = {}
    for dir in dirs:
        subset = get_subset(df, dir)
        one_percent = int(len(subset)*0.01)+1
        top1 = subset.iloc[:one_percent]
        n_true_actives, n_true_total = n_df.loc[dir]
        normalizer = (min(one_percent, n_true_actives) / one_percent) / (n_true_actives / n_true_total)
        ef1 = (top1["active_groundtruth"].sum()/len(top1)) / (n_true_actives / n_true_total)

        output_data[dir] = dict(EF1 = ef1,
                                NEF1 = ef1/normalizer,
                                n_actives_found = top1["active_groundtruth"].sum(),
                                length = one_percent,
                                gt_active_proportion = (n_true_actives / n_true_total),
                                total_succeeded = len(subset),
                                total_ligands = n_true_total
                                )
    out_df = pd.DataFrame(output_data).T.astype({"length": int, "n_actives_found": int, "total_succeeded": int, "total_ligands": int})
    out_df.rename(columns = {"length": r"1% is equal to"}, inplace = True)
    out_df.sort_index(inplace = True)
    out_df.index.name = "Target"
    return out_df
