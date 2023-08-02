import sys
import pandas as pd
import os


def get_rRNAs(cds_file):
    rRNA_genes = []
    with open(cds_file, "r") as cds:
        cds.readline()
        cds.readline()
        for l in cds:
            space_break = l.split("\t")
            if len(space_break) > 3 and space_break[2] in ["rRNA"]:
                colon_break = space_break[8].split(";")
                rRNA_genes.extend(
                    i[i.index("=") + 1 :] for i in colon_break if i.startswith("gene=")
                )
    return rRNA_genes


def count_percent_rRNAs(rRNA_genes, dataframe_file):
    experiment_df = pd.read_csv(dataframe_file)
    num_rRNA = sum(experiment_df["Counts"].loc[experiment_df["Name"].isin(rRNA_genes)])
    num_mRNA = sum(experiment_df["Counts"].loc[~experiment_df["Name"].isin(rRNA_genes)])
    return num_rRNA, num_mRNA


if __name__ == "__main__":
    """
    Program starts here when called from command line.
    Expected arguments (in order):
        -cds_file = which cds to include
        -data_folder = string path where datafiles to analyze are kept
        -data_file_end = string with end of df file. (eg 'dataframe.txt')
    """

    cds_file = sys.argv[1]
    data_folder = sys.argv[2]
    data_file_end = sys.argv[3]
    rRNA_genes = get_rRNAs(cds_file)
    print(rRNA_genes)
    rRNAs_in_exp = {}
    for f in os.listdir(data_folder):
        if f.endswith(data_file_end):
            rRNA, mRNA = count_percent_rRNAs(rRNA_genes, "".join([data_folder, f]))
            name = f[: f.rfind(data_file_end)]
            rRNAs_in_exp[name] = {"mRNA": mRNA, "rRNA": rRNA}

    rRNA_file_name = "".join([data_folder, "breakdown_of_rnas.txt"])
    with open(rRNA_file_name, "w") as f:
        for exp, rna_breakdown in rRNAs_in_exp.items():
            f.write(f"Experiment: {exp}\n")
            f.write(
                "\n".join(
                    [
                        f"{rna_type}:  {rna_count}"
                        for rna_type, rna_count in rna_breakdown.items()
                    ]
                )
            )
            f.write("\n")
            f.write(
                f'Percent rRNA: {100*rna_breakdown["rRNA"]/(rna_breakdown["rRNA"] + rna_breakdown["mRNA"]):.2f}'
            )
            f.write("\n\n")
