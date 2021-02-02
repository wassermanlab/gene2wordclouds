#!/usr/bin/env python

import click
import gzip
import itertools
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
# from scipy.stats import zscore
# import seaborn as sns
import sys

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i", "--input-file",
    help="Input (i.e. entrezid2pmids.json.gz).",
    type=click.File("rt"),
    required=True
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(),
    default="./",
    show_default=True
)

def gene2pmid_distribution(input_file, output_dir):

    # Create output dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Parse input
    with gzip.open(input_file.name, "rt") as handle:
        mappings = json.load(handle)    
    entrezids, pmids, pmids_orthologs = list(map(list, zip(*mappings)))
    input_file.close()

    # Count genes per PMIDs
    data = []
    counts_pmids = __count_genes_per_paper(entrezids, pmids)
    counts_pmids_orthologs = __count_genes_per_paper(entrezids,
        pmids_orthologs)
    counts_pmids_all = __count_genes_per_paper(entrezids,
        list(zip(pmids, pmids_orthologs)))
    for pmid, counts in counts_pmids_all.items():
        booleans = [False, False, True]
        if pmid in counts_pmids:
            booleans[0] = True
        if pmid in counts_pmids_orthologs:
            booleans[1] = True
        data.append([pmid, counts, booleans[0], booleans[1], booleans[2]])
    columns = ["PMID", "Genes", "PMID ∈ gene", "PMID ∈ orthologs",
        "PMID ∈ gene ∪ orthologs"]
    df = pd.DataFrame(data, columns=columns)

    # Compute Z-scores
    mean = np.mean(df["Genes"].tolist())
    std = np.std(df["Genes"].tolist())
    df["Z-score"] = [__Z_score(x, mean, std) for x in df["Genes"].tolist()]

    # Plot
    fig, ax = plt.subplots()
    labels = ["PMID ∈ gene", "PMID ∈ orthologs", "PMID ∈ gene ∪ orthologs"]
    colors = ["#1965B0", "#F7F056", "#DC050C"]
    for l, c in zip(labels, colors):
        x = []; y = []
        for i, v in df[df[l] == True]["Genes"].value_counts().items():
            x.append(i); y.append(v)
        ax.plot(x, y, "o", label=l, color=c)
    x = 2*std+mean
    kwargs = {"ls" : "--", "lw" : 1}
    ax.axvline(x=x, color="black", label="Z-score = 2 (%s)" % round(x, 3),
        **kwargs)
    ax.legend(frameon=False)
    ax.set(xscale="log", yscale="log")
    ax.set_xlabel("Genes")
    ax.set_ylabel("PMIDs")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig("test2.png")

def __count_genes_per_paper(genes, pmids):

    # Initialize
    counts = {}

    for i in range(len(genes)):
        try:
            iterator = set(itertools.chain.from_iterable(pmids[i]))
        except:
            iterator = set(pmids[i])
        for pmid in iterator:
            counts.setdefault(pmid, 0)
            counts[pmid] += 1

    return(counts)

def __Z_score(x, mean, std):
    return((x-mean)/std)

#     cwd = os.getcwd()
#     outdir = os.path.join(cwd, 'genevspmid')
#     if not os.path.exists(outdir):
#         os.mkdir(outdir)
#     os.chdir(outdir)

#     parseelements()
#     countgeneperpmid()
#     getstats()
#     getvisual()

# def parseelements():
#     entrezids = []
#     pmids = []
#     homolopmids = []
#     for line in entrezid2pmids:
#         entrezids.append(line[0])
#         pmids.append(line[1])
#         homolopmids.append(line[2])

#     global flat_pmids
#     flat_pmids = []
#     for e in pmids:
#         flat_pmids += e

#     global flat_homolopmids
#     flat_homolopmids = []
#     for e in homolopmids:
#         flat_homolopmids += e

# def countgeneperpmid():
#     genecount = dict()
#     for pmid in flat_pmids:
#         genecount[pmid] = genecount.get(pmid, 0) + 1

#     df_pmid_genecount = pd.DataFrame.from_dict(genecount, orient="index").rename(columns={0 : "pmid"})
#     df_pmid_count = pd.DataFrame(df_pmid_genecount["pmid"].value_counts())

#     for pmid in flat_homolopmids:
#         genecount[pmid] = genecount.get(pmid, 0) + 1

#     hompmid_genecount = pd.DataFrame.from_dict(genecount, orient="index").rename(columns={0 : "homo+pmid"})
#     df_hompmid_count = pd.DataFrame(hompmid_genecount["homo+pmid"].value_counts())

#     total_genecount = df_pmid_genecount.join(hompmid_genecount)

#     global df_total_sum
#     df_total_sum = df_pmid_count.join(df_hompmid_count).sort_index()

# def getstats():
#     df_total_sum["zscore_pmid"] = (df_total_sum['pmid'] - df_total_sum['pmid'].mean())/df_total_sum['pmid'].std()
#     df_total_sum["zscore_homo+pmid"] = (df_total_sum['homo+pmid'] - df_total_sum['homo+pmid'].mean())/df_total_sum['homo+pmid'].std()
#     df_total_sum.to_csv('numgenes_numpmid_zscore.txt', sep="\t", index=True)

#     mean = [df_total_sum['pmid'].mean(), df_total_sum['homo+pmid'].mean()]
#     median = [df_total_sum['pmid'].median(), df_total_sum['homo+pmid'].median()]
#     stdev = [df_total_sum['pmid'].std(), df_total_sum['homo+pmid'].std()]

#     df_stats = pd.DataFrame(list(zip(mean, median, stdev)),
#                      columns = ['Mean', 'Median', 'Stdev'],
#                      index = ['pmid', 'homo+pmid'])
#     df_stats.to_csv('summarystats.txt', sep="\t", index=True)

# def getvisual():
#     plt.plot(df_total_sum.index.tolist(), df_total_sum['pmid'], label='pmid', color="blue");
#     plt.plot(df_total_sum.index.tolist(), df_total_sum['homo+pmid'], label='homo+pmid', color="red");
#     plt.axvline(x=1, ymin=0, color="lightblue", label="zscore_pmid > 2")
#     plt.axvline(x=2, ymin=0, color="lightpink", label="zscore_homo+pmid > 2")
#     plt.xscale('log')
#     plt.xlabel("# of genes")
#     plt.yscale('log')
#     plt.ylabel("# of pmids")
#     plt.legend()

#     plt.savefig('plot.png')


if __name__ == "__main__":
    gene2pmid_distribution()
