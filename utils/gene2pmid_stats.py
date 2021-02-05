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
#import sys

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
@click.option(
    "--zscore",
    help="Filter pmids with |zscore| higher than threshold.",
    type=float,
    default=None
)

def cli(**params):
    # Create output dir
    if not os.path.isdir(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Parse input
    with gzip.open(params["input_file"].name, "rt") as handle:
        mappings = json.load(handle)    
    entrezids, pmids, pmids_orthologs = list(map(list, zip(*mappings)))
    params["input_file"].close()    

    __get_genes4pmids_stats(entrezids, pmids, pmids_orthologs, params["output_dir"])   
    __get_pmids4genes_stats(entrezids, pmids, pmids_orthologs, params["output_dir"]) 

def __get_genes4pmids_stats(entrezids, pmids, pmids_orthologs, outdir):
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
    
    df.to_csv(os.path.join(outdir, "genesperpmids_table.tsv.gz"), sep="\t", index=False)
    
    # Plot
    __get_genes4pmids_plot(df, mean, std, outdir)
    
    return(df)


def __get_pmids4genes_stats(entrezids, pmids, pmids_orthologs, outdir):
    columns = ["Genes", "PMIDs", "Ortholog PMIDs", "All PMIDs"]
    
    data = []
    for i in range(len(entrezids)):
        data.append([entrezids[i], len(pmids[i]), len(pmids_orthologs[i]), len(pmids[i]) + len(pmids_orthologs[i])])
    df = pd.DataFrame(data, columns=columns)
    
    # Compute Z-scores
    mean = np.mean(df["All PMIDs"].tolist())
    std = np.std(df["All PMIDs"].tolist())
    df["Z-score"] = [__Z_score(x, mean, std) for x in df["All PMIDs"].tolist()]
    
    df.to_csv(os.path.join(outdir, "pmidspergene_table.tsv.gz"), sep="\t", index=False)
    
    # Plot
    __get_pmids4genes_plot(df, mean, std, outdir)

 
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

def __get_genes4pmids_plot(df, mean, std, outdir):
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
    fig.savefig(os.path.join(outdir, "geneVSpmid_distribution.png"))
    plt.clf()
    
def __get_pmids4genes_plot(df, mean, std, outdir):
    fig, ax = plt.subplots()
    labels = ["PMIDs", "Ortholog PMIDs", "All PMIDs"]
    colors = ["#1965B0", "#F7F056", "#DC050C"]
    for l, c in zip(labels, colors):
        x = []; y = []
        for i, v in df[l].value_counts().items():
            x.append(i); y.append(v)
        ax.plot(x, y, "o", label=l, color=c)
    
    kwargs = {"ls" : "--", "lw" : 1}
    ax.axvline(x=mean, color="black", label="Mean (%s)" % round(mean,3), **kwargs)
    ax.legend(frameon=False)
    ax.set(xscale="log", yscale="log")
    ax.set_xlabel("PMIDs")
    ax.set_ylabel("Genes")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(os.path.join(outdir, "pmidVSgene_distribution.png"))
    plt.clf()

if __name__ == "__main__":
    cli()
