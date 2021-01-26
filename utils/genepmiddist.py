import sys
import os
import numpy as np
import pandas as pd
import json
import gzip
import matplotlib.pyplot as plt

def genepmiddist():
    with gzip.open(sys.argv[1], 'rt') as file:
        global entrezid2pmids
        entrezid2pmids = json.load(file)
        file.close

    cwd = os.getcwd()
    outdir = os.path.join(cwd, 'genevspmid')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)

    parseelements()
    countgeneperpmid()
    getstats()
    getvisual()

def parseelements():
    entrezids = []
    pmids = []
    homolopmids = []
    for line in entrezid2pmids:
        entrezids.append(line[0])
        pmids.append(line[1])
        homolopmids.append(line[2])

    global flat_pmids
    flat_pmids = []
    for e in pmids:
        flat_pmids += e

    global flat_homolopmids
    flat_homolopmids = []
    for e in homolopmids:
        flat_homolopmids += e

def countgeneperpmid():
    genecount = dict()
    for pmid in flat_pmids:
        genecount[pmid] = genecount.get(pmid, 0) + 1

    df_pmid_genecount = pd.DataFrame.from_dict(genecount, orient="index").rename(columns={0 : "pmid"})
    df_pmid_count = pd.DataFrame(df_pmid_genecount["pmid"].value_counts())

    for pmid in flat_homolopmids:
        genecount[pmid] = genecount.get(pmid, 0) + 1

    hompmid_genecount = pd.DataFrame.from_dict(genecount, orient="index").rename(columns={0 : "homo+pmid"})
    df_hompmid_count = pd.DataFrame(hompmid_genecount["homo+pmid"].value_counts())

    total_genecount = df_pmid_genecount.join(hompmid_genecount)

    global df_total_sum
    df_total_sum = df_pmid_count.join(df_hompmid_count).sort_index()

def getstats():
    df_total_sum["zscore_pmid"] = (df_total_sum['pmid'] - df_total_sum['pmid'].mean())/df_total_sum['pmid'].std()
    df_total_sum["zscore_homo+pmid"] = (df_total_sum['homo+pmid'] - df_total_sum['homo+pmid'].mean())/df_total_sum['homo+pmid'].std()
    df_total_sum.to_csv('numgenes_numpmid_zscore.txt', sep="\t", index=True)

    mean = [df_total_sum['pmid'].mean(), df_total_sum['homo+pmid'].mean()]
    median = [df_total_sum['pmid'].median(), df_total_sum['homo+pmid'].median()]
    stdev = [df_total_sum['pmid'].std(), df_total_sum['homo+pmid'].std()]

    df_stats = pd.DataFrame(list(zip(mean, median, stdev)),
                     columns = ['Mean', 'Median', 'Stdev'],
                     index = ['pmid', 'homo+pmid'])
    df_stats.to_csv('summarystats.txt', sep="\t", index=True)

def getvisual():
    plt.plot(df_total_sum.index.tolist(), df_total_sum['pmid'], label='pmid', color="blue");
    plt.plot(df_total_sum.index.tolist(), df_total_sum['homo+pmid'], label='homo+pmid', color="red");
    plt.axvline(x=1, ymin=0, color="lightblue", label="zscore_pmid > 2")
    plt.axvline(x=2, ymin=0, color="lightpink", label="zscore_homo+pmid > 2")
    plt.xscale('log')
    plt.xlabel("# of genes")
    plt.yscale('log')
    plt.ylabel("# of pmids")
    plt.legend()

    plt.savefig('plot.png')


if __name__ == "__main__":
    genepmiddist()
