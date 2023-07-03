#!/usr/bin/env python

import click
from click_option_group import (optgroup, RequiredMutuallyExclusiveOptionGroup)
import gzip
import itertools
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
import pandas as pd
import numpy as np
import ast
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@optgroup.group("Input", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "-e", "--entrezid", "entrezids",
    help="Entrez IDs (e.g. 10664 or P49711).",
    multiple=True
)
@optgroup.option(
    "--entrezid-file",
    help="List of entrez-ids.",
    type=click.File("rt")
)
@click.option(
    "--class-file",
    help="TF Class file.",
    type=click.Path(),
    required=True
)
@click.option(
    "--classname-file",
    help="File with subclass names of TF Class.",
    type=click.Path(),
    required=True
)
@click.option(
    "--geneinfo-file",
    help="Gene info file.",
    type=click.Path(),
    required=True
)
@click.option(
    "--tfidf-dir",
    help="TF-IDFs directory.",
    type=click.Path(),
    default="./tf-idfs",
    show_default=True
)
@click.option(
    "--output-dir",
    help="Output directory.",
    type=click.Path(),
    default="./",
    show_default=True
)
@click.option(
    "--zscore",
    help="Filter out stems with zscore higher than threshold.",
    type=float,
    default=None
)

def cli(**params):
    # Create output dir
    if not os.path.isdir(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get list of entrezids
    if params["entrezid_file"] is not None:
        params["entrezids"] = list(map(int, pd.read_csv(params["entrezid_file"], header=None)[0]))
    
    # Get df from files
    params["geneinfo_file"] = pd.read_table(params["geneinfo_file"], sep="\t", 
                                            names=["entrezid", "symbol"],
                                            usecols=[1, 2], skiprows=1)
    
    params["class_file"] = pd.read_table(params["class_file"], sep="\t", names=['symbol', 'class'])
    params["classname_file"] = pd.read_table(params["classname_file"], sep="\t", names=['subclass', 'full_name'])
       
    tfclass = __get_tfclass_df(params["class_file"], params["classname_file"])
    
    # Compute Statistics (Genes per Stem)
    df, entrezids2stems = __get_genes4stems_stats(params["entrezids"], params["tfidf_dir"], params["output_dir"])
    
    # Comparison
    __get_pairwisecomparison(df, params["entrezids"], tfclass, params["geneinfo_file"], entrezids2stems, params["zscore"], params["tfidf_dir"], params["output_dir"])
    
def __get_tfclass_df(tfclass, classname):
    tfclass = tfclass.astype({'class':'str'})
    tfclass['subclass'] = [ a[0:a.rfind('.')] for a in [ s[0:s.rfind('.')] for s in list(tfclass['class'])] ]
    
    classname = classname.astype({'subclass':'str'})
    
    tfclass = tfclass.merge(classname)
    
    return(tfclass)
    
def __get_genes4stems_stats(entrezids, tfidfs_dir, outdir):
    
    geneperstem = os.path.join(outdir, "genesperstem_table.tsv.gz")
    json_file = os.path.join(outdir, "entrezid2stems.json.gz")
    counts = dict()
    entrezids_stems = dict()
    
    for id in entrezids:
        tsv_file = os.path.join(tfidfs_dir, "%s.tsv.gz" % id)   
        df = pd.read_csv(tsv_file, sep="\t", header=0,
                         converters={"Stem": ast.literal_eval})
        
        if not df.empty:
            stems_df = pd.DataFrame(list(itertools.chain(*df['Stem'])), 
                                    columns=['stem']).groupby('stem', as_index=False).aggregate(list)
            
            entrezids_stems[str(id)] = list(stems_df['stem'])
        
            for s in stems_df['stem']:
                counts[s] = counts.get(s, 0) + 1
            
    df_count = pd.DataFrame.from_dict(counts, orient="index",columns=['Genes']).reset_index().rename(columns={"index":"Stem"})
    
    mean = np.mean(df_count["Genes"].tolist())
    std = np.std(df_count["Genes"].tolist())
    df_count['Zscore'] = [__Z_score(x, mean, std) for x in df_count["Genes"].tolist()]
    
    # Export Files
    df_count.to_csv(geneperstem, sep="\t", index=False, compression="gzip")
    
    with gzip.open(json_file, "wt") as handle:
        handle.write(json.dumps(entrezids_stems, indent=4))
    
            
    # Plot
    __get_genes4stems_plot(df_count, mean, std, outdir)
    
    return(df_count, entrezids_stems)

def __get_genes4stems_plot(df, mean, std, outdir):
    fig, ax = plt.subplots()
    
    x = []; y = []
    for i, v in df["Genes"].value_counts().items():
        x.append(i); y.append(v)
    ax.plot(x, y, "o")
    
    z2 = 2*std+mean
    kwargs = {"ls" : "--", "lw" : 1}
    ax.axvline(x=z2, color="black", label="Z-score = 2 (%s)" % round(z2, 3),
        **kwargs)
    ax.legend(frameon=False)
    ax.set(xscale="log", yscale="log")
    ax.set_title("Number of Gene(s) per Stem")
    ax.set_xlabel("Genes")
    ax.set_ylabel("Stems")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(os.path.join(outdir, "genesperstem_distribution.png"))
    plt.clf()

def __Z_score(x, mean, std):
    return((x-mean)/std)

def __get_pairwisecomparison(df, entrezids, df_tfclass, gene_info, entrezids2stems, zscore, tfidfs_dir, outdir):    
    
    order_entrezids, symbols, colors, legend, patch = __get_ordered_entrezids(entrezids, df_tfclass, gene_info)

    # Filter Stems to use for comparison
    if zscore is not None:
        df = df[(df["Genes"] > 1) & (df["Zscore"] < zscore)]
    else:
        df = df[df["Genes"] > 1]
    
    stems2compare = df['Stem']
    
    # Vectors
    abspres_vec = __get_vectors_abspres(stems2compare, order_entrezids, entrezids2stems)
 
    # Pairwise Comparison
    __get_abspres_cosine_class(order_entrezids, colors, abspres_vec, legend, patch, outdir)
#    __get_abspres_cosine_cluster(order_entrezids, colors, abspres_vec, outdir)
    
    __get_manual_matrix_symbol(order_entrezids, symbols, colors, tfidfs_dir, legend, patch, outdir)


def __get_ordered_entrezids(entrezids, tfclass, gene_info):    
    aliases_df = gene_info[['entrezid', 'symbol']][gene_info['entrezid'].isin(entrezids)]
        
    aliases2cluster = tfclass.merge(aliases_df, how='left')
    aliases2cluster = aliases2cluster[aliases2cluster['entrezid'].isin(entrezids)].sort_values(['class', 'entrezid'])
  
    uniq_sub = list(aliases2cluster['full_name'].unique())
    new_cmap = list(dict.fromkeys(sns.color_palette("Set1") + sns.color_palette("tab20") +  sns.color_palette("Accent") + sns.color_palette("Dark2")))
    aliases2cluster['colors'] = [ new_cmap[uniq_sub.index(subclass)] for subclass in aliases2cluster['full_name'] ]

    entrezids = list(map(int, aliases2cluster['entrezid']))
    symbols = list(aliases2cluster['symbol'])
    colors = list(aliases2cluster['colors'])
    
    legend = dict(zip(uniq_sub, new_cmap[0:len(uniq_sub)]))
    patch = [Patch(facecolor=legend[name]) for name in legend]
    
    return(entrezids, symbols, colors, legend, patch)
    

def __get_vectors_abspres(stems, entrezids, entrezid_stems):
    vectors = []

    for id in entrezids:
        abs_pres = []
        set_stems = set(entrezid_stems[str(id)])

        for s in stems:
            if s in set_stems:
                abs_pres.append(1)
            else:
                abs_pres.append(0)
        
        vectors.append(np.array(abs_pres))
         
    return(vectors)


def __get_manual_matrix_symbol(entrezids, symbols, colors, tfidfs_dir, legend, patch, outdir):
    matrix_man = []
    for i in range(len(entrezids)):
        tfidf_vec = []
        
        tfidf_file = os.path.join(tfidfs_dir, "%s.tsv.gz" % entrezids[i]) 
        tfidf_df = pd.read_table(tfidf_file)

        for s in symbols:
            match = tfidf_df[tfidf_df['Word'] == s.lower()]
#            match = tfidf_df[tfidf_df['Word'].str.contains(s.lower())]
            if not match.empty:
                tfidf_vec.append(max(match['Combo TF-IDF']))
            else:
                tfidf_vec.append(0)
        
        matrix_man.append(np.array(tfidf_vec))
        
#    pd.DataFrame(entrezids).to_csv(os.path.join(outdir, "filtentrezids.txt"), sep="\t", header=False, index=False)
     
    matrix_exact_o = pd.DataFrame(matrix_man)
    matrix_exact_r = matrix_exact_o
    matrix_exact_r = pd.DataFrame(np.rot90(np.fliplr(matrix_exact_r)))
    
    matrix_man = (matrix_exact_o * matrix_exact_r) ** (1/2)
    
#    matrix_man.to_csv(os.path.join(outdir, "matrix_man_exact.tsv"), sep="\t", header=False, index=False)    

    plt.subplots(figsize=(15,15))
    heat_file = os.path.join(outdir, "manual_tfidf_pairwise_class.png")
    heatmap_man = sns.clustermap(matrix_man, cmap='viridis', cbar_pos=(0.2, -0.05, 0.6, .025),
                             cbar_kws = {'orientation':'horizontal', 'label': 'Cosine similarity'},
                             xticklabels=[], yticklabels=[],
                             row_cluster=False, col_cluster=False,
                             dendrogram_ratio=[0.01,0.01],
                             row_colors=colors, col_colors=colors)
    plt.legend(patch, legend, bbox_to_anchor=(1.1, 0.9), title='Class',
               bbox_transform=plt.gcf().transFigure, loc='upper left')
    heatmap_man.savefig(heat_file, dpi=400)
    plt.clf()


def __get_abspres_cosine_class(entrezids, colors, vectors, legend, patch, outdir):

    heatmap_file = os.path.join(outdir, "tf_pairwise_class.png")
    matrix = cosine_similarity(vectors)

#    df_matrix = pd.DataFrame(matrix)
#    df_matrix.to_csv(os.path.join(outdir, "matrix_cor.tsv"), sep="\t", header=False, index=False)    

    heatmap = sns.clustermap(matrix, cmap='viridis', cbar_pos=(0.2, -0.05, 0.6, .025),
                             cbar_kws = {'orientation':'horizontal', 'label': 'TF-IDF'},
                             xticklabels=[], yticklabels=[],
                             vmin=0, vmax=1,
                             row_cluster=False, col_cluster=False,
                             dendrogram_ratio=[0.01,0.01],
                             row_colors=colors, col_colors=colors)
    plt.legend(patch, legend, bbox_to_anchor=(1.1, 0.9), title='Class',
               bbox_transform=plt.gcf().transFigure, loc='upper left')
    heatmap.savefig(heatmap_file, dpi=400)
    plt.clf()


def __get_abspres_cosine_cluster(entrezids, colors, vectors, outdir):
    clustermap_file = os.path.join(outdir, "tf_pairwise_clusters.png")
    matrix = cosine_similarity(vectors)

    clustermap = sns.clustermap(matrix, cmap='viridis', cbar_pos=(0.25, -0.1, 0.50, .025),
                                cbar_kws = {'orientation':'horizontal', 'label': 'Cosine similarity'},
                                xticklabels=[], yticklabels=[],
                                vmin=0, vmax=1,
                                row_cluster=False,
                                row_colors=colors, col_colors=colors)
    clustermap.savefig(clustermap_file, dpi=400)
    plt.clf()


if __name__ == "__main__":
    cli()
