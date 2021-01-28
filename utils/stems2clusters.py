#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import gzip
import pandas as pd

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@optgroup.group("Words", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--word", "words",
    help="Word (e.g. \"11-zinc-finger\").",
    multiple=True
)
@optgroup.option(
    "--words-file",
    help="File with words.",
    type=click.File("rt")
)
@optgroup.group("Stems", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--stem", "stems",
    help="Stem (e.g. \"11,finger,zinc\").",
    multiple=True
)
@optgroup.option(
    "--stems-file",
    help="File with stems per word.",
    type=click.File("rt")
)
@optgroup.group("Stem Counts", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--count", "counts",
    help="Stem counts/occurrances based on the word associated (e.g. 4).",
    type=int,
    multiple=True
)
@optgroup.option(
    "--counts-file",
    help="File with stem counts.",
    type=click.File("rt")
)
@optgroup.group("Weights", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--weight", "weights",
    help="Weight/TF-IDFs (e.g. 1.0).",
    type=float,
    multiple=True
)
@optgroup.option(
    "--weights-file",
    help="File with weights/tf-idfs.",
    type=click.File("rt")
)

def cli(**params):

    # Parse words
    if params["words_file"] is not None:
        params["words"] = []
        if params["words_file"].name.endswith(".gz"):
            handle = gzip.open(params["words_file"].name, "rt")
            params["words_file"].close()
        else:
            handle = params["words_file"]
        for line in handle:
            params["words"].append(line.strip("\n"))
        handle.close()
        
    # Parse stems
    if params["stems_file"] is not None:
        params["stems"] = []
        if params["stems_file"].name.endswith(".gz"):
            handle = gzip.open(params["stems_file"].name, "rt")
            params["stems_file"].close()
        else:
            handle = params["stems_file"]
        for line in handle:
            params["stems"].append(line.strip("\n"))
        handle.close()

    # Parse stem counts
    if params["counts_file"] is not None:
        params["counts"] = []
        if params["counts_file"].name.endswith(".gz"):
            handle = gzip.open(params["counts_file"].name, "rt")
            params["counts_file"].close()
        else:
            handle = params["counts_file"]
        for line in handle:
            params["counts"].append(int(line.strip("\n")))
        handle.close()
    
    # Parse weights
    if params["weights_file"] is not None:
        params["weights"] = []
        if params["weights_file"].name.endswith(".gz"):
            handle = gzip.open(params["weights_file"].name, "rt")
            params["weights_file"].close()
        else:
            handle = params["weights_file"]
        for line in handle:
            params["weights"].append(int(float(line.strip("\n"))))
        handle.close()

    for w, we, c in __get_word_clusters(params["words"], params["stems"], params["counts"], params["weights"]):
        print("%s\t%s\t%s" % (w, we, c))

def __get_word_clusters(words, stems, counts, weights):

    # Initialize
    df = pd.DataFrame(list(zip(words, stems, counts, weights)),
                    columns = ['word', 'stem', 'stem_count', 'weight'])
    df['stem'] = df['stem'].apply(lambda row: row.split(","))
    
    counts = dict()
    clusters = []
    
    # Rank stems based on frequency of occurances
    for index, row in df.iterrows():
        for s in row['stem']:
            counts[s] = counts.get(s, 0) + row['stem_count']

    stem_cluster = pd.DataFrame.from_dict(counts, orient="index", columns=['count']).sort_values(by='count').index.tolist()

    # Assign each word to their highest rank stem cluster

    for index, row in df.iterrows():
        rank = -1
        for s in row['stem']:
            if stem_cluster.index(s) > rank:
                rank = stem_cluster.index(s)
                cluster = s            
        clusters.append(cluster)
    
    df['cluster'] = clusters

    # Group words based on cluster
    df_filter = df.groupby('cluster').head(1)
    
    # Return the word and the weight
    return([[r['word'], r['weight'], r['cluster']] for i, r in df_filter.iterrows()])

if __name__ == "__main__":
    cli()