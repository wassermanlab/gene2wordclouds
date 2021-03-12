#!/usr/bin/env python

import click
import pandas as pd

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "protein_info",
    type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "protein_links",
    type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-f", "--filter-genes",
    help="List of genes to filter by.",
    type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-o", "--output-file",
    help="Output file.  [default: STDOUT]",
    type=click.Path(writable=True, readable=False, resolve_path=True,
        allow_dash=True),
    default="-"
)

def cli(**params):

    # Get info
    df = pd.read_csv(params["protein_info"], sep="\t", usecols=[0, 1])
    info = dict(
        zip(df.protein_external_id.to_list(), df.preferred_name.to_list())
    )

    # Get genes
    if params["filter_genes"] is not None:
        df = pd.read_csv(params["filter_genes"], header=None)
        genes = set(df[0].to_list())
    else:
        genes = set(info.values())


    # For each row...
    interactions = []
    chunks = pd.read_csv(params["protein_links"], sep="\t", chunksize=1000)
    for chunk in chunks:
        for _, row in chunk.iterrows():
            protein1, protein2, score = row.values[0].split()
            if info[protein1] not in genes:
                    continue
            if info[protein2] not in genes:
                    continue
            if int(score) >= 900:
                confidence = "highest"
            elif int(score) >= 700:
                confidence = "high"
            elif int(score) >= 400:
                confidence = "medium"
            elif int(score) >= 150:
                confidence = "low"
            else:
                continue
            interaction = [info[protein1], info[protein2], score, confidence]
            interactions.append(interaction)

    # Write
    column_names = ["Gene1", "Gene2", "Score", "Confidence"]
    df = pd.DataFrame(interactions, columns=column_names)
    with click.open_file(params["output_file"], "wt") as f:
        df.to_csv(f, sep="\t")

if __name__ == "__main__":
    cli()