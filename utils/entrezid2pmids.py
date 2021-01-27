#!/usr/bin/env python

import click
import os
import pandas as pd
from urllib import request

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.option(
    "-e", "--entrezid", "entrezids",
    help="Entrez Gene ID (e.g. 10664).",
    type=int,
    multiple=True,
    required=True,
)
@click.option(
    "-o", "--orthologs", "add_orthologs",
    help="Also add PubMed ID(s) of ortholog(s).",
    is_flag=True
)

def cli(**params):

    print_str = "Entrez Gene ID\tPubMed ID(s)"
    if params["add_orthologs"]:
        print_str += "\tPubMed ID(s) of ortholog(s)"
    print(print_str)
    for entrezid, pmids, orthologs_pmids in __get_entrezids_pmids(
        params["entrezids"], params["add_orthologs"]
    ):
        print_str = "%s\t%s" % (entrezid, ",".join(map(str, pmids)))
        if params["add_orthologs"]:
            print_str += "\t%s" % ",".join(map(str, orthologs_pmids))
        print(print_str)

def __get_entrezids_pmids(entrezids, add_orthologs=False):

    # Initialize
    entrezids_pmids = []
    gene2pubmed_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "gene2pubmed.gz")
    homologene_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "homologene.data")

    # Get PubMed IDs
    if not os.path.exists(gene2pubmed_file):
        request.urlretrieve(
            "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz",
            gene2pubmed_file)
    gene2pubmed = pd.read_csv(gene2pubmed_file, sep="\t",
        names=["entrezid", "pmid"], usecols=[1, 2], skiprows=1)
    if not os.path.exists(homologene_file):
        request.urlretrieve(
            "ftp://ftp.ncbi.nlm.nih.gov/HomoloGene/current/homologene.data",
            homologene_file)
    if add_orthologs:
        homologene = pd.read_csv(homologene_file, sep="\t",
            names=["cluster", "entrezid"], usecols=[0, 2])
    for eid in entrezids:
        pmids = gene2pubmed[gene2pubmed["entrezid"] == eid]["pmid"]
        entrezids_pmids.append([eid, pmids.tolist()])
        if add_orthologs:
            try:
                c = int(homologene[homologene["entrezid"] == eid]["cluster"])
                eids = homologene[homologene["cluster"] == c]["entrezid"].tolist()
                eids.remove(eid)
                pmids = gene2pubmed[gene2pubmed["entrezid"].isin(eids)]["pmid"]
                entrezids_pmids[-1].append(pmids.tolist())
            except:
                entrezids_pmids[-1].append([])
        else:
            entrezids_pmids[-1].append([])

    return(entrezids_pmids)

if __name__ == "__main__":
    cli()