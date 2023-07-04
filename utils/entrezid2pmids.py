#!/usr/bin/env python

import click
import os
import pandas as pd
from urllib import request

from .entrezid2aliases import __is_file_7_days_old

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
    "-o", "--orthologs",
    help="Also add PubMed ID(s) of ortholog(s).",
    is_flag=True
)

def cli(**params):

    gene2pubmed, homologene = __load_datasets(params["orthologs"])

    print_str = "Entrez Gene ID\tPubMed ID(s)"
    if params["orthologs"]:
        print_str += "\tPubMed ID(s) of ortholog(s)"
    print(print_str)
    for entrezid in params["entrezids"]:
        pmids, orthologs_pmids = __get_entrezid_pmids(entrezid, gene2pubmed,
            homologene)
        print_str = "%s\t%s" % (entrezid, ",".join(map(str, pmids)))
        if params["orthologs"]:
            print_str += "\t%s" % ",".join(map(str, orthologs_pmids))
        print(print_str)

def __load_datasets(orthologs=False):

    # Initialize
    homologene = None
    utils_dir = os.path.dirname(os.path.realpath(__file__))
    gene2pubmed_file = os.path.join(utils_dir, "data", "gene2pubmed.gz")
    homologene_file = os.path.join(utils_dir, "data", "homologene.data")

    # Get PubMed IDs
    if not os.path.exists(gene2pubmed_file) \
       or __is_file_7_days_old(gene2pubmed_file):
        request.urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz",
            gene2pubmed_file)
    gene2pubmed = pd.read_csv(gene2pubmed_file, sep="\t",
        names=["entrezid", "pmid"], usecols=[1, 2], skiprows=1)
    if not os.path.exists(homologene_file) \
       or __is_file_7_days_old(homologene_file):
        request.urlretrieve(
            "https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",
            homologene_file)
    if orthologs:
        homologene = pd.read_csv(homologene_file, sep="\t",
            names=["cluster", "entrezid"], usecols=[0, 2])        

    return(gene2pubmed, homologene)

def __get_entrezid_pmids(entrezid, gene2pubmed, homologene=None):

    # Get Pubmed IDs
    s = gene2pubmed[gene2pubmed["entrezid"] == entrezid]["pmid"]
    pmids = s.tolist()
    if homologene is not None:
        try:
            c = int(homologene[homologene["entrezid"] == entrezid]["cluster"])
            l = homologene[homologene["cluster"] == c]["entrezid"].tolist()
            l.remove(entrezid)
            s = \
                    gene2pubmed[gene2pubmed["entrezid"].isin(l)]["pmid"]
            orthologs_pmids = s.tolist()
        except:
            orthologs_pmids = []
    else:
        orthologs_pmids = []

    return(pmids, orthologs_pmids)

if __name__ == "__main__":
    cli()