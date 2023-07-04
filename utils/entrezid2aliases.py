#!/usr/bin/env python

import click
import datetime
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
    "-o", "--orthologs",
    help="Also add PubMed ID(s) of ortholog(s).",
    is_flag=True
)

def cli(**params):

    gene_info, homologene = __load_datasets(params["orthologs"])
    print_list = ["Entrez Gene ID", "Gene Symbol", "Full Name",
        "Synonym(s)"]
    if params["orthologs"]:
        print_list += ["Gene Symbol(s) of ortholog(s)", 
            "Full Name(s) of ortholog(s)", "Synonym(s) of ortholog(s)"]
    print("\t".join(print_list))
    for entrezid in params["entrezids"]:
        (gene_symbol, full_name, synonyms, gene_symbol_orthologs,
            full_name_orthologs, synonyms_orthologs) = \
        __get_entrezid_aliases(entrezid, gene_info, homologene)
        print_list = [gene_symbol, full_name, ",".join(map(str, synonyms))]
        if params["orthologs"]:
            print_list += [",".join(map(str, gene_symbol_orthologs))]
            print_list += [",".join(map(str, full_name_orthologs))]
            print_list += [",".join(map(str, synonyms_orthologs))]
        print("%s\t%s" % (entrezid, "\t".join(print_list)))

def __load_datasets(orthologs=False):

    # Initialize
    homologene = None
    utils_dir = os.path.dirname(os.path.realpath(__file__))
    gene_info_file = os.path.join(utils_dir, "data", "gene_info.gz")
    homologene_file = os.path.join(utils_dir, "data", "homologene.data")

    # Get PubMed IDs
    if not os.path.exists(gene_info_file) \
       or __is_file_7_days_old(gene_info_file):
        request.urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
            gene_info_file)
    names = ["entrezid", "symbol", "synonyms", "official_symbol",
        "official_full_name"]
    gene_info = pd.read_csv(gene_info_file, sep="\t", names=names,
        usecols=[1, 2, 4, 10, 11], skiprows=1)
    if not os.path.exists(homologene_file) \
       or __is_file_7_days_old(homologene_file):
        request.urlretrieve(
            "https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",
            homologene_file)
    if orthologs:
        homologene = pd.read_csv(homologene_file, sep="\t",
            names=["cluster", "entrezid"], usecols=[0, 2])        

    return(gene_info, homologene)

def __get_entrezid_aliases(entrezid, gene_info, homologene=None):

    # Get aliases
    df = gene_info[gene_info["entrezid"] == entrezid]
    gene_symbol = df["symbol"].to_string(index=False)
    full_name  = df["official_full_name"].to_string(index=False)
    synonyms = set(df["synonyms"].to_string(index=False).split("|"))
    synonyms.add(df["official_symbol"].to_string(index=False))
    if gene_symbol in synonyms:
        synonyms.remove(gene_symbol)
    synonyms = list(synonyms)
    if homologene is not None:
        try:
            c = int(homologene[homologene["entrezid"] == entrezid]["cluster"])
            l = homologene[homologene["cluster"] == c]["entrezid"].tolist()
            l.remove(entrezid)
            df = gene_info[gene_info["entrezid"].isin(l)]
            gene_symbol_orthologs = list(set(df["symbol"].tolist()))
            full_name_orthologs = list(set(df["official_full_name"].tolist()))
            synonyms_orthologs = set()
            for ortholog_synonyms in df["synonyms"].tolist():
                for synonym in ortholog_synonyms.split("|"):
                    synonyms_orthologs.add(synonym)
            for ortholog_official_symbol in df["official_symbol"].tolist():
                synonyms_orthologs.add(ortholog_official_symbol)
            for ortholog_gene_symbol in gene_symbol_orthologs:
                if ortholog_gene_symbol in synonyms_orthologs:
                    synonyms_orthologs.remove(ortholog_gene_symbol)
            synonyms_orthologs = list(synonyms_orthologs)
        except:
            gene_symbol_orthologs = []
            full_name_orthologs = []
            synonyms_orthologs = []
    else:
        gene_symbol_orthologs = []
        full_name_orthologs = []
        synonyms_orthologs = []

    return(gene_symbol, full_name, synonyms, gene_symbol_orthologs,
        full_name_orthologs, synonyms_orthologs)

def __is_file_7_days_old(file_name):

    # Get today's date
    today = datetime.datetime.today()

    # Get last modified date for file
    modified_date = datetime.datetime.fromtimestamp(
        os.path.getmtime(file_name)
    )

    # Compute age of the file
    age = today - modified_date

    return age.days >= 7


if __name__ == "__main__":
    cli()