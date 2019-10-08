#!/usr/bin/env python

import argparse
from copy import deepcopy
from fuzzywuzzy import fuzz
import gzip
from itertools import product
from numpy import log10 as log, nan
import os
import pandas
import pickle
import pyreadr
import re
import ssl
import subprocess
from unidecode import unidecode
from urllib import parse, request

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", default="./", help="idfs directory (i.e. output directory from get_idfs.py)", metavar="DIR")
    parser.add_argument("-p", default="./", help="pmids directory (i.e. output directory from get_pmids.py)", metavar="DIR")
    parser.add_argument("-o", default="./", help="output directory (default = ./)", metavar="DIR")
    parser.add_argument("-w", default=50, type=int, help="max. number of words (default = 50)", metavar="INT")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make files
    get_tf_idfs(os.path.abspath(args.i), os.path.abspath(args.p), os.path.abspath(args.o), args.w)

def get_tf_idfs(idfs_dir, pmids_dir, output_dir="./", max_words=50):

    # Initialize
    regexp = re.compile("^\w.+\w$")
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]
    thresh = 89.0
    aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    # Get all aminoacid pairs
    aminoacid_pairs = set(["".join(c) for c in product([a.lower() for a in aa], repeat=2)])

    # Globals
    global cwd
    cwd = os.getcwd()

    tfs = {

        "P01100": "FOS",
        "P05412": "JUN",
        "O75364": "PITX3",
        "P04637": "TP53",
        "P19838": "NFKB1",
        "P49711": "CTCF",
        "P41235": "HNF4A",
    }

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Move to output directory
    os.chdir(output_dir)

    # Load pickle file
    pickle_file = os.path.join(pmids_dir, ".uniacc2entrezid.pickle")
    with open(pickle_file, "rb") as f:
        uniacc2entrezid = pickle.load(f)

    # Load pickle file
    pickle_file = os.path.join(pmids_dir, ".entrezid2pmid.pickle")
    with open(pickle_file, "rb") as f:
        entrezid2pmid = pickle.load(f)

    # Get entrezid to gene name mappings
    entrezid2genename = _get_uniacc_to_genename_mappings(uniacc2entrezid)

    # For each taxon...
    for taxon in taxons:

        if taxon != "vertebrates":
            continue

        # Initialize
        taxon_dir = os.path.join(output_dir, taxon)
        rds_dir = os.path.join(pmids_dir, taxon)

        # Create taxon directory
        if not os.path.exists(taxon_dir):
            os.makedirs(taxon_dir)

        # Move to taxon directory
        os.chdir(taxon_dir)

        # Load pickle file
        with open(os.path.join(pmids_dir, ".uniaccs.%s.pickle" % taxon), "rb") as f:
            uniaccs = pickle.load(f)

        # Load pickle file
        with open(os.path.join(idfs_dir, "idfs.%s.pickle" % taxon), "rb") as f:
            idfs = pickle.load(f)

        for uniacc in sorted(uniaccs):

            if uniacc not in tfs:
                continue

            # Skip if RDS file already exists
            # out_file = os.path.join(taxon_dir, "%s.rds" % uniacc)
            out_file = os.path.join(taxon_dir, "%s.rds" % tfs[uniacc])
            if not os.path.exists(out_file):

                # Initialize
                df = None
                tfidfs = {}
                word_count = 0

                # Skip if uniacc not mapped to an entrezid
                if uniacc not in uniacc2entrezid:
                    continue

                # Skip if entrezid not mapped to a pmid
                if uniacc2entrezid[uniacc] not in entrezid2pmid:
                    continue

                # For each pmid
                for pmid in entrezid2pmid[uniacc2entrezid[uniacc]]:
                    
                    # Skip if no RDS file
                    rds_file = os.path.join(pmids_dir, taxon, "%s.rds" % pmid)
                    if not os.path.exists(rds_file):
                        continue

                    # Read RDS
                    result = pyreadr.read_r(os.path.join(rds_dir, rds_file))

                    # Extract data frame
                    pmid_df = result[None]

                    # Append data frame
                    if df is None:
                        df = pmid_df
                    else:
                        # Append at the end
                        df = df.append(pmid_df, ignore_index=True)

                # Skip empty data frames
                if df is None:
                    continue

                # Set word freq per PMID to 1
                df["Freq"] = 1

                # Group data frame
                df = df.groupby("Var1").sum().reset_index()

                # Get TFs
                df["tf"] = df["Freq"] / df.shape[0]

                # Get IDFs
                df["idf"] = nan
                for idx, row in df.iterrows():                   
                    df["idf"][idx] = idfs[row["Var1"]]

                # Get TF-IDFs
                df["tfidf"] = df["idf"] * df["tf"]

                # Sort data frame by TF-IDFs
                df = df.sort_values(by="tfidf", ascending=False).reset_index()

                # Filter weird words
                for idx, row in df.iterrows():

                    if word_count >= max_words:
                        df["tfidf"][idx] = -1

                    if df["tfidf"][idx] != -1:

                        # Filter weird words
                        if unidecode(row["Var1"]) in aminoacid_pairs:
                            df["tfidf"][idx] = -1
                            continue

                        if not regexp.match(unidecode(row["Var1"])):
                            df["tfidf"][idx] = -1
                            continue

                        # Filter redundant words
                        for nextidx, nextrow in df.iterrows():

                            # i.e. time saver
                            if nextidx > 250 or nextidx <= idx or df["tfidf"][nextidx] == -1:
                                continue

                            ratio = fuzz.ratio(unidecode(row["Var1"]), unidecode(nextrow["Var1"]))
                            partial_ratio = fuzz.partial_ratio(unidecode(row["Var1"]), unidecode(nextrow["Var1"]))

                            if ratio >= thresh or partial_ratio >= thresh:
                                df["tfidf"][nextidx] = -1

                        # Tokenize
                        tkns = re.findall(r'\w+', unidecode(row["Var1"]))

                        # Filter TF-containing words
                        for genename in entrezid2genename[uniacc2entrezid[uniacc]]:

                            global_ratio = fuzz.ratio(genename, unidecode(row["Var1"]))

                            # Initialize
                            gene_tkns = re.findall(r'\w+', genename.lower())

                            for tkn in tkns:

                                if len(tkn) == 1:
                                    continue

                                if df["tfidf"][idx] == -1:
                                    break

                                for gene_tkn in gene_tkns:

                                    if len(gene_tkn) == 1:
                                        continue

                                    ratio = fuzz.ratio(gene_tkn, tkn)
                                    partial_ratio = fuzz.partial_ratio(gene_tkn, tkn)

                                    if (ratio + partial_ratio) / 2.0 >= thresh:
                                        df["tfidf"][idx] = -1
                                        break

                    if df["tfidf"][idx] != -1:
                        word_count += 1

                # Write RDS file
                pyreadr.write_rds(out_file, df)

            # Make word cloud
            # _make_word_cloud(uniacc)
            _make_word_cloud(tfs[uniacc], max_words)

        # Return to original directory
        os.chdir(cwd)

def _get_uniacc_to_genename_mappings(uniacc2entrezid):

    # Initialize
    entrezids = set()
    entrezid2genename = {}
    url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"
    file_name = "gene_info.gz"

    # Skip if pickle file already exists
    pickle_file = ".entrezid2genename.pickle"
    if not os.path.exists(pickle_file):

        # For each uniacc, entrezid pair...
        for uniacc, entrezid in uniacc2entrezid.items():
            entrezids.add(entrezid)

        # Download gene to pubmed mappings
        request.urlretrieve(os.path.join(url, file_name), file_name)

        with gzip.open(file_name, "r") as f:

            # For each line
            for line in f:

                # Get gene info
                info = line.decode("utf-8").split("\t")
                entrezid = info[1]
                genename = info[2]
                synonyms = info[4].split("|")

                # If valid entrezid...
                if info[1] not in entrezids:
                    continue

                # Map entrezid to gene names
                entrezid2genename.setdefault(entrezid, [])
                entrezid2genename[entrezid].append(genename)
                for synonym in synonyms:
                    entrezid2genename[entrezid].append(synonym)

        if os.path.exists(file_name):
            os.remove(file_name)

        # Write pickle file
        with open(pickle_file, "wb") as f:
            pickle.dump(entrezid2genename, f)

    # Load pickle file
    with open(pickle_file, "rb") as f:
        entrezid2genename = pickle.load(f)

    return(entrezid2genename)

def _make_word_cloud(uniacc, max_words=50):

    # Skip if already created
    png_file = "%s.png" % uniacc
    if not os.path.exists(png_file):

        # Get pmid
        cmd = "Rscript ../make_word_cloud.R %s %s" % (uniacc, max_words)
        process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()