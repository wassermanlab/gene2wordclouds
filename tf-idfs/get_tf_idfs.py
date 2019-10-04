#!/usr/bin/env python

import argparse
from copy import deepcopy
from enchant.tokenize import get_tokenizer
from fuzzywuzzy import fuzz
from numpy import log10 as log, nan
import os
import pandas
import pickle
import math
import numpy as  np
import pyreadr
import re
import ssl
import subprocess
from urllib import parse, request

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", default="./../idfs", help="idfs directory (i.e. output directory from get_idfs.py)", metavar="DIR")
    parser.add_argument("-p", default="./../pmids", help="pmids directory (i.e. output directory from get_pmids.py)", metavar="DIR")
    parser.add_argument("-o", default="./", help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make files
    get_tf_idfs(os.path.abspath(args.i), os.path.abspath(args.p), os.path.abspath(args.o))

def get_tf_idfs(idfs_dir, pmids_dir, output_dir="./"):

    # Initialize
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]
    regexp = re.compile("^\w.+\w$")

    # Globals
    global cwd
    cwd = os.getcwd()

    tfs = {
        #"P04637": "TP53",
        #"P19838": "NFKB1",
        "P49711": "CTCF",
    }

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load pickle file
    pickle_file = os.path.join(pmids_dir, ".uniacc2entrezid.pickle")
    with open(pickle_file, "rb") as f:
        uniacc2entrezid = pickle.load(f)

    # Load pickle file
    pickle_file = os.path.join(pmids_dir, ".entrezid2pmid.pickle")
    with open(pickle_file, "rb") as f:
        entrezid2pmid = pickle.load(f)

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

        # Get uniacc to gene name mappings
        uniacc2genename = _get_uniacc_to_genename_mappings(uniaccs)

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

                # Filter weird words + TF containing words
                tknzr = get_tokenizer("en")
                gene_name = uniacc2genename[uniacc]
                thresh = 85
                for idx, row in df.iterrows():
                    if not regexp.match(row["Var1"]):
                        tokens = [w for w in tknzr(row["Var1"])]
                        if len(tokens) == 0:
                            df["tfidf"][idx] = -1
                        elif tokens[0][1] > 0:
                            df["tfidf"][idx] = -1
                        elif (len(tokens[-1][0]) + tokens[-1][1]) < len(row["Var1"]):
                            df["tfidf"][idx] = -1
                    # Fuzzy string matching
                    if fuzz.partial_ratio(gene_name.lower(), row["Var1"]) > 50:
                        #print(gene_name.lower(), row["Var1"], fuzz.partial_ratio(gene_name.lower(), row["Var1"]))
                        df["tfidf"][idx] = -1
                    #filter redundant words
                    # if word2 is > threshold to word1 with higher tfidf score than word2[tfidf] = -1
                df = df.sort_values(by="tfidf", ascending=False)
                df = df.iloc[:75,:].reset_index()
                for idx, row in df.iterrows():
                    for idx_current, row_current in df.iloc[:idx,:].iterrows():
                        if fuzz.partial_ratio(row_current["Var1"], row["Var1"]) > thresh:
                            df["tfidf"][idx] = -1
                df = df.sort_values(by="tfidf", ascending=False)

                # for token in tokens:
                    #     if fuzz.partial_ratio(gene_name.lower(), token) > thresh or fuzz.partial_ratio(token, gene_name.lower()) > thresh:
                            
                    #         df["tfidf"][idx] = -1
                    #         break

                # Write RDS file
                pyreadr.write_rds(out_file, df)

            # Make word cloud
            # _make_word_cloud(uniacc)
            _make_word_cloud(tfs[uniacc])

        # Return to original directory
        os.chdir(cwd)

def _get_uniacc_to_genename_mappings(uniaccs):

    # Initialize
    uniacc2genename = {}
    gcontext = ssl.SSLContext()
    url = "https://www.uniprot.org/uploadlists/"

    # # Skip if pickle file already exists
    # pickle_file = ".uniacc2genename.pickle"
    # if not os.path.exists(pickle_file):

    # Set query
    params = {
        "from": "ACC+ID",
        "to": "GENENAME",
        "format": "tab",
        "query": " ".join(list(uniaccs))
    }

    # Encode parameters
    params = parse.urlencode(query=params).encode("utf-8")

    # Make request
    req = request.Request(url, params)

    # Get response
    with request.urlopen(req, context=gcontext) as f:
        response = f.read().decode("utf-8")

    # For each line...
    for line in response.split("\n"):

        # i.e. a mapping
        m = re.search("(\w+)\t(\w+)", line)
        if m:

            # Map uniacc to genename
            uniacc, genename = line.split("\t")
            uniacc2genename.setdefault(uniacc, genename)

    #     # Write pickle file
    #     with open(pickle_file, "wb") as f:
    #         pickle.dump(uniacc2genename, f)

    # # Load pickle file
    # with open(pickle_file, "rb") as f:
    #     uniacc2genename = pickle.load(f)

    return(uniacc2genename)

def _make_word_cloud(uniacc):

    # Skip if already created
    png_file = "%s.png" % uniacc
    if not os.path.exists(png_file):

        # Get pmid
        cmd = "Rscript ../make_word_cloud.R %s" % uniacc
        process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
