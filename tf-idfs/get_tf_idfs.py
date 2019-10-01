#!/usr/bin/env python

import argparse
from copy import deepcopy
from numpy import log10 as log, nan
import os
import pandas
import pickle
import pyreadr

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

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make files
    get_tf_idfs(os.path.abspath(args.i), os.path.abspath(args.p), os.path.abspath(args.o))

def get_tf_idfs(idfs_dir, pmids_dir, output_dir="./"):

    # Initialize
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Move to taxon directory
    os.chdir(output_dir)

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

        # if taxon != "vertebrates":
        #     continue

        # Skip if taxon directory already exists
        taxon_dir = os.path.join(output_dir, taxon)
        if not os.path.exists(taxon_dir):

            # Initialize
            rds_dir = os.path.join(pmids_dir, taxon)

            # Create taxon directory
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

                # if uniacc != "P49711":
                #     continue

                # Skip if RDS file already exists
                rds_file = "tf-idfs.%s.rds" % taxon
                if not os.path.exists(rds_file):

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

                    # Group data frame
                    df = df.groupby("Var1").sum().reset_index()

                    # Get TFs
                    df["tf"] = df["Freq"] / df.shape[0]

                    # Initialize IDFs
                    df["idf"] = nan

                    # For each idx, row...
                    for idx, row in df.iterrows():

                        # Get IDFs
                        df["idf"][idx] = idfs[row["Var1"]]

                    # Get TF-IDFs
                    df["tf-idf"] = df["idf"] * df["tf"]

                # Write RDS file
                pyreadr.write_rds(rds_file, df)

        # Return to original directory
        os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()