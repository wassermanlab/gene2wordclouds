#!/usr/bin/env python

import argparse
from copy import deepcopy
from numpy import log10 as log
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

    parser.add_argument("-i", default="./", help="input directory (i.e. output directory from get_pmids.py)", metavar="DIR")
    parser.add_argument("-o", default="./", help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make files
    get_idfs(args.i, args.o)

def get_idfs(input_dir, output_dir="./"):

    # Initialize
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output directory
    if not os.path.exists(os.path.abspath(output_dir)):
        os.makedirs(os.path.abspath(output_dir))

    # Move to taxon directory
    os.chdir(output_dir)

    # For each taxon...
    for taxon in taxons:

        # Initialize
        rds_dir = os.path.join(os.path.abspath(input_dir), taxon)

    # Skip if pickle file already exists
    pickle_file = "idfs.%s.pickle" % taxon
    if not os.path.exists(pickle_file):

        # Initialize
        df = None
        idfs = {}
        pmids = 0

        # For each RDS file...
        for rds_file in os.listdir(rds_dir):

            # Read RDS
            result = pyreadr.read_r(os.path.join(rds_dir, rds_file))

            # Extract data frame
            next_df = result[None]

            # Append data frame
            if df is None:
                df = next_df
            else:
                # Append at the end
                df = df.append(next_df, ignore_index=True)

            # +1
            pmids += 1

        # Group data frame
        df["Freq"] = 1
        df = df.groupby("Var1").sum().reset_index()
        word2pmids = dict(zip(df.Var1, df.Freq))

        # For each word...
        for wrd in word2pmids:

            # Calculate idf
            idfs.setdefault(wrd, log(float(pmids) / word2pmids[wrd]))

        # Write pickle file
        with open(pickle_file, "wb") as f:
            pickle.dump(idfs, f)

    # Return to original directory
    os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()