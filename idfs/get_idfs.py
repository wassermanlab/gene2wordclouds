#!/usr/bin/env python

import argparse
from collections import Counter
from numpy import log10 as log
import os
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
    get_idfs(os.path.abspath(args.i), os.path.abspath(args.o))

def get_idfs(input_dir, output_dir="./"):

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
    pickle_file = os.path.join(input_dir, ".uniacc2entrezid.pickle")
    with open(pickle_file, "rb") as f:
        uniacc2entrezid = pickle.load(f)

    # Load pickle file
    pickle_file = os.path.join(input_dir, ".entrezid2pmid.pickle")
    with open(pickle_file, "rb") as f:
        entrezid2pmid = pickle.load(f)

    # For each taxon...
    for taxon in taxons:

        # Initialize
        rds_dir = os.path.join(input_dir, taxon)

        # Skip if pickle file already exists
        pickle_file = "idfs.%s.pickle" % taxon
        if not os.path.exists(pickle_file):

            # Initialize
            idfs = {}
            words = []
            counter = 0

            # Load pickle file
            with open(os.path.join(input_dir, ".uniaccs.%s.pickle" % taxon), "rb") as f:
                uniaccs = pickle.load(f)

            # For each uniacc...
            for uniacc in sorted(uniaccs):

                # Initialize
                df = None
                word_set = set()

                # Skip if uniacc not mapped to an entrezid
                if uniacc not in uniacc2entrezid:
                    continue

                # Skip if entrezid not mapped to a pmid
                if uniacc2entrezid[uniacc] not in entrezid2pmid:
                    continue

                # For each pmid
                for pmid in entrezid2pmid[uniacc2entrezid[uniacc]]:

                    # Skip if no RDS file
                    rds_file = os.path.join(input_dir, taxon, "%s.rds" % pmid)
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

                # Add words
                for word in df["Var1"].to_list():
                    words.append(word)

                # +1
                counter += 1

            # For each word...
            for word, occurrences in Counter(words).items():

                # Calculate idf
                idfs.setdefault(word, log(float(counter) / occurrences))

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