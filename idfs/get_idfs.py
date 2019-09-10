#!/usr/bin/env python

import argparse
import coreapi
import gzip
import json
from multiprocessing import Pool
import pickle
import os
import re
import ssl
import subprocess
from tqdm import tqdm
from urllib import parse, request

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
        taxon_dir = os.path.join(os.path.abspath(output_dir), taxon)

        print(taxon_dir)

    # Return to original directory
    os.chdir(cwd)

def _get_idf(taxon, taxon_dir):

    # Skip if already downloaded
    rds_file = "%s.rds" % taxon
    if not os.path.exists(rds_file):

        # Get pmid
        cmd = "Rscript ../get_idf.R %s" % taxon_dir
        process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()