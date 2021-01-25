#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import hashlib
import os
import sys

# Import from lib
sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "lib")
)
from entrezid2pmid import __get_entrezids_pmids
from uniacc2entrezid import __get_uniaccs_entrezids

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@optgroup.group("Input", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "-i", "--identifier", "identifiers",
    help="Identifier (e.g. 10664 or P49711).",
    multiple=True
)
@optgroup.option(
    "--input-file",
    help="List of identifiers.",
    type=click.File("rt")
)
@click.option(
    "--id-type",
    help="Identifier type.",
    type=click.Choice(["entrezid", "uniacc"], case_sensitive=False),
    default="entrezid",
    show_default=True
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    default="./",
    show_default=True
)
@click.option(
    "-p", "--prefix",
    help="Prefix.  [default: md5 digest]",
)

def wordcloud(identifiers, input_file, id_type, output_dir, prefix):

    # Parse identifiers
    if input_file is not None:
        identifiers = []
        for line in input_file:
            if id_type == "uniacc":
                identifiers.append(line.strip("\n"))
            else:
                identifiers.append(int(line.strip("\n")))
        input_file.close()

    # Get MD5
    if prefix is None:
        h = hashlib.md5()
        h.update(",".join(identifiers).encode("utf-8"))
        prefix = h.hexdigest()

    # Create output dir
    output_dir = os.path.join(output_dir, prefix)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    exit(0)

    # UniAcc to EntrezID
    if id_type == "uniacc":
        mappings = __get_uniaccs_entrezids(identifiers)
        uniaccs, identifiers = list(map(list, zip(*mappings)))

    # EntrezID to PMID
    mappings = __get_entrezids_pmids(identifiers, add_orthologs=True)
    entrezids, pmids, pmids_orthologs = list(map(list, zip(*mappings)))
    print(entrezids, pmids, pmids_orthologs)

if __name__ == "__main__":
    wordcloud()