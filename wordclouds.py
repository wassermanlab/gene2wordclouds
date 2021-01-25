#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import gzip
import hashlib
from itertools import chain, zip_longest
import json
import os
import pandas as pd
import sys

# Import from lib
sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "lib")
)
from abstract2words import __get_abstract_words
from entrezid2pmids import __get_entrezids_pmids
from pmid2abstract import __get_pmids_abstracts
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
    "-e", "--email",
    help="E-mail address.",
    required=True
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

def wordcloud(identifiers, input_file, email, id_type, output_dir, prefix):

    # Parse identifiers
    if input_file is not None:
        identifiers = []
        if input_file.name.endswith(".gz"):
            with gzip.open(input_file.name, "rt") as handle:
                for line in handle:
                    if id_type == "uniacc":
                        identifiers.append(line.strip("\n"))
                    else:
                        identifiers.append(int(line.strip("\n")))
        else:
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

    # Create dirs
    output_dir = os.path.join(output_dir, prefix)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    abstracts_dir = os.path.join(output_dir, "abstracts")
    if not os.path.isdir(abstracts_dir):
        os.makedirs(abstracts_dir)
    words_dir = os.path.join(output_dir, "words")
    if not os.path.isdir(words_dir):
        os.makedirs(words_dir)

    # UniAcc to EntrezID
    if id_type == "uniacc":
        json_file = os.path.join(output_dir, "uniacc2entrezid.json.gz")
        if not os.path.exists(json_file):
            mappings = __get_uniaccs_entrezids(identifiers)
            with gzip.open(json_file, "wt") as handle:
                handle.write(json.dumps(mappings, sort_keys=True, indent=4))
        with gzip.open(json_file, "rt") as handle:
            mappings = json.load(handle)
        uniaccs, identifiers = list(map(list, zip(*mappings)))

    # EntrezID to PMIDs
    json_file = os.path.join(output_dir, "entrezid2pmids.json.gz")
    if not os.path.exists(json_file):
        mappings = __get_entrezids_pmids(identifiers, add_orthologs=True)
        with gzip.open(json_file, "wt") as handle:
            handle.write(json.dumps(mappings, sort_keys=True, indent=4))
    with gzip.open(json_file, "rt") as handle:
        mappings = json.load(handle)    
    entrezids, pmids, pmids_orthologs = list(map(list, zip(*mappings)))

    # PMID to Abstract
    all_pmids = set(list(chain.from_iterable(pmids+pmids_orthologs)))
    for abstract_file in os.listdir(abstracts_dir):
        if abstract_file.endswith(".txt.gz"):
            all_pmids.remove(int(abstract_file[:-7]))
    for chunk in __get_chunks(all_pmids, 1000):
        for pmid, abstract in __get_pmids_abstracts(chunk, email):
            abstract_file = os.path.join(abstracts_dir, "%s.txt.gz" % pmid)
            with gzip.open(abstract_file, "wt") as handle:
                if abstract is None:
                    continue
                handle.write(abstract)

    # Abstract to Words
    for abstract_file in os.listdir(abstracts_dir):
        if not abstract_file.endswith(".txt.gz"):
            continue
        words_file = os.path.join(words_dir, "%s.txt.gz" % abstract_file[:-7])
        if not os.path.exists(words_file):
            abstract_file = os.path.join(abstracts_dir, abstract_file)
            with gzip.open(abstract_file, "rt") as handle:
                abstract = "".join([line for line in handle])
            data = [[w, t, c] for w, t, c in __get_abstract_words(abstract)]
            df = pd.DataFrame(data, columns=["Word", "Tag", "Counts"])
            df.to_csv(words_file, sep="\t", index=False, compression="gzip")


def __get_chunks(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks.
    https://docs.python.org/3/library/itertools.html#itertools.zip_longest
    """
    args = [iter(iterable)] * n
    for g in zip_longest(*args, fillvalue=fillvalue):
        yield(list(filter((fillvalue).__ne__, g)))

if __name__ == "__main__":
    wordcloud()