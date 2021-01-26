#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
from functools import partial
import gzip
import hashlib
from itertools import chain, zip_longest
import json
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import sys
from tqdm import tqdm

# Import from lib
sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "utils")
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
    "-i", "--id", "identifiers",
    help="Identifier (e.g. 10664 or P49711).",
    multiple=True
)
@optgroup.option(
    "--input-file",
    help="List of identifiers.",
    type=click.File("rt")
)
@click.option(
    "--add-orthologs",
    help="Include PubMed ID(s) of ortholog(s).",
    is_flag=True
)
@click.option(
    "--background-file",
    help="Background list of identifiers.",
    type=click.File("rt")
)
@click.option(
    "--email",
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
    "--output-dir",
    help="Output directory.",
    default="./",
    show_default=True
)
@click.option(
    "--prefix",
    help="Prefix.  [default: md5 digest]",
)
@click.option(
    "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)

def wordcloud(identifiers, input_file, add_orthologs, background_file, email,
    id_type, output_dir, prefix, threads):

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
        mappings = __get_entrezids_pmids(identifiers, add_orthologs)
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

    # IDFs (i.e. Inverse Document Frequencies)
    # From https://en.wikipedia.org/wiki/Tf-idf
    # IDF = log10(N/nt), where N is total number of documents in the corpus
    # and nt is the number of documents featuring the term t
    tsv_file = os.path.join(output_dir, "idfs.tsv.gz")
    if not os.path.exists(tsv_file):
        data = []; idfs = []
        all_pmids = set(list(chain.from_iterable(pmids+pmids_orthologs)))
        pool = Pool(threads)
        p = partial(__load_PMID_words, words_dir=words_dir)
        for pmid_words in tqdm(pool.imap(p, all_pmids), total=len(all_pmids)):
            for pmid_word in pmid_words:
                data.append(pmid_word)
        df = pd.DataFrame(data, columns=["PMID", "Word"])
        N = df["PMID"].nunique()
        for index, value in df["Word"].value_counts().items():
            idfs.append([index, np.log10(N/float(value))])
        df = pd.DataFrame(idfs, columns=["Word", "IDF"])
        df.to_csv(tsv_file, sep="\t", index=False, compression="gzip")

def __get_chunks(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks.
    https://docs.python.org/3/library/itertools.html#itertools.zip_longest
    """
    args = [iter(iterable)] * n
    for g in zip_longest(*args, fillvalue=fillvalue):
        yield(list(filter((fillvalue).__ne__, g)))

def __load_PMID_words(pmid, words_dir):

    # Initialize
    pmid_words = []
    words_file = os.path.join(words_dir, "%s.txt.gz" % pmid)

    # Load PMID words
    df = pd.read_csv(
        words_file, sep="\t", header=0, usecols=[0, 1]
    )

    return([[pmid, w] for w in set(df["Word"].tolist())])

if __name__ == "__main__":
    wordcloud()