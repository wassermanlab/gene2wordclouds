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
from gene2pmid_distribution import gene2pmid_distribution
from pmid2abstract import __get_pmids_abstracts
from uniacc2entrezid import __get_uniaccs_entrezids
from words2cloud import __make_word_cloud

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
    tfidfs_dir = os.path.join(output_dir, "tf-idfs")
    if not os.path.isdir(tfidfs_dir):
        os.makedirs(tfidfs_dir)
    figs_dir = os.path.join(output_dir, "figs")
    if not os.path.isdir(figs_dir):
        os.makedirs(figs_dir)

    # Identifiers
    if input_file is not None:
        identifiers = []
        if input_file.name.endswith(".gz"):
            handle = gzip.open(input_file.name, "rt")
            input_file.close()
        else:
            handle = input_file
        for line in handle:
            if id_type == "uniacc":
                identifiers.append(line.strip("\n"))
            else:
                identifiers.append(int(line.strip("\n")))
        handle.close()

    # Get MD5
    if prefix is None:
        h = hashlib.md5()
        h.update(",".join(identifiers).encode("utf-8"))
        prefix = h.hexdigest()

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

    # Word to IDF (i.e. Inverse Document Frequency)
    # From https://en.wikipedia.org/wiki/Tf-idf recommendations
    #  * IDF = log10(N/n(t))
    # Where N is the total number of documents and nt is the number of
    # documents that include term t.
    tsv_file = os.path.join(output_dir, "word2idf.tsv.gz")
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
    df = pd.read_csv(tsv_file, sep="\t", header=0)
    idfs = dict(zip(df["Word"].tolist(), df["IDF"].tolist()))

    # Word to TF-IDF (i.e. Term Frequency–IDF)
    # From https://en.wikipedia.org/wiki/Tf-idf recommendations
    #  * TF = log10(1+f(t,d)) and TF-IDF = TF * IDF
    # Where f(t,d) is the frequency of term t in document d. However, we
    # focus on the set of documents assigned to gene g (i.e. ng). We have
    # modified the computation of TF accordingly:
    #  * TF = log10(1+ng(t))
    # Where ng(t) is the number of documents of assigned to gene g that
    # include term t.
    if id_type == "uniacc":
        iterator = zip(entrezids, uniaccs, pmids, pmids_orthologs)
    else:
        iterator = zip(entrezids, entrezids, pmids, pmids_orthologs)
    for iteration in iterator:
        if iteration[1] != "P49711":
            continue
        __compute_TFIDFs(iteration, idfs, tfidfs_dir, words_dir)

    # Word Cloud
    for f in os.listdir(tfidfs_dir):
        df = pd.read_csv(os.path.join(tfidfs_dir, f), sep="\t", header=0)
        fig_file = os.path.join(figs_dir, "%s.png" % f[:-7])
        __make_word_cloud(
            df["Word"].tolist(), df["TF-IDF"].tolist(), fig_file
        )   

def __get_chunks(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks.
    https://docs.python.org/3/library/itertools.html#itertools.zip_longest
    """
    args = [iter(iterable)] * n
    for g in zip_longest(*args, fillvalue=fillvalue):
        yield(list(filter((fillvalue).__ne__, g)))

def __load_PMID_words(pmid, words_dir, unique=True):

    # Initialize
    pmid_words = []
    words_file = os.path.join(words_dir, "%s.txt.gz" % pmid)

    # Load PMID words
    df = pd.read_csv(words_file, sep="\t", header=0, usecols=[0, 1])

    if unique:
        return([[pmid, w] for w in set(df["Word"].tolist())])
    else:
        return([[pmid, w] for w in df["Word"].tolist()])

def __compute_TFIDFs(zipped_values, idfs, tfidfs_dir, words_dir):

    # Initialize
    data = []; tfidfs = []
    entrezid = zipped_values[0]
    gene = str(zipped_values[1])
    pmids = zipped_values[2] + zipped_values[3]

    # Compute TF-IDF
    tsv_file = os.path.join(tfidfs_dir, "%s.tsv.gz" % gene)
    if not os.path.exists(tsv_file):
        for pmid in pmids:
            pmid_words = __load_PMID_words(pmid, words_dir, False)
            for pmid_word in pmid_words:
                data.append(pmid_word)
        df = pd.DataFrame(data, columns=["PMID", "Word"])
        N = df["PMID"].nunique()
        df["n(t,g)"] = 1 # set the number of documents that term t occurs to 1
        df = df.groupby("Word").sum().reset_index()
        for idx, row in df.iterrows():
            word = row["Word"]
            ntg = row["n(t,g)"]
            tfidfs.append([word, np.log10(1+float(ntg))*idfs[word]])
        # for index, value in df["Word"].value_counts().items():
        #     tfidfs.append([index, np.log10(1+float(value)/N)*idfs[index]])
        tfidfs.sort(key=lambda x: x[-1], reverse=True)
        df = pd.DataFrame(tfidfs, columns=["Word", "TF-IDF"])
        df.to_csv(tsv_file, sep="\t", index=False, compression="gzip")

if __name__ == "__main__":
    wordcloud()