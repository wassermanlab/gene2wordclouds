#!/usr/bin/env python

import click
from click_option_group import (optgroup, MutuallyExclusiveOptionGroup,
    RequiredMutuallyExclusiveOptionGroup)
import copy
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

# Add utils to path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
    "utils"))

from abstract2words import __get_abstract_words
from entrezid2pmids import __load_datasets, __get_entrezid_pmids
from gene2pmid_distribution import gene2pmid_distribution
from pmid2abstract import __get_pmids_abstracts
from uniacc2entrezid import __get_uniaccs_entrezids
from words2cloud import __get_word_cloud

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
@optgroup.group("Background", cls=MutuallyExclusiveOptionGroup)
@optgroup.option(
    "--background-file",
    help="Background list of identifiers.",
    type=click.File("rt")
)
@optgroup.option(
    "--taxid",
    help="Taxonomic identifier (e.g. 9606).",
    type=int
)
@click.option(
    "-e", "--email",
    help="E-mail address.",
    required=True
)
@click.option(
    "--input-type",
    help="Input type.",
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
    "-p", "--prefix",
    help="Prefix.  [default: md5 digest]",
)
@click.option(
    "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)

def main(**params):

    # Identifiers
    identifiers = __get_identifiers(params["identifiers"],
        params["input_file"], params["input_type"])

    # Directories
    prefix = __get_prefix(identifiers, params["prefix"])
    out_dir = __create_directories(prefix, params["output_dir"])

    # UniAcc to EntrezID
    if params["input_type"] == "uniacc":
        uniaccs, identifiers = __uniacc2entrezid(identifiers, out_dir)  

    # EntrezID to PMIDs
    entrezids, pmids, pmids_orthologs = __entrezid2pmids(identifiers, out_dir)

    ##################
    # Backgroud here #
    ##################

    # PMID to Abstract
    pmids_set = set(list(chain.from_iterable(pmids+pmids_orthologs)))
    __pmid2abstract(copy.copy(pmids_set), params["email"], out_dir)

    # Abstract to Words
    __abstract2words(out_dir, params["threads"])

    # Word to IDF (i.e. Inverse Document Frequency)
    # word2idf, stem2idf = __compute_IDFs(copy.copy(pmids_set), out_dir,
    idfs = __compute_IDFs(copy.copy(pmids_set), out_dir, params["threads"])

    # Word to TF-IDF (i.e. Term Frequency–IDF)
    iterator = __get_iterator(entrezids, pmids, pmids_orthologs, out_dir)
    __compute_TFIDFs(iterator, idfs, out_dir, params["threads"])

    # Word Cloud
    for f in os.listdir(os.path.join(out_dir, "tf-idfs")):
        with gzip.open(os.path.join(out_dir, "tf-idfs", f), "rt") as handle:
            words_tfidfs = json.load(handle)
        words = [word for word, _ in words_tfidfs]
        weights = [tfidf for _, tfidf in words_tfidfs]
        fig_file = os.path.join(out_dir, "figs", "%s.png" % f[:-8])
        __get_word_cloud(words, weights, fig_file)   

def __get_identifiers(identifiers, input_file, input_type):

    # Get identifiers
    if input_file is not None:
        identifiers = []
        if input_file.name.endswith(".gz"):
            handle = gzip.open(input_file.name, "rt")
            input_file.close()
        else:
            handle = input_file
        for line in handle:
            identifiers.append(line.strip("\n"))
        handle.close()
    if input_type == "entrezid":
        identifiers = list(map(int, identifiers))

    return(identifiers)

def __get_prefix(identifiers, prefix=None):

    # Get prefix
    if prefix is None:
        h = hashlib.md5()
        h.update(",".join(identifiers).encode("utf-8"))
        prefix = h.hexdigest()

    return(prefix)

def __create_directories(prefix, output_dir="./"):

    # Create directories
    output_dir = os.path.join(output_dir, prefix)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # if not os.path.isdir(os.path.join(output_dir, "abstracts")):
    #     os.makedirs(os.path.join(output_dir, "abstracts"))
    # if not os.path.isdir(os.path.join(output_dir, "words")):
    #     os.makedirs(os.path.join(output_dir, "words"))
    # if not os.path.isdir(os.path.join(output_dir, "tf-idfs")):
    #     os.makedirs(os.path.join(output_dir, "tf-idfs"))
    # if not os.path.isdir(os.path.join(output_dir, "figs")):
    #     os.makedirs(os.path.join(output_dir, "figs"))

    return(output_dir)

def __uniacc2entrezid(uniaccs, output_dir="./"):

    # Get Entrez Gene IDs
    json_file = os.path.join(output_dir, "uniacc2entrezid.json.gz")
    if not os.path.exists(json_file):
        uniaccs_entrezids = __get_uniaccs_entrezids(uniaccs)
        with gzip.open(json_file, "wt") as handle:
            handle.write(json.dumps(uniaccs_entrezids, indent=4))
    with gzip.open(json_file, "rt") as handle:
        uniaccs_entrezids = json.load(handle)
    uniaccs, entrezids = list(map(list, zip(*uniaccs_entrezids)))

    return(uniaccs, entrezids)

def __entrezid2pmids(entrezids, output_dir="./"):

    # Get PubMed IDs 
    json_file = os.path.join(output_dir, "entrezid2pmids.json.gz")
    if not os.path.exists(json_file):
        entrezids_pmids = []
        kwargs = {"total": len(entrezids), "ncols": 100}
        gene2pubmed, homologene = __load_datasets(True)
        for entrezid in tqdm(entrezids, **kwargs):
            pmids, orthologs_pmids = __get_entrezid_pmids(entrezid,
                gene2pubmed, homologene)
            entrezids_pmids.append([entrezid, pmids, orthologs_pmids])
        with gzip.open(json_file, "wt") as handle:
            handle.write(json.dumps(entrezids_pmids, indent=4))
    with gzip.open(json_file, "rt") as handle:
        entrezids_pmids = json.load(handle)    
    entrezids, pmids, pmids_orthologs = list(map(list, zip(*entrezids_pmids)))

    return(entrezids, pmids, pmids_orthologs)

def __pmid2abstract(pmids, email, output_dir="./"):

    # Initialize
    abstracts_dir = os.path.join(output_dir, "abstracts")

    # Get abstracts
    for abstract_file in os.listdir(abstracts_dir):
        if abstract_file.endswith(".txt.gz"):
            pmids.remove(int(abstract_file[:-7]))
    if len(pmids) > 0:
        chunks = [chunk for chunk in __get_chunks(pmids, 1000)]
        kwargs = {"total": len(chunks), "ncols": 100}
        for chunk in tqdm(chunks, **kwargs):
            for p, a in __get_pmids_abstracts(chunk, email):
                abstract_file = os.path.join(abstracts_dir, "%s.txt.gz" % p)
                with gzip.open(abstract_file, "wt") as handle:
                    if a is not None:
                        handle.write(a)

def __abstract2words(output_dir="./", threads=1):

    # Initialize
    abstract_files = []
    abstracts_dir = os.path.join(output_dir, "abstracts")
    words_dir = os.path.join(output_dir, "words")

    # Get words
    for abstract_file in os.listdir(abstracts_dir):
        if not abstract_file.endswith(".txt.gz"):
            continue
        json_file = os.path.join(words_dir, "%s.json.gz" % abstract_file[:-7])
        if not os.path.exists(json_file):
            abstract_file = os.path.join(abstracts_dir, abstract_file)
            abstract_files.append(abstract_file)
    if len(abstract_files) > 0:
        pool = Pool(threads)
        kwargs = {"total": len(abstract_files), "ncols": 100}
        for abstract_file, abstract_words in tqdm(
            pool.imap(__get_abstract_words_wrapper, abstract_files), **kwargs
        ):
            json_file = os.path.join(words_dir, "%s.json.gz" % \
                os.path.basename(abstract_file)[:-7])
            with gzip.open(json_file, "wt") as handle:
                handle.write(json.dumps(abstract_words, sort_keys=True,
                    indent=4))

def __get_abstract_words_wrapper(abstract_file):

    # Initialize
    abstract_words = []

    # Get words
    with gzip.open(abstract_file, "rt") as handle:
        abstract = "".join([line for line in handle])
    for word, stems, count in __get_abstract_words(abstract, True, True):
        abstract_words.append([word, stems, count])

    return(abstract_file, abstract_words)

def __get_chunks(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks.
    https://docs.python.org/3/library/itertools.html#itertools.zip_longest
    """
    args = [iter(iterable)] * n
    for g in zip_longest(*args, fillvalue=fillvalue):
        yield(list(filter((fillvalue).__ne__, g)))

def __compute_IDFs(pmids, output_dir="./", threads=1):
    """
    From https://en.wikipedia.org/wiki/Tf-idf recommendations:
     * IDF = log10(N/n(t))
    Where N is the total number of documents and nt is the number of
    documents that include term t.
    """

    # Initialize
    # pmids_words = []
    # words_idfs = []
    stem2idf = {}
    pmids_stems = []
    stem2pmids = {}
    stems_idfs = []
    words_dir = os.path.join(output_dir, "words")

    # # Compute word IDFs
    # tsv_file = os.path.join(output_dir, "word2idf.tsv.gz")
    # if not os.path.exists(tsv_file):
    #     pool = Pool(threads)
    #     p = partial(__load_PMID_words, words_dir=words_dir)
    #     for pmid_words in tqdm(pool.imap(p, pmids), **kwargs):
    #         for pmid_word in pmid_words:
    #             pmids_words.append(pmid_word)
    #     df = pd.DataFrame(pmids_words, columns=["PMID", "Word"])
    #     N = df["PMID"].nunique()
    #     for t, nt in df["Word"].value_counts().items():
    #         words_idfs.append([t, np.log10(N/float(nt))])
    #     df = pd.DataFrame(words_idfs, columns=["Word", "IDF"])
    #     df.to_csv(tsv_file, sep="\t", index=False, compression="gzip")
    # df = pd.read_csv(tsv_file, sep="\t", header=0)
    # word2idf = dict(zip(df["Word"].tolist(), df["IDF"].tolist()))

    # Compute stem IDFs
    json_file = os.path.join(output_dir, "stem2idf.json.gz")
    if not os.path.exists(json_file):
        pool = Pool(threads)
        p = partial(__load_PMID_words_and_stems, words_dir=words_dir)
        kwargs = {"total": len(pmids), "ncols": 100}
        for l in tqdm(pool.imap(p, pmids), **kwargs):
            for pmid, _, stems in l:
                pmids_stems.append([pmid, stems])
        unique_stems = set([stems for _, stems in pmids_stems])
        for pmid, stems in pmids_stems:
            for stem in stems:
                stem2pmids.setdefault(stem, set())
                stem2pmids[stem].add(pmid)
        N = len(set([pmid for pmid, _ in pmids_stems]))
        kwargs = {"total": len(unique_stems), "ncols": 100}
        for stems in tqdm(unique_stems, **kwargs):
            nt = __get_stem_nt(stems, stem2pmids)
            stems_idfs.append([list(stems), np.log10(N/float(nt))])
        with gzip.open(json_file, "wt") as handle:
            handle.write(json.dumps(stems_idfs, indent=4))
    with gzip.open(json_file, "rt") as handle:
        stems_idfs = json.load(handle)
    for stems, idf in stems_idfs:
        stem2idf.setdefault(tuple(stems), idf)

    # return(word2idf, stem2idf)
    return(stem2idf)

def __load_PMID_words_and_stems(pmid, words_dir):

    # Initialize
    json_file = os.path.join(words_dir, "%s.json.gz" % pmid)

    # Load PMID words
    with gzip.open(json_file, "rt") as handle:
        pmid_words = json.load(handle)

    # Get unique stems
    unique_words_stems = set(
        [tuple([word, tuple(stems)]) for word, stems, c in pmid_words]
    )

    return([[pmid, word, stems] for word, stems in unique_words_stems])

def __get_stem_nt(stems, stem2pmids):

    # Initialize
    pmids = set()

    # Get stem PMIDs
    for stem in stems:
        pmids = pmids.union(stem2pmids[stem])

    return(len(pmids))

def __get_iterator(entrezids, pmids, pmids_orthologs, output_dir="./"):

    # Initialize
    iterator = []
    tfidfs_dir = os.path.join(output_dir, "tf-idfs")

    # Get iterator
    for i in range(len(entrezids)):
        json_file = os.path.join(tfidfs_dir, "%s.json.gz" % entrezids[i])
        if os.path.exists(json_file):
            continue
        iterator.append([entrezids[i], set(pmids[i] + pmids_orthologs[i])])

    return(iterator)

def __compute_TFIDFs(iterator, idfs, output_dir="./", threads=1):
    """
    From https://en.wikipedia.org/wiki/Tf-idf recommendations
     * TF = log10(1+f(t,d)) and TF-IDF = TF * IDF
    Where f(t,d) is the frequency of term t in document d. However, we
    focus on the set of documents assigned to gene g (i.e. ng). We have
    modified the computation of TF accordingly:
     * TF = log10(1+ng(t))
    Where ng(t) is the number of documents of assigned to gene g that
    include term t.
    """

    # Initialize
    kwargs = {"total": len(iterator), "ncols": 100}
    examples = set([3725, 7157, 6908, 10664])

    for iteration in iterator:
        if iteration[0] in examples:
            __compute_gene_TFIDF(iteration, idfs, output_dir, threads)

    # # Get TF-IDFs
    # if len(iterator) > 0:
    #     pool = Pool(threads)
    #     p = partial(__compute_gene_TFIDF, word2idf=word2idf, stem2idf=stem2idf,
    #         output_dir=output_dir)
    #     for _ in tqdm(pool.imap(p, iterator), **kwargs):
    #         pass

    # # Initialize
    # data = []; tfidfs = []
    # entrezid = zipped_values[0]
    # gene = str(zipped_values[1])
    # words_dir = os.path.join(output_dir, "words")
    # tfidfs_dir = os.path.join(output_dir, "tf-idfs")

    # if id_type == "uniacc":
    #     pass
    # for iteration in iterator:
    #     if iteration[1] != "P49711":
    #         continue
    #     __compute_TFIDFs(iteration, idfs, tfidfs_dir, words_dir)

def __compute_gene_TFIDF(iteration, idfs, output_dir="./", threads=1):

    # Initialize
    done = set()
    stem2pmids = {}
    word2pmids = {}
    word2stems = {}
    pmids_words_stems = []
    words_tfidfs = []
    entrezid = iteration[0]
    pmids = list(iteration[1])
    words_dir = os.path.join(output_dir, "words")
    tfidfs_dir = os.path.join(output_dir, "tf-idfs")

    # Compute word TF-IDFs
    json_file = os.path.join(tfidfs_dir, "%s.json.gz" % entrezid)
    if not os.path.exists(json_file):
        pool = Pool(threads)
        p = partial(__load_PMID_words_and_stems, words_dir=words_dir)
        for l in pool.imap(p, pmids):
            for pmid_word_stems in l:
                pmids_words_stems.append(pmid_word_stems)
        for pmid, word, stems in pmids_words_stems:
            for stem in stems:
                stem2pmids.setdefault(stem, set())
                stem2pmids[stem].add(pmid)
            word2pmids.setdefault(word, set())
            word2pmids[word].add(pmid)
            word2stems.setdefault(word, stems)
        for word in sorted(
            word2pmids, key=lambda x: len(word2pmids[x]), reverse=True
        ):
            if done.intersection(set([s for s in word2stems[word]])):
                continue
            ngt = __get_stem_nt(word2stems[word], stem2pmids)
            idf = idfs[word2stems[word]]
            words_tfidfs.append([word, np.log10(1+ngt)*idf])
            done = done.union(set([s for s in word2stems[word]]))
        with gzip.open(json_file, "wt") as handle:
            handle.write(json.dumps(words_tfidfs, indent=4))

if __name__ == "__main__":
    main()