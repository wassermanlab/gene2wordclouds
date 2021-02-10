#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import collections
import gzip
from nltk import sent_tokenize, word_tokenize
from nltk.stem.snowball import SnowballStemmer
import os
import re
from urllib import request

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@optgroup.group("Input", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "-a", "--abstract",
    help="Abstract (e.g. \"Once upon a time...\")."
)
@optgroup.option(
    "--input-file",
    help="File containing the abstract.",
    type=click.File("rt")
)
@optgroup.group("Word filters")
@optgroup.option(
    "--non-letters",
    help="Filter words without letters.",
    is_flag=True,
)
@optgroup.option(
    "--stop-words",
    help="Filter stop words (i.e. most common English words).",
    is_flag=True,
)

def cli(**params):

    # Parse abstract
    if params["abstract"] is None:
        params["abstract"] = ""
        if params["input_file"].name.endswith(".gz"):
            handle = gzip.open(params["input_file"].name, "rt")
            params["input_file"].close()
        else:
            handle = params["input_file"]
        for line in handle:
            params["abstract"] += line
        handle.close()      

    print("Word\tStem(s)\tCount")
    for word, stems, count in __get_abstract_words(
        params["abstract"], params["non_letters"], params["stop_words"]
    ):
        print("%s\t%s\t%s" % (word, ",".join(stems), count))

def __get_abstract_words(abstract, non_letters=False, stop_words=False):

    # Initialize
    words = []
    stemmer = SnowballStemmer("english", ignore_stopwords=True)
    stopwords = __get_stop_words()

    # Get words
    for s in [word_tokenize(t) for t in sent_tokenize(abstract.lower())]:
        for token in s:
            if non_letters:
                if re.search(r"^[^a-z]+$", token):
                    continue
            if stop_words:
                if __remove_non_alphanumeric_characters(token) in stopwords:
                    continue
            word = __strip_non_alphanumeric_characters(token)
            nr_words = set(__get_alphanumeric_words(word))
            stems = sorted([stemmer.stem(w) for w in nr_words])
            if len(stems) > 0:
                words.append((word, tuple(stems)))

    return([[k[0], k[1], v] for k, v in collections.Counter(words).items()])

def __get_stop_words():
    """
    https://github.com/stopwords-iso/stopwords-enletter_words
    """

    # Intialize
    url = "https://raw.githubusercontent.com/stopwords-iso/stopwords-en/master/stopwords-en.txt"
    stopwords_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        "data", "stopwords-en.txt")

    # Stop words
    if not os.path.exists(stopwords_file):
        request.urlretrieve(url, stopwords_file)
    with open(stopwords_file) as handle:
        return(set([line.strip("\n") for line in handle]))

def __remove_non_alphanumeric_characters(string):
    return(re.sub(r"[^0-9a-z]", "", string))

def __strip_non_alphanumeric_characters(string):
    return(re.sub(r"^[^0-9a-z]+|[^0-9a-z]+$", "", string))

def __get_letter_words(string):
    return(re.findall(r"[a-z]+", string))

def __get_alphanumeric_words(string):
    return(re.findall(r"[0-9a-z]+", string))

if __name__ == "__main__":
    cli()