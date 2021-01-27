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
    "--non-alphanumeric",
    help="Filter non-alphanumeric words (e.g. punctuation characters).",
    is_flag=True,
)
@optgroup.option(
    "--non-letters",
    help="Filter non-letter words (e.g. numbers).",
    is_flag=True,
)
@optgroup.option(
    "--stop-words",
    help="Filter stop words (i.e. most common English words).",
    is_flag=True,
)

def abstract2words(abstract, input_file, non_alphanumeric, non_letters, 
    stop_words):

    # Parse abstract
    if abstract is None:
        abstract = ""
        if input_file.name.endswith(".gz"):
            handle = gzip.open(input_file.name, "rt")
            input_file.close()
        else:
            handle = input_file
        for line in handle:
            abstract += line
        handle.close()      

    print("Word\tStem(s)\tCounts")
    for word, stems, counts in __get_abstract_words(abstract, non_alphanumeric,
        non_letters, stop_words):
        print("%s\t%s\t%s" % (word, ";".join(stems), counts))

def __get_abstract_words(abstract, non_alphanumeric, non_letters, stop_words):

    # Initialize
    words = []
    stemmer = SnowballStemmer("english", ignore_stopwords=True)

    # Filter stop words
    if stop_words:
        stopwords = __stop_words()

    # Get words
    for s in [word_tokenize(t) for t in sent_tokenize(abstract.lower())]:
        for token in s:
            if non_alphanumeric:
                if re.search(r"^[^0-9a-z]+$", token):
                    continue
            if non_letters:
                if re.search(r"^[^a-z]+$", token):
                    continue
            if stop_words:
                if __remove_non_alphanumeric(token) in stopwords:
                    continue
            word = __strip_non_alphanumeric(token)
            stems = sorted([stemmer.stem(w) for w in __extract_words(word)])
            words.append((word, tuple(stems)))

    return([[k[0], k[1], v] for k, v in collections.Counter(words).items()])

def __strip_non_alphanumeric(string):
    return(re.sub(r"^[^0-9a-z]+|[^0-9a-z]+$", "", string))

def __extract_words(string):
    return(re.findall(r"\w+", string))

def __remove_non_alphanumeric(string):
    return(re.sub(r"[^0-9a-z]", "", string))

def __stop_words():
    """
    https://github.com/stopwords-iso/stopwords-en
    """

    # Intialize
    url = "https://raw.githubusercontent.com/stopwords-iso/stopwords-en/master/stopwords-en.txt"
    stopwords_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        "stopwords-en.txt")

    # Stop words
    if not os.path.exists(stopwords_file):
        request.urlretrieve(url, stopwords_file)
    with open(stopwords_file) as handle:
        return(set([line.strip("\n") for line in handle]))

if __name__ == "__main__":
    abstract2words()