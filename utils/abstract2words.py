#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import collections
import gzip
import nltk
import os
import re

# Compile
global regexp
regexp = re.compile("[^a-zA-Z]")

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
    help="File with the abstract.",
    type=click.File("rt")
)

def abstract2words(abstract, input_file):

    # Parse abstract
    if abstract is None:
        abstract = ""
        if input_file.name.endswith(".gz"):
            with gzip.open(input_file.name, "rt") as handle:
                for line in handle:
                    abstract += line
        else:
            for line in input_file:
                abstract += line
        input_file.close()

    print("Word\tTag\tCounts")
    for word, tag, counts in __get_abstract_words(abstract):
        print("%s\t%s\t%s" % (word, tag, counts))

def __get_abstract_words(abstract):

    # Get words
    words = __filter_words(nltk.word_tokenize(abstract.lower()))
    words = collections.Counter(words)

    return([[k[0], k[1], v] for k, v in words.items()])

def __filter_words(words):

    # Intialize
    # allowed_tags = {"ADJ", "ADV", "NOUN", "VERB"}
    filtered_words = []
    stopwords_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "stopwords.txt")

    # Filter stop words
    if not os.path.exists(stopwords_file):
        request.urlretrieve(
            "http://ir.dcs.gla.ac.uk/resources/linguistic_utils/stop_words",
            stopwords_file)
    with open(stopwords_file) as handle:
        stopwords = set([line.strip("\n") for line in handle])
    for word, tag in nltk.pos_tag(words, tagset="universal"):
        # if tag not in allowed_tags:
        #     continue
        clean_word = regexp.sub("", word)
        if clean_word.lower() in stopwords:
            continue
        if len(clean_word) == 0:
            continue
        filtered_words.append((word, tag))

    return(filtered_words)

if __name__ == "__main__":
    abstract2words()